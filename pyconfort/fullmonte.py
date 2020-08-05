#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	  used for genrating details for fullmonte      #
#####################################################.

import os
import sys
import subprocess
import time
import numpy as np
import random
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski
from rdkit.Geometry import Point3D
from progress.bar import IncrementalBar

from pyconfort.argument_parser import possible_atoms
from pyconfort.filter import ewin_filter,pre_E_filter,RMSD_and_E_filter,set_metal_atomic_number
from pyconfort.filter import get_conf_RMS


possible_atoms = possible_atoms()

#steps to realign mol
def realign_mol(mol,conf,coord_Map, alg_Map, mol_template,args,log):
	num_atom_match = mol.GetSubstructMatch(mol_template)
	GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
	for k, idxI in enumerate(num_atom_match):
		for l in range(k + 1, len(num_atom_match)):
			idxJ = num_atom_match[l]
			d = coord_Map[idxI].Distance(coord_Map[idxJ])
			GetFF.AddDistanceConstraint(idxI, idxJ, d, d, 10000)
	GetFF.Initialize()
	GetFF.Minimize(maxIts=args.opt_steps_RDKit)
	# rotate the embedded conformation onto the core_mol:
	rdMolAlign.AlignMol(mol, mol_template, prbCid=conf,refCid=-1,atomMap=alg_Map,reflect=True,maxIters=100)
	return mol,GetFF

#function to get rdkit energy
def minimize_rdkit_energy(mol,conf,args,log):
	if args.ff == "MMFF":
		GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
	elif args.ff == "UFF":
		GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
	else:
		log.write(' Force field {} not supported!'.format(args.ff))
		sys.exit()
	try:
		GetFF.Initialize()
		GetFF.Minimize(maxIts=args.opt_steps_RDKit)
	except:
		#if not MMFF use UFF
		GetFF = Chem.UFFGetMoleculeForceField(mol,confId=conf)
		GetFF.Initialize()
		GetFF.Minimize(maxIts=args.opt_steps_RDKit)
	return GetFF

def rotate_dihedral(mol_rot,dih_rot,args,conf):
	for i,_ in enumerate(dih_rot):
		rdMolTransforms.SetDihedralRad(mol_rot.GetConformer(conf),*dih_rot[i],value=args.ang_fullmonte)
	return mol_rot


def generating_conformations_fullmonte(name,args,rotmatches,log,selectedcids_rdkit,outmols,sdwriter,dup_data,dup_data_idx,coord_Map,alg_Map, mol_template):

	##wroking with fullmonte
	log.write("\n\no  Generation of confomers using FULLMONTE using {0} uniques conformer as starting points".format(len(selectedcids_rdkit)))

	#Writing the conformers as mol objects to sdf
	sdtemp = Chem.SDWriter(name+'_'+'rdkit'+args.output)
	for conf in selectedcids_rdkit:
		sdtemp.write(outmols[conf],conf)
	sdtemp.close()

	fmmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
	if fmmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	#array of all unique mols after fukllmonte for each unique from rdkit
	all_unique_mol,all_unique_mol_fin,all_unique_ene_fin = [],[],[]
	# array for each each unique from rdkit

	bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(fmmols))
	# no=1
	for mol_fm in fmmols:

		unique_mol,c_energy = [],[]

		#STEP 1: Use start conformation for and append to unique list #DO this uniformly by using a uniform sampling tecnique
		nsteps = 1
		unique_mol.append(mol_fm)
		c_energy.append(float(mol_fm.GetProp("Energy")))


		while nsteps < args.nsteps_fullmonte:
			#STEP 2: Choose mol object form unique_mol:
			mol_rot = random.sample(unique_mol,1)[0]

			#STEP 3: Choose randmon subset of dihedral from rotmatches
			if len(rotmatches) > args.nrot_fullmonte:
				dih_rot = random.sample(rotmatches, k=args.nrot_fullmonte) #k can be varied base on user definitions
			else:
				dih_rot = random.choices(rotmatches, k=len(rotmatches)) #k can be varied base on user definitions
			#STEP 4: for the given conformation, then apply a random rotation to each torsion in the subset
			rot_mol = rotate_dihedral(mol_rot,dih_rot,args,-1)
			#STEP 4: Optimize geometry rot_mol
			GetFF = minimize_rdkit_energy(rot_mol,-1,args,log)
			energy = float(GetFF.CalcEnergy())

			rot_mol.SetProp("Energy", str(energy))

			#STEP 5 : Check for DUPLICATES - energy and rms filter (reuse)
					 #  if the conformer is unique then save it the list
			exclude_conf= False
			# compare against allprevious conformers located
			for j,seenmol in enumerate(unique_mol):
				if abs(energy - c_energy[j]) < args.energy_threshold:
					rms = get_conf_RMS(rot_mol,seenmol,-1,-1, args.heavyonly, args.max_matches_RMSD,log)
					if rms < args.rms_threshold:
						exclude_conf = True
						break

			if not exclude_conf:
				unique_mol.append(rot_mol)
				c_energy.append(energy)

			#STEP 6: ANALYSE THE UNIQUE list for lowest energy, reorder the uniques if greater the given thershold remove
			globmin = min(c_energy)
			for ind,_ in enumerate(c_energy):
				if  globmin - c_energy[ind] > args.ewin_fullmonte:
					c_energy.pop(ind)
					unique_mol.pop(ind)
			nsteps += 1

		# sdtemp = Chem.SDWriter(str(no)+'_'+'temp'+args.output)
		# for mol in unique_mol:
		# 	sdtemp.write(mol)
		# sdtemp.close()
		# no+=1

		#appending to all_uniq Mol
		all_unique_mol.append(unique_mol)
		bar.next()
	bar.finish()

	#STEP 7: from a list of list do all three FILTERS
	for mol_list in all_unique_mol:
		all_unique_mol_fin += mol_list

	dup_data.at[dup_data_idx, 'FullMonte-conformers'] = len(all_unique_mol_fin)

	for i, ro_mol_i in enumerate(all_unique_mol_fin):
		all_unique_ene_fin.append(float(ro_mol_i.GetProp('Energy')))

	fm_cids = list(range(len(all_unique_mol_fin)))
	sorted_fm_cids = sorted(fm_cids, key = lambda cid: all_unique_ene_fin[cid])

	# filter based on energy window ewin_csearch
	sortedcids_fm = ewin_filter(sorted_fm_cids,all_unique_ene_fin,args,dup_data,dup_data_idx,log,'fullmonte')
	# pre-filter based on energy only
	selectedcids_initial_fm = pre_E_filter(sortedcids_fm,all_unique_ene_fin,args,dup_data,dup_data_idx,log,'fullmonte')
	# filter based on energy and RMSD
	selectedcids_fm = RMSD_and_E_filter(all_unique_mol_fin,selectedcids_initial_fm,all_unique_ene_fin,args,dup_data,dup_data_idx,log,'fullmonte')

	#STEP 9: WRITE FINAL uniques to sdf for xtb or ani
	for i, cid in enumerate(selectedcids_fm):
		mol_ro = Chem.RWMol(all_unique_mol_fin[cid])
		if coord_Map is None and alg_Map is None and mol_template is None:
			if args.metal_complex:
				set_metal_atomic_number(mol_ro,args)
			mol_ro.SetProp('_Name',mol_rot.GetProp("_Name")+'_'+str(i))
			sdwriter.write(mol_ro)
		else:
			mol_ro_realigned,_ = realign_mol(mol_ro,-1,coord_Map, alg_Map, mol_template,args,log)
			if args.metal_complex:
				set_metal_atomic_number(mol_ro_realigned,args)
			mol_ro.SetProp('_Name',mol_rot.GetProp("_Name")+'_'+str(i+1))
			sdwriter.write(mol_ro_realigned)

	status = 1
	#removes the rdkit file
	os.remove(name+'_'+'rdkit'+args.output)

	return status
