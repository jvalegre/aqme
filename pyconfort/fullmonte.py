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
import math
import random
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, rdMolAlign, Lipinski,rdchem
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

def rotate_dihedral(mol_rot,dih_rot,args,conf,nsteps):
	rad_range = np.arange(0.0, 360.0,args.ang_fullmonte)
	for i,_ in enumerate(dih_rot):
		random.seed(nsteps)
		rad_ang = random.choice(rad_range)
		rad = math.pi*rad_ang/180.0
		rdMolTransforms.SetDihedralRad(mol_rot.GetConformer(conf),*dih_rot[i],value=rad)
	return mol_rot


def generating_conformations_fullmonte(name,args,rotmatches,log,selectedcids_rdkit,outmols,sdwriter,dup_data,dup_data_idx,coord_Map,alg_Map, mol_template):

	##wroking with fullmonte
	log.write("\n\no  Generation of confomers using FULLMONTE using {0} unique conformer(s) as starting point(s)".format(len(selectedcids_rdkit)))

	#Writing the conformers as mol objects to sdf
	sdtemp = Chem.SDWriter(name+'_'+'rdkit'+args.output)
	for conf in selectedcids_rdkit:
		sdtemp.write(outmols[conf],conf)
	sdtemp.close()

	fmmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
	if fmmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)


	# array for each each unique from rdkit
	unique_mol,c_energy,unique_mol_sample  = [],[],[]

	#STEP 1: Use start conformation for and append to unique list
	nsteps = 1
	for mol_fm in fmmols:
		unique_mol.append(mol_fm)
		c_energy.append(float(mol_fm.GetProp("Energy")))

	#defining unique mol sample for choosing
	globmin = min(c_energy)
	for ene in reversed(c_energy):
		if abs(globmin-ene) < args.ewin_sample_fullmonte:
			unique_mol_sample.append(unique_mol[c_energy.index(ene)])

	# bar = IncrementalBar('o  Generating conformations for Full Monte', max = args.nsteps_fullmonte)
	while nsteps < args.nsteps_fullmonte+1:

		#STEP 2: Choose mol object form unique_mol:
		random.seed(nsteps)
		mol_rot = random.choices(unique_mol_sample,k=1)[0]

		#updating the location of mol object i.e., the hexadecimal locaiton to a new one so the older one isnt affected
		mol = Chem.RWMol(mol_rot)
		rot_mol = mol.GetMol()

		#STEP 3: Choose randmon subset of dihedral from rotmatches
		if len(rotmatches) > args.nrot_fullmonte:
			random.seed(nsteps)
			dih_rot = random.choices(rotmatches, k=args.nrot_fullmonte) #k can be varied base on user definitions
		else:
			random.seed(nsteps)
			dih_rot = random.choices(rotmatches, k=len(rotmatches)) #k can be varied base on user definitions
		#STEP 4: for the given conformation, then apply a random rotation to each torsion in the subset
		rot_mol = rotate_dihedral(rot_mol,dih_rot,args,-1,nsteps)

		#STEP 4: Optimize geometry rot_mol
		if coord_Map is None and alg_Map is None and mol_template is None:
			GetFF = minimize_rdkit_energy(rot_mol,-1,args,log)
			energy = float(GetFF.CalcEnergy())
		else:
			mol,GetFF = realign_mol(rot_mol,-1,coord_Map, alg_Map, mol_template,args,log)
			energy = float(GetFF.CalcEnergy())

		#STEP 5 : Check for DUPLICATES - energy and rms filter (reuse)
				 #  if the conformer is unique then save it the list
		exclude_conf= False
		# compare against allprevious conformers located
		for j,seenmol in enumerate(unique_mol):
			if abs(energy - c_energy[j]) < args.initial_energy_threshold:
				exclude_conf = True
				break
			if  abs(energy - c_energy[j]) < args.energy_threshold:
				rms = get_conf_RMS(rot_mol,seenmol,-1,-1, args.heavyonly, args.max_matches_RMSD,log)
				if rms < args.rms_threshold:
					exclude_conf = True
					break
		if not exclude_conf:
			unique_mol.append(rot_mol)
			c_energy.append(energy)
			unique_mol[c_energy.index(energy)].SetProp("Energy", str(energy))

		unique_mol_sample = []
		#STEP 6: ANALYSE THE UNIQUE list for lowest energy, reorder the uniques if greater the given thershold remove
		globmin = min(c_energy)
		for ene in reversed(c_energy):
			indx = c_energy.index(ene)
			if abs(globmin-ene) > args.ewin_fullmonte:
				unique_mol.pop(indx)
				c_energy.pop(indx)
			if abs(globmin-ene) < args.ewin_sample_fullmonte:

				unique_mol_sample.append(unique_mol[indx])

		nsteps += 1
	# 	bar.next()
	#
	# bar.finish()

	dup_data.at[dup_data_idx, 'FullMonte-Unique-conformers'] = len(unique_mol)

	if args.verbose:
		log.write("o  "+ str(len(unique_mol))+" unique conformers remain")

	cids = list(range(len(unique_mol)))
	sorted_all_cids = sorted(cids,key = lambda cid: c_energy[cid])

	#STEP 9: WRITE FINAL uniques to sdf for xtb or ani
	for i, cid in enumerate(sorted_all_cids):
		unique_mol[cid].SetProp('_Name',name+' '+str(i))
		if coord_Map is None and alg_Map is None and mol_template is None:
			if args.metal_complex:
				set_metal_atomic_number(unique_mol[cid],args)
			sdwriter.write(unique_mol[cid])
		else:
			mol_realigned,_ = realign_mol(unique_mol[cid],-1,coord_Map, alg_Map, mol_template,args,log)
			if args.metal_complex:
				set_metal_atomic_number(mol_realigned,args)
			sdwriter.write(mol_realigned)

	status = 1

	return status
