#!/usr/bin/python
from __future__ import print_function
import argparse, math, os, sys, traceback
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, PropertyMol
from rdkit.Geometry import Point3D

import DBGEN.db_gen_functions

### TORCHANI IMPORTS
import ase
import ase.optimize
import torch, torchani
os.environ['KMP_DUPLICATE_LIB_OK']='True'
device = torch.device('cpu')
model = torchani.models.ANI1ccx()
from ase.units import kJ,mol,Hartree,kcal

#import xtb
#from xtb import GFN2
output = '.sdf'
final_output = '_confs.sdf'

def get_conf_RMS(mol, c1,c2, heavy):
	'''generate RMS distance between two molecules (ignoring hydrogens)'''
	if heavy == True: mol = Chem.RemoveHs(mol)
	rms = Chem.GetBestRMS(mol,mol,c1,c2)
	return rms

def getPMIDIFF(mol1, mol2):
	pmi1 = Chem.CalcPMI1(mol1) - Chem.CalcPMI1(mol2)
	pmi2 = Chem.CalcPMI2(mol1) - Chem.CalcPMI2(mol2)
	pmi3 = Chem.CalcPMI3(mol1) - Chem.CalcPMI3(mol2)
	diff = (pmi1 **2 + pmi2 **2 + pmi3 **2) ** 0.5
	return diff

def getDihedralMatches(mol, heavy):
	'''return list of atom indices of dihedrals'''
	#this is rdkit's "strict" pattern
	pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
	qmol = Chem.MolFromSmarts(pattern)
	matches = mol.GetSubstructMatches(qmol);

	#these are all sets of 4 atoms, uniquify by middle two
	uniqmatches = []
	seen = set()
	for (a,b,c,d) in matches:
		if (b,c) not in seen:
			if heavy == True:
				if mol.GetAtomWithIdx(a).GetSymbol() != 'H' and mol.GetAtomWithIdx(d).GetSymbol() != 'H':
					seen.add((b,c))
					uniqmatches.append((a,b,c,d))
			if heavy == False:
				seen.add((b,c))
				uniqmatches.append((a,b,c,d))
	return uniqmatches

def genConformer_r(mol, conf, i, matches, degree, sdwriter):
	'''recursively enumerate all angles for rotatable dihedrals.  i is
	which dihedral we are enumerating by degree to output conformers to out'''
	if i >= len(matches): #base case, torsions should be set in conf
		sdwriter.write(mol,conf)
		return 1
	else:
		#incr = math.pi*degree / 180.0
		total = 0
		deg = 0
		while deg < 360.0:
			rad = math.pi*deg / 180.0
			rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
			total += genConformer_r(mol, conf, i+1, matches, degree, sdwriter)
			deg += degree
		return total

def summ_search(mol, name,args):
	'''embeds core conformers, then optimizes and filters based on RMSD. Finally the rotatable torsions are systematically rotated'''

	sdwriter = Chem.SDWriter(name+output)

	Chem.SanitizeMol(mol)
	mol = Chem.AddHs(mol)
	mol.SetProp("_Name",name);

	if args.nodihedrals == False: rotmatches = getDihedralMatches(mol, args.heavyonly)
	else: rotmatches = []

	if len(rotmatches) > args.max_torsions:
		print("x  Too many torsions (%d). Skipping %s" %(len(rotmatches),(name+output)))
	else:
		if args.etkdg:
			cids = Chem.EmbedMultipleConfs(mol, args.sample, Chem.ETKDG(),randomSeed=args.seed)
		else:
			cids = Chem.EmbedMultipleConfs(mol, args.sample,randomSeed=args.seed)
		if args.verbose:
			print("o ", len(cids),"conformers initially sampled")

		#energy minimize all to get more realistic results
		if args.verbose: print("o  Optimizing", len(cids), "initial conformers with", args.ff)
		if args.verbose:
			if args.nodihedrals == False: print("o  Found", len(rotmatches), "rotatable torsions")
			else: print("o  Systematic torsion rotation is set to OFF")

		cenergy = []
		for i, conf in enumerate(cids):

			#identify the atoms and decide Force Field
			for atom in mol.GetAtoms():
				if atom.GetAtomicNum() > 36: #upto Kr for MMFF, if not use UFF
					args.ff = "UFF"

					#print("UFF is used because there are atoms that MMFF doesn't recognise")

			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol),confId=conf)
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol)
			else: print('   Force field {} not supported!'.format(args.ff)); sys.exit()

			GetFF.Initialize()
			converged = GetFF.Minimize()
			cenergy.append(GetFF.CalcEnergy())
			#if args.verbose:
			#    print("-   conformer", (i+1), "optimized: ", args.ff, "energy", GetFF.CalcEnergy())

		#reduce to unique set
		if args.verbose: print("o  Removing duplicate conformers ( RMSD <", args.rms_threshold, ")")
		sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
		selectedcids = []
		for i, conf in enumerate(sortedcids):

			#set torsions to zero
			if len(rotmatches) != 0:
				for m in rotmatches:
					rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*m,value=0)
			#check rmsd
			for seenconf in selectedcids:
				rms = get_conf_RMS(mol,seenconf,conf, args.heavyonly)
				if rms < args.rms_threshold:
					break
			else: #loop completed normally - no break, included empty
				selectedcids.append(conf)

		#now exhaustively drive torsions of selected conformers
		if args.verbose: print("o ", len(selectedcids),"unique conformers remain")
		n_confs = int(len(selectedcids) * (360 / args.degree) ** len(rotmatches))
		if args.verbose and len(rotmatches) != 0: print("o  Systematic generation of", n_confs, "confomers")

		total = 0
		for conf in selectedcids:
			total += genConformer_r(mol, conf, 0, rotmatches, args.degree, sdwriter)
		if args.verbose and len(rotmatches) != 0: print("o  %d total conformations generated"%total)

	sdwriter.close()

def mult_min(mol, name,args):
	'''optimizes a bunch of molecules and then checks for unique conformers and then puts in order of energy'''

	opt = True # switch to off for single point only
	opt_precision = 0.005 # toggle for optimization convergence

	#adjust opt convergence criteria (args.convergence defaults to 1.0)
	opt_precision = opt_precision * args.convergence

	inmols = Chem.SDMolSupplier(name+output, removeHs=False)
	if inmols is None:
		print("Could not open ", name+output)
		sys.exit(-1)

	c_converged, c_energy, outmols = [], [], []
	ani_energy,xtb_energy = 0,0
	if args.ANI1ccx == True or args.xtb == True: SQM_energy, SQM_cartesians = [], []

	globmin = None

	for i,mol in enumerate(inmols):
		conf = 1
		if mol is not None:

			if args.ff == "MMFF":
				GetFF = Chem.MMFFGetMoleculeForceField(mol, Chem.MMFFGetMoleculeProperties(mol))
			elif args.ff == "UFF":
				GetFF = Chem.UFFGetMoleculeForceField(mol)
			else: print(('   Force field {} not supported!'.format(args.ff))); sys.exit()

			GetFF.Initialize()
			converged = GetFF.Minimize(maxIts=1000)
			energy = GetFF.CalcEnergy()
			# append to list
			#if args.verbose: print("   conformer", (i+1), energy)
			if globmin == None: globmin = energy
			if energy < globmin: globmin = energy

			if converged == 0 and (energy - globmin) < args.ewin:
				#if args.verbose: print('   minimization converged!')
				unique, dup_id = 0, None
				#print("Conformer", (i+1), "optimized with", args.ff, "Energy:", energy)
				for j,seenmol in enumerate(outmols):
					if abs(energy - c_energy[j]) < args.energy_threshold:
						#print((i+1), energy, (j+1), c_energy[j], getPMIDIFF(mol,seenmol))
						if getPMIDIFF(mol, seenmol) < args.rms_threshold * 25:
							#print("o  Conformer", (i+1), "matches conformer", (j+1))
							unique += 1
							dup_id = (j+1)

				if unique == 0:
					if args.verbose == True: print("-  Conformer", (i+1), "is unique")

					if args.ANI1ccx == True or args.xtb == True:
						cartesians = mol.GetConformers()[0].GetPositions()
						elements = ''
						for atom in mol.GetAtoms(): elements += atom.GetSymbol()

						coordinates = torch.tensor([cartesians.tolist()], requires_grad=True, device=device)

						if args.ANI1ccx == True:
							species = model.species_to_tensor(elements).to(device).unsqueeze(0)
							_, ani_energy = model((species, coordinates))
							if args.verbose: print("ANI Initial E:",ani_energy.item(),'eH') #Hartree

							if opt == True:
								ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=model.ase())
								### make a function for constraints and optimization
								if constraints != None:
									fb = ase.constraints.FixBondLength(0, 1)
									ase_molecule.set_distance(0,1,2.0)
									ase_molecule.set_constraint(fb)

								optimizer = ase.optimize.BFGS(ase_molecule)
								optimizer.run(fmax=float(opt_precision))
								species_coords = ase_molecule.get_positions().tolist()
								coordinates = torch.tensor([species_coords], requires_grad=True, device=device)

							###############################################################################
							# Now let's compute energy:
							_, ani_energy = model((species, coordinates))
							aniE = ani_energy.item() #Hartree
							if args.verbose: print("ANI Final E:", aniE,'eH', ase_molecule.get_potential_energy(),'eV') #Hartree, eV
							###############################################################################
### INCLUDE THE OPTIONS TO SOTRE MOLECULAR Descriptors
### CHECK THIS WEBPAGE: https://github.com/grimme-lab/xtb/tree/master/python
						elif args.xtb == True:
							ase_molecule = ase.Atoms(elements, positions=coordinates.tolist()[0], calculator=GFN2()) #define ase molecule using GFN2 Calculator
							if opt == True:
								if args.verbose: print("Initial XTB energy", ase_molecule.get_potential_energy()/Hartree,'Eh',ase_molecule.get_potential_energy(),'eV') #Hartree, eV
								optimizer = ase.optimize.BFGS(ase_molecule)
								optimizer.run(fmax=float(opt_precision))
								species_coords = ase_molecule.get_positions().tolist()
								coordinates = torch.tensor([species_coords], requires_grad=True, device=device)
							###############################################################################
							# Now let's compute energy:
							xtb_energy = ase_molecule.get_potential_energy()
							if args.verbose: print("Final XTB E:",xtb_energy/Hartree,'Eh',xtb_energy,'eV') #Hartree, eV
							###############################################################################

						if args.ANI1ccx == True or args.xtb == True:#save Eh and coordinates to write to SDF
							if args.xtb == True:SQM_energy.append(xtb_energy/Hartree)
							else:SQM_energy.append(ani_energy.item())
							cartesians = np.array(coordinates.tolist()[0])
							SQM_cartesians.append(cartesians)


					pmol = PropertyMol.PropertyMol(mol)
					outmols.append(pmol); c_converged.append(converged); c_energy.append(energy)
					conf += 1

				else: print("x  Conformer", (i+1), "is a duplicate of", dup_id)
			else:
				print("x  Minimization of conformer", (i+1), " not converged / energy too high!", converged, (energy - globmin), args.ewin)
			#pass
		else:
			pass #print("No molecules to optimize")


	# if SQM energy exists, overwrite RDKIT energies and geometries
	cids = list(range(len(outmols)))
	sortedcids = sorted(cids, key = lambda cid: c_energy[cid])

	if args.ANI1ccx == True or args.xtb == True:
		for conf in cids:
			c_energy[conf] = SQM_energy[conf]
			c = outmols[conf].GetConformer()
			for j in range(outmols[conf].GetNumAtoms()):
				#print(cartesians[i])
				[x,y,z] = SQM_cartesians[conf][j]
				c.SetAtomPosition(j,Point3D(x,y,z))

			for j in range(0,conf):
				if abs(c_energy[conf] - c_energy[j]) < args.energy_threshold / 2625.5 and getPMIDIFF(outmols[conf], outmols[j]) <  args.rms_threshold:
					print("It appears ",conf, "is the same as", j)

	for i, cid in enumerate(sortedcids):
		outmols[cid].SetProp('_Name', name + ' conformer ' + str(i+1))
		outmols[cid].SetProp('Energy', c_energy[cid])

	return outmols, c_energy
