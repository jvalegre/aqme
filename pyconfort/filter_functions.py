#!/usr/bin/env python

#####################################################.
# 		  This file stores all the functions 	    #
# 	             used for filtering			        #
#####################################################.

import glob
import sys
from progress.bar import IncrementalBar
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, Descriptors
from pyconfort.argument_parser import possible_atoms


possible_atoms = possible_atoms()

def set_metal_atomic_number(mol,args):
	for atom in mol.GetAtoms():
		if atom.GetIdx() in args.metal_idx:
			re_symbol = args.metal_sym[args.metal_idx.index(atom.GetIdx())]
			atomic_number = possible_atoms.index(re_symbol)
			atom.SetAtomicNum(atomic_number)

# RULES TO GET EXPERIMENTAL CONFORMERS
def exp_rules_output(mol, args,log):
	passing = True
	ligand_links = []
	atom_indexes = []
	for atom in mol.GetAtoms():
		# Finds the Ir atom and gets the atom types and indexes of all its neighbours
		if atom.GetSymbol() in args.metal:
			atomic_number = possible_atoms.index(atom.GetSymbol())
			atom.SetAtomicNum(atomic_number)
	for atom in mol.GetAtoms():
		if atom.GetAtomicNum() == atomic_number:
			metal_idx = atom.GetIdx()
			for x in atom.GetNeighbors():
				ligand_links.append(x.GetSymbol())
				atom_indexes.append(x.GetIdx())
	# I need to get the only 3D conformer generated in that mol object for rdMolTransforms
	mol_conf = mol.GetConformer(0)
	# This part will identify the pairs of C and N atoms that are part of the same Ph_Py ligand.
	# The shape of the atom pairs is '[[C1_ATOM_NUMBER, N1_ATOM_NUMBER],[C2, N2],...]'.
	# This information is required for the subsequent filtering process based on angles
	if len(atom_indexes) == args.complex_coord:
		ligand_atoms = []

		for i,_ in enumerate(atom_indexes):
			# This is a filter that excludes molecules that fell apart during DFT geometry
			# optimization (i.e. a N atom from one of the ligands separated from Ir). The
			# max distance allowed can be tuned in length_filter
			bond_length = rdMolTransforms.GetBondLength(mol_conf,metal_idx,atom_indexes[i])
			if ligand_links[i] == 'P':
				length_filter = 2.60
			else:
				length_filter = 2.25
			if bond_length > length_filter:
				passing = False
				break
			for j,_ in enumerate(atom_indexes):
				# Avoid combinations of the same atom with itself
				if atom_indexes[i] != atom_indexes[j]:
					# We know that the ligands never have 2 carbon atoms bonding the Ir atom. We
					# only use atom_indexes[i] for C atoms, and atom_indexes[j] for the potential
					# N atoms that are part of the same Ph_Py ligand
					if ligand_links[i] == 'C':
						# This part detects the Ir-C bond and breaks it, breaking the Ph_Py ring
						bond = mol.GetBondBetweenAtoms(atom_indexes[i], metal_idx)
						new_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()],addDummies=True, dummyLabels=[(atom_indexes[i], metal_idx)])
						if new_mol.GetAtomWithIdx(atom_indexes[i]).IsInRingSize(5):
							five_mem = True
						else:
							five_mem = False
						# identify whether or not the initial 5-membered ring formed between [-Ir-C-C-C-N-] is broken when we break the Ir-C bond. This works
						# because Ph_Py units bind Ir in the same way always, through 1 C and 1 N that are in the same position, forming a 5-membered ring.
						# If this ring is broken, atom_indexes[j] will not be part of a 5-membered ring (atom.IsInRingSize(5) == False) which means that
						# this atom was initially inside the same ligand as the parent C of atom_indexes[i])
						if not five_mem:
							if not new_mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
								bond_2 = mol.GetBondBetweenAtoms(atom_indexes[j], metal_idx)
								new_mol_2 = Chem.FragmentOnBonds(mol, [bond_2.GetIdx()],addDummies=True, dummyLabels=[(atom_indexes[j], metal_idx)])
								#doing backwards as well eg. Ir N bond
								if not new_mol_2.GetAtomWithIdx(atom_indexes[i]).IsInRingSize(5):
									ligand_atoms.append([atom_indexes[i],atom_indexes[j]])
									break
						else:
							if not new_mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
								ligand_atoms.append([atom_indexes[i],atom_indexes[j]])
								break
		if passing:
			# This stop variable and the breaks inside the inner loops will break the nested loop if there
			# is one angle that does not meet the criteria for valid conformers
			stop = False
			# For complexes with 3 Ph_Py ligands:
			if len(ligand_atoms) == 3:
				for i,_ in enumerate(ligand_atoms):
					if not stop:
						for j,_ in enumerate(ligand_atoms):
							# the i<=j part avoids repeating atoms, the i != j part avoid angles
							# containing the same number twice (i.e. 4-16-4, this angle will fail)
							if i <= j and i != j:
								# Calculate the angle between 2 N atoms from different Ph_Py ligands.
								# When there are 3 Ph_Py ligands, no 2 N atoms must be in 180 degrees
								angle = rdMolTransforms.GetAngleDeg(mol_conf,ligand_atoms[i][1],metal_idx,ligand_atoms[j][1])
								if (180 - args.angle_off) <= angle <= (180 + args.angle_off):
									passing = False
									break
			# For complexes with 2 Ph_Py ligands + 1 ligand that is not Ph_Py
			if len(ligand_atoms) == 2:
				# Since there are only 2 N atoms, we do not need to include a nested loop
					angle = rdMolTransforms.GetAngleDeg(mol_conf,ligand_atoms[0][1],metal_idx,ligand_atoms[1][1])
					# Calculate the angle between 2 N atoms from different Ph_Py ligands.
					# When there are 2 Ph_Py ligands, the 2 N atoms from the 2 Ph_Py ligands must be in 180 degrees
					if (180 - args.angle_off) <= angle <= (180 + args.angle_off):
						pass
					else:
						passing = False
	# it filters off molecules that the SDF only detects 5 Ir neighbours
	else:
		passing = False
	return passing

# MAIN OPTION FOR DISCARDING MOLECULES BASED ON USER INPUT DATA (REFERRED AS EXPERIMENTAL RULES)
def exp_rules_main(args,log):
	if args.verbose:
		log.write("   ----- Applying experimental rules to write the new confs file -----")
	### do 2 cases, for RDKit only and RDKIt+xTB
	#grab all the gaussian files
	if not args.xtb and not args.ANI1ccx:
		conf_files =  glob.glob('*_rdkit.sdf')
	elif args.xtb:
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.ANI1ccx:
		conf_files =  glob.glob('*_ani.sdf')

	for file in conf_files:
		allmols = Chem.SDMolSupplier(file, removeHs=False)
		if allmols is None:
			log.write("Could not open "+ file)
			sys.exit(-1)

		sdwriter = Chem.SDWriter(file.split('.')[0]+args.exp_rules_output_ext)

		for mol in allmols:
			check_mol = True
			check_mol = exp_rules_output(mol,args,log)
			if check_mol:
				sdwriter.write(mol)
		sdwriter.close()

# FILTER TO BE APPLIED FOR SMILES
def filters(mol,args,log):
	valid_structure = True
	# Second filter: molecular weight
	if Descriptors.MolWt(mol) < args.max_MolWt:
		# Third filter: this filters salts off (2 separated components)
		#if len(Chem.MolToSmiles(mol).split('.')) == 1:
		for atom in mol.GetAtoms():
			#Fourth filter: atoms outside the scope chosen in 'possible_atoms'
			if atom.GetSymbol() not in possible_atoms:
				valid_structure = False
				if args.verbose:
					log.write(" Exiting as atom isn't in atoms in the periodic table")
	else:
		valid_structure = False
		if args.verbose:
			log.write(" Exiting as total molar mass > {0}".format(args.max_MolWt))
	return valid_structure

# CALCULATES RMSD between two molecules
def get_conf_RMS(mol1, mol2, c1, c2, heavy, max_matches_RMSD,log):
	if heavy:
		 mol1 = Chem.RemoveHs(mol1)
		 mol2 = Chem.RemoveHs(mol2)
	rms = Chem.GetBestRMS(mol1,mol2,c1,c2,maxMatches=max_matches_RMSD)
	return rms

def filter_after_rotation(args,name,log,dup_data,dup_data_idx):
	rdmols = Chem.SDMolSupplier(name+'_'+'rdkit'+args.output, removeHs=False)
	if rdmols is None:
		log.write("Could not open "+ name+args.output)
		sys.exit(-1)

	writer_mol_objects = []
	bar = IncrementalBar('o  Filtering based on energy and rms after rotation of dihedrals', max = len(rdmols))

	rd_count = 0
	rd_selectedcids,rd_dup_energy,rd_dup_rms_eng =[],-1,0
	for i, rd_mol_i in enumerate(rdmols):
		mol_rd = Chem.RWMol(rd_mol_i)
		mol_rd.SetProp('_Name',rd_mol_i.GetProp('_Name')+' '+str(i))
		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process
		if rd_count == 0:
			rd_selectedcids.append(rd_mol_i)
			if args.metal_complex:
				set_metal_atomic_number(mol_rd,args)
			writer_mol_objects.append([mol_rd,float(mol_rd.GetProp('Energy'))])
		# Only the first ID gets included
		rd_count = 1
		# check rmsd
		for rd_mol_j in rd_selectedcids:
			if abs(float(rd_mol_i.GetProp('Energy')) - float(rd_mol_j.GetProp('Energy'))) < args.initial_energy_threshold: # comparison in kcal/mol
				excluded_conf = True
				rd_dup_energy += 1
				break
			if abs(float(rd_mol_i.GetProp('Energy')) - float(rd_mol_j.GetProp('Energy'))) < args.energy_threshold: # in kcal/mol
				rms = get_conf_RMS(mol_rd,rd_mol_j,-1,-1, args.heavyonly, args.max_matches_RMSD,log)
				if rms < args.rms_threshold:
					excluded_conf = True
					rd_dup_rms_eng += 1
					break
		if not excluded_conf:
			if rd_mol_i not in rd_selectedcids:
				rd_selectedcids.append(rd_mol_i)
				if args.metal_complex:
					set_metal_atomic_number(mol_rd,args)
				writer_mol_objects.append([mol_rd,float(mol_rd.GetProp('Energy'))] )
		bar.next()
	bar.finish()

	# writing sorted mol objects
	sdwriter_rd = Chem.SDWriter(name+'_'+'rdkit'+'_'+'rotated'+args.output)
	sortedmols = sorted(writer_mol_objects,key=lambda x: x[1])
	for i, write_mol in enumerate(sortedmols):
		sdwriter_rd.write(write_mol[0])

	sdwriter_rd.close()

	if args.verbose:
		log.write("o  "+str(rd_dup_energy)+ " Duplicates removed initial energy (E < "+str(args.initial_energy_threshold)+" kcal/mol)")
	if args.verbose:
		log.write("o  "+str(rd_dup_rms_eng)+ " Duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol) after rotation")
	if args.verbose:
		log.write("o  "+str(len(rd_selectedcids) )+ " unique conformers remain")

	# filtering process after rotations
	dup_data.at[dup_data_idx, 'RDKIT-Rotated-Unique-conformers'] = len(rd_selectedcids)

	status = 1

	return status
