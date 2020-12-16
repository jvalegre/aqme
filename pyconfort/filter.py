#!/usr/bin/env python

#####################################################.
# 		  This file stores all the functions 	    #
# 	             used for filtering			        #
#####################################################.

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
def exp_rules_output(mol,args,log,file,print_error_exp_rules,ob_compat,rdkit_compat):
	passing = True
	for rule in args.exp_rules:
		if rule == 'Ir_bidentate_x3':
			ligand_links = []
			atom_indexes = []
			# Finds the Ir atom and gets the atom types and indexes of all its neighbours
			# the filter is compatible with molecules that do not contain Ir (always passing)
			metal_idx = None
			for atom in mol.GetAtoms():
				if atom.GetAtomicNum() == 77:
					metal_idx = atom.GetIdx()
					for x in atom.GetNeighbors():
						ligand_links.append(x.GetSymbol())
						atom_indexes.append(x.GetIdx())

			# I need to get the only 3D conformer generated in that mol object for rdMolTransforms
			mol_conf = mol.GetConformer(0)
			# This part will identify the pairs of C and N atoms that are part of the same Ph_Py ligand.
			# The shape of the atom pairs is '[[C1_ATOM_NUMBER, N1_ATOM_NUMBER],[C2, N2],...]'.
			# This information is required for the subsequent filtering process based on angles
			if len(atom_indexes) == 6:
				ligand_atoms = []
				for i,_ in enumerate(atom_indexes):
					# This is a filter that excludes molecules that fell apart during DFT geometry
					# optimization (i.e. a N atom from one of the ligands separated from Ir). The
					# max distance allowed can be tuned in length_filter
					bond_length = rdMolTransforms.GetBondLength(mol_conf,metal_idx,atom_indexes[i])
					length_filter = 2.60 # based on observation from DFT optimized geometries
					if bond_length > length_filter:
						passing = False
						break
					for j,_ in enumerate(atom_indexes):
						# Avoid combinations of the same atom with itself
						if atom_indexes[i] != atom_indexes[j]:
							# We know that the ligands never have 2 carbon atoms bonding the Ir atom except
							# for carbenes. We only use atom_indexes[i] for C atoms, and atom_indexes[j] for the potential
							# N atoms that are part of the same Ph_Py ligand
							if ligand_links[i] == 'C':
								carbene_like = False
								for inside_neighbour in mol.GetAtoms()[atom_indexes[i]].GetNeighbors():
									N_group = ['N', 'P', 'As']
									if inside_neighbour.GetSymbol() in N_group:
										if inside_neighbour.GetTotalValence() == 4:
											carbene_like = True
								if not carbene_like:
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
										bond_2 = mol.GetBondBetweenAtoms(atom_indexes[j], metal_idx)
										new_mol_2 = Chem.FragmentOnBonds(mol, [bond_2.GetIdx()],addDummies=True, dummyLabels=[(atom_indexes[j], metal_idx)])
										#doing backwards as well eg. Ir N bond
										if not new_mol_2.GetAtomWithIdx(atom_indexes[i]).IsInRingSize(5):
											ligand_atoms.append([atom_indexes[i],atom_indexes[j]])
											break

									else:
										if not new_mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
											if mol.GetAtomWithIdx(atom_indexes[j]).IsInRingSize(5):
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
				if metal_idx is not None:
					passing = False

		else:
			try:
				atoms_filter = rule.split(',')[0].split('-')
				angle_rules = rule.split(',')[1]
			except IndexError:
				log.write('x  The exp_rules parameter(s) was not correctly defined, this filter will be turned off')
				break
			# the elements of this initial list will be replaced by the corresponding atom id numebrs
			atom_idx = ['ATOM1','ATOM2','ATOM3']

			find_angle = 0
			incompatibility = False
			for atom in mol.GetAtoms():
				# count matches
				neigh_count_first = 0
				neigh_count_second = 0
				# Finds the metal atom and gets the atom types and indexes of all its neighbours
				if atom.GetSymbol() == atoms_filter[1]:
					# idx of the central atom
					atom_idx[1] = atom.GetIdx()
					for x in atom.GetNeighbors():
						if x.GetSymbol() == atoms_filter[0] or x.GetSymbol() == atoms_filter[2]:
							# this ensures that both neighbours are used
							if x.GetSymbol() == atoms_filter[0] and x.GetSymbol() == atoms_filter[2]:
								if neigh_count_first <= neigh_count_second:
									neigh_count_first += 1
									atom_idx[0] = x.GetIdx()
								else:
									neigh_count_second += 1
									atom_idx[2] = x.GetIdx()
							elif x.GetSymbol() == atoms_filter[0]:
								neigh_count_first += 1
								atom_idx[0] = x.GetIdx()
							elif x.GetSymbol() == atoms_filter[2]:
								neigh_count_second += 1
								atom_idx[2] = x.GetIdx()
							# count matches
							matches = neigh_count_first + neigh_count_second
							if matches > 2:
								if print_error_exp_rules == 0:
									log.write('x  There are multiple options in exp_rules for '+ file + ', this filter will be turned off')
									incompatibility = True
									break
					if neigh_count_first == 1 and neigh_count_second == 1:
						find_angle += 1
			if find_angle == 0 and not incompatibility:
				if print_error_exp_rules == 0:
					log.write('x  No angles matching the description from exp_rules in '+ file + ', this filter will be turned off')
			elif find_angle > 1:
					log.write('x  '+ file + ' contain more than one atom that meets the exp_rules criteria, this filter will be turned off')
			elif find_angle == 1:
				# I need to get the only 3D conformer generated in that mol object for rdMolTransforms
				mol_conf = mol.GetConformer(0)

				# Calculate the angle between the 3 elements
				angle = rdMolTransforms.GetAngleDeg(mol_conf,atom_idx[0],atom_idx[1],atom_idx[2])
				if (int(angle_rules) - args.angle_off) <= angle <= (int(angle_rules) + args.angle_off):
					passing = False
			if not passing:
				break

	return passing

# COMPARES THE BOND LENGTHS BETWEEN AN INPUT AND AN OPTIMIZED OUTPUT FILES
def check_geom_filter(mol,mol2,args):
	passing_geom = True
	mol_bonds = []
	for atom in mol.GetAtoms():
		for x in atom.GetNeighbors():
			indiv_bond = [atom.GetIdx(),x.GetIdx()]
			mol_bonds.append(indiv_bond)

	# I need to get the only 3D conformer generated in that mol object for rdMolTransforms
	mol_conf = mol.GetConformer(0)
	mol_conf2 = mol2.GetConformer(0)
	for bond in mol_bonds:
		bond_length = rdMolTransforms.GetBondLength(mol_conf,bond[0],bond[1])
		try:
			bond_length2 = rdMolTransforms.GetBondLength(mol_conf2,bond[0],bond[1])
		except RuntimeError:
			passing_geom = False
			break

		if bond_length < bond_length2:
			smaller_bond = bond_length
			bigger_bond = bond_length2
		else:
			smaller_bond = bond_length2
			bigger_bond = bond_length

		if bigger_bond > args.length_criteria*smaller_bond:
			passing_geom = False

	return passing_geom

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

# filter based on energy window (ewin_csearch)
def ewin_filter(sorted_all_cids,cenergy,args,dup_data,dup_data_idx,log,calc_type):
	sortedcids,nhigh_csearch,nhigh=[],0,0
	#if calc_type == 'rdkit' or calc_type == 'summ' or calc_type == 'fullmonte':
	if calc_type == 'rdkit' or calc_type == 'summ':
		for i,cid in enumerate(sorted_all_cids):
			if i == 0:
				cenergy_min = cenergy[cid]
			if abs(cenergy[cid] - cenergy_min) < args.ewin_csearch:
				sortedcids.append(cid)
			else:
				nhigh_csearch +=1
	if calc_type == 'rdkit':
		if args.verbose:
			log.write("o  "+str(nhigh_csearch)+ " conformers rejected based on energy window ewin_csearch (E > "+str(args.ewin_csearch)+" kcal/mol)")
		dup_data.at[dup_data_idx, 'RDKit-energy-window'] = nhigh_csearch
	if calc_type == 'summ':
		if args.verbose:
			log.write("o  "+str(nhigh_csearch)+ " conformers rejected after rotation based on energy window ewin_csearch (E > "+str(args.ewin_csearch)+" kcal/mol)")
		dup_data.at[dup_data_idx, 'summ-energy-window'] = nhigh_csearch
	# if calc_type == 'fullmonte':
	# 	if args.verbose:
	# 		log.write("o  "+str(nhigh_csearch)+ " conformers rejected based on energy window ewin_csearch (E > "+str(args.ewin_csearch)+" kcal/mol)")
	# 	dup_data.at[dup_data_idx, 'FullMonte-energy-window'] = nhigh_csearch

	if calc_type == 'xtb' or calc_type == 'ani':
		for i,cid in enumerate(sorted_all_cids):
			if i == 0:
				cenergy_min = cenergy[cid]
			if abs(cenergy[cid] - cenergy_min) < args.ewin_cmin:
				sortedcids.append(cid)
			else:
				nhigh +=1
		if args.verbose:
			log.write("o  "+str(nhigh)+ " conformers rejected based on energy window ewin_cmin (E > "+str(args.ewin_cmin)+" kcal/mol)")
		if calc_type == 'ani':
			dup_data.at[dup_data_idx, 'ANI-energy-window'] = nhigh
		elif calc_type == 'xtb':
			dup_data.at[dup_data_idx, 'xTB-energy-window'] = nhigh
	return sortedcids

# pre-energy filter for RDKit
def pre_E_filter(sortedcids,cenergy,args,dup_data,dup_data_idx,log,calc_type):
	selectedcids_initial, eng_dup =[],-1
	# bar = IncrementalBar('o  Filtering based on energy (pre-filter)', max = len(sortedcids))
	for i, conf in enumerate(sortedcids):
		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process
		if i == 0:
			selectedcids_initial.append(conf)
		# check rmsd
		for seenconf in selectedcids_initial:
			E_diff = abs(cenergy[conf] - cenergy[seenconf]) # in kcal/mol
			if E_diff < args.initial_energy_threshold:
				eng_dup += 1
				excluded_conf = True
				break
		if not excluded_conf:
			if conf not in selectedcids_initial:
				selectedcids_initial.append(conf)
		# bar.next()
	# bar.finish()

	if args.verbose:
		log.write("o  "+str(eng_dup)+ " duplicates removed  pre-energy filter (E < "+str(args.initial_energy_threshold)+" kcal/mol)")

	if calc_type == 'rdkit':
		dup_data.at[dup_data_idx, 'RDKit-initial_energy_threshold'] = eng_dup
	if calc_type == 'summ':
		dup_data.at[dup_data_idx, 'summ-initial_energy_threshold'] = eng_dup
	# if calc_type == 'fullmonte':
	# 	dup_data.at[dup_data_idx, 'FullMonte-initial_energy_threshold'] = eng_dup
	elif calc_type == 'xtb' or calc_type == 'ani':
		if calc_type == 'ani':
			dup_data.at[dup_data_idx, 'ANI-initial_energy_threshold'] = eng_dup
		elif calc_type == 'xtb':
			dup_data.at[dup_data_idx, 'xTB-initial_energy_threshold'] = eng_dup

	return selectedcids_initial

# filter based on energy and RMSD
def RMSD_and_E_filter(outmols,selectedcids_initial,cenergy,args,dup_data,dup_data_idx,log,calc_type):
	if args.verbose:
		log.write("o  Removing duplicate conformers (RMSD < "+ str(args.rms_threshold)+ " and E difference < "+str(args.energy_threshold)+" kcal/mol)")
	# bar = IncrementalBar('o  Filtering based on energy and RMSD', max = len(selectedcids_initial))

	selectedcids,eng_rms_dup = [],-1
	for i, conf in enumerate(selectedcids_initial):
		# This keeps track of whether or not your conformer is unique
		excluded_conf = False
		# include the first conformer in the list to start the filtering process

		if i == 0:
			selectedcids.append(conf)
		# check energy and rmsd
		for seenconf in selectedcids:
			E_diff = abs(cenergy[conf] - cenergy[seenconf]) # in kcal/mol
			if  E_diff < args.energy_threshold:
				if calc_type == 'rdkit':
					rms = get_conf_RMS(outmols[seenconf],outmols[conf],seenconf,conf, args.heavyonly, args.max_matches_RMSD,log)
				#elif calc_type == 'summ' or calc_type == 'fullmonte' or calc_type =='xtb' or calc_type =='ani':
				elif calc_type == 'summ' or calc_type =='xtb' or calc_type =='ani':
					rms = get_conf_RMS(outmols[conf],outmols[seenconf],-1,-1, args.heavyonly, args.max_matches_RMSD,log)
				if rms < args.rms_threshold:
					excluded_conf = True
					eng_rms_dup += 1
					break
		if not excluded_conf:
			if conf not in selectedcids:
				selectedcids.append(conf)
		# bar.next()
	# bar.finish()

	if args.verbose:
		log.write("o  "+str(eng_rms_dup)+ " duplicates removed (RMSD < "+str(args.rms_threshold)+" / E < "+str(args.energy_threshold)+" kcal/mol)")
	if args.verbose:
		log.write("o  "+ str(len(selectedcids))+" unique conformers remain")

	if calc_type == 'rdkit':
		dup_data.at[dup_data_idx, 'RDKit-RMSD-and-energy-duplicates'] = eng_rms_dup
		dup_data.at[dup_data_idx, 'RDKIT-Unique-conformers'] = len(selectedcids)
	# if calc_type == 'fullmonte':
	# 	dup_data.at[dup_data_idx, 'FullMonte-RMSD-and-energy-duplicates'] = eng_rms_dup
	# 	dup_data.at[dup_data_idx, 'FullMonte-Unique-conformers'] = len(selectedcids)
	if calc_type == 'summ':
		dup_data.at[dup_data_idx, 'summ-RMSD-and-energy-duplicates'] = eng_rms_dup
		dup_data.at[dup_data_idx, 'summ-Unique-conformers'] = len(selectedcids)
	elif calc_type == 'ani' or calc_type == 'xtb':
		if calc_type == 'ani':
			dup_data.at[dup_data_idx, 'ANI-RMSD-and-energy-duplicates'] = eng_rms_dup
			dup_data.at[dup_data_idx, 'ANI-Unique-conformers'] = len(selectedcids)
		elif calc_type == 'xtb':
			dup_data.at[dup_data_idx, 'xTB-RMSD-and-energy-duplicates'] = eng_rms_dup
			dup_data.at[dup_data_idx, 'xTB-Unique-conformers'] = len(selectedcids)

	return selectedcids
