#####################################################.
#            This file stores all the functions     #
#               used in conformer generation        #
#####################################################.

import math
import os
import sys
import subprocess
import time
from pathlib import Path
import concurrent.futures as futures  # RAUL: This is for the main

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem

from rdkit.Chem import rdMolTransforms, PropertyMol, rdDistGeom, Lipinski
from progress.bar import IncrementalBar  # RAUL: This is for the main

from aqme.filter import filters, ewin_filter, pre_E_filter, RMSD_and_E_filter
from aqme.tmbuild import template_embed
from aqme.cmin import rules_get_charge, substituted_mol
from aqme.fullmonte import (
	generating_conformations_fullmonte,
	minimize_rdkit_energy,
	realign_mol,
)
from aqme.utils import (
	Logger,
	set_metal_atomic_number,
	com_2_xyz_2_sdf,
	getDihedralMatches,
	load_from_yaml,
	creation_of_dup_csv_csearch,
	check_charge_smi,
	check_for_pieces,
	nci_ts_mol,
)
from aqme.crest import crest_opt

from aqme.argument_parser import set_options

SUPPORTED_INPUTS = [
	".smi",
	".sdf",
	".cdx",
	".csv",
	".com",
	".gjf",
	".mol",
	".mol2",
	".xyz",
	".txt",
	".yaml",
	".yml",
	".rtf",
]

# get pre-determined geometries for metal complexes
accepted_complex_types = ["squareplanar", "squarepyramidal", "linear", "trigonalplanar"]


def process_csearch(smi_, name_, dir_, varfile_, charge_default_, constraints_dist_, constraints_angle_, constraints_dihedral_):
	obj = csearch(smi_, name_, dir_, varfile_, charge_default_, constraints_dist_, constraints_angle_, constraints_dihedral_)
	total_data = obj.compute_confs()
	return total_data


class csearch:
	"""
	Representation of the neccesary information related with csearch.

	Parameters
	----------
	mol : RDKit Mol object, neededSMILES string necessary for setting up csearch object
	name : Name of smiles, neededstring representing the code name for the object.
	"""

	def __init__(
		self,
		smi=None,
		name=None,
		w_dir_initial=os.getcwd(),
		varfile=None,
		charge_defualt=0,
		constraints_dist=[],
		constraints_angle=[],
		constraints_dihedral=[],
		**kwargs,
	):
		self.smi = smi
		self.name = name
		self.w_dir_initial = w_dir_initial
		self.charge_default = charge_defualt

		if "options" in kwargs:
			self.args = kwargs["options"]
		else:
			self.args = set_options(kwargs)
		self.args.varfile = varfile

		csearch_dir = Path(self.w_dir_initial) / "CSEARCH"
		dat_dir = csearch_dir / "dat_files"
		dat_dir.mkdir(parents=True, exist_ok=True)

		self.log = Logger(dat_dir / self.name, self.args.output_name)

		if varfile is not None:
			self.args, self.log = load_from_yaml(self.args, self.log)

		self.mol , self.args.constraints_dist, self.args.constraints_angle, self.args.constraints_dihedral  = smi_to_mol(self.smi, self.args, self.name, constraints_dist, constraints_angle, constraints_dihedral)

		self.args.charge_default = self.charge_default
		self.args.charge = rules_get_charge(self.mol, self.args)

		self.csearch_folder = Path(self.w_dir_initial).joinpath(f"CSEARCH/{self.args.CSEARCH}")
		self.csearch_folder.mkdir(exist_ok=True)
		self.csearch_file = self.csearch_folder.joinpath(self.name + "_" + self.args.CSEARCH + self.args.output)
		self.sdwriter = Chem.SDWriter(str(self.csearch_file))

	def compute_confs(self):

		"""
		function to start conf generation

		Parameters
		----------
		w_dir_initial : [type]    [description]
		mol : rdkit.Chem.Mol    [description]
		name : [type]    [description]
		i : [type]    [description]

		Returns
		-------
		pandas.Dataframe    total_data
		"""

		# Converts each line to a rdkit mol object
		if self.args.verbose:
			self.log.write(f"   -> Input Molecule {Chem.MolToSmiles(self.mol)}")

		if self.args.metal_complex:
			for _ in self.args.metal:
				self.args.metal_idx.append(None)
				self.args.complex_coord.append(None)
				self.args.metal_sym.append(None)

			(
				self.args.mol,
				self.args.metal_idx,
				self.args.complex_coord,
				self.args.metal_sym,
			) = substituted_mol(self.mol, self.args)

			# get pre-determined geometries for metal complexes
			accepted_complex_types = [
				"squareplanar",
				"squarepyramidal",
				"linear",
				"trigonalplanar",
			]
			if self.args.complex_type in accepted_complex_types:
				count_metals = 0
				for metal_idx_ind in self.args.metal_idx:
					if metal_idx_ind is not None:
						count_metals += 1
				if count_metals == 1:
					os.chdir(self.w_dir_initial)
					template_kwargs = dict()
					template_kwargs["complex_type"] = self.args.complex_type
					template_kwargs["metal_idx"] = self.args.metal_idx
					template_kwargs["maxsteps"] = self.args.opt_steps_RDKit
					template_kwargs["heavyonly"] = self.args.heavyonly
					template_kwargs["maxmatches"] = self.args.max_matches_RMSD
					items = template_embed(self.mol, self.name, self.log, **template_kwargs)

					total_data = creation_of_dup_csv_csearch(self.args.CSEARCH)

					for mol_obj, name_in, coord_map, alg_map, template in zip(*items):
						data = self.conformer_generation(
							mol_obj,
							name_in,
							self.args,
							self.log,
							coord_map,
							alg_map,
							template,
						)
						frames = [total_data, data]
						total_data = pd.concat(frames, sort=True)
				else:
					log.write(
						"x  Cannot use templates for complexes involving more than 1 metal or for organic molecueles."
					)
					total_data = None
			else:
				total_data = self.conformer_generation(
					self.mol, self.name, self.args, self.log
				)
		else:
			total_data = self.conformer_generation(
				self.mol, self.name, self.args, self.log
			)
		return total_data

	def conformer_generation(
		self, mol, name, args, log, coord_Map=None, alg_Map=None, mol_template=None
	):

		"""
		FUNCTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS

		Parameters
		----------
		mol : rdkit.Chem.Mol    [description]
		name : [type]    [description]
		args : argparse.args    [description]
		log : [type]    [description]
		coord_Map : [type], optional    [description], by default None
		alg_Map : [type], optional    [description], by default None
		mol_template : [type], optional    [description], by default None

		Returns
		-------
		pd.Dataframe    dup_data
		"""
		dup_data = creation_of_dup_csv_csearch(args.CSEARCH)

		dup_data_idx = 0
		start_time = time.time()
		valid_structure = filters(mol, log, args.max_MolWt, args.verbose)
		if valid_structure:
			if args.verbose:
				log.write(f"\n   ----- {name} -----")
			try:
				# the conformational search for RDKit
				status, update_to_rdkit = self.summ_search(
					mol,
					name,
					args,
					log,
					dup_data,
					dup_data_idx,
					coord_Map,
					alg_Map,
					mol_template,
				)
				dup_data.at[dup_data_idx, "status"] = status
				dup_data.at[dup_data_idx, "update_to_rdkit"] = update_to_rdkit
			except (KeyboardInterrupt, SystemExit):  # RAUL: This try-except is useless.
				raise
		else:
			log.write("\nx  ERROR: The structure is not valid")

		if args.time:
			n_seconds = round(time.time() - start_time, 2)
			log.write(f"\n Execution time CSEARCH: {n_seconds} seconds")
			dup_data.at[dup_data_idx, "CSEARCH time (seconds)"] = n_seconds
		return dup_data

	def summ_search(
		self,
		mol,
		name,
		args,
		log,
		dup_data,
		dup_data_idx,
		coord_Map=None,
		alg_Map=None,
		mol_template=None,
	):

		"""
		EMBEDS, OPTIMIZES AND FILTERS RDKIT CONFORMERS

		Parameters
		----------
		mol : [type]    [description]
		name : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]
		dup_data : [type]    [description]
		dup_data_idx : [type]    [description]
		coord_Map : [type], optional    [description], by default None
		alg_Map : [type], optional    [description], by default None
		mol_template : [type], optional    [description], by default None

		Returns
		-------
		tuple    status, update_to_rdkit
		"""

		# writes sdf for the first RDKit conformer generation
		status, rotmatches, update_to_rdkit = self.rdkit_to_sdf(
			mol,
			name,
			args,
			log,
			dup_data,
			dup_data_idx,
			self.sdwriter,
			coord_Map,
			alg_Map,
			mol_template

		)
		# reads the initial SDF files from RDKit and uses dihedral scan if selected
		if status != -1 or status != 0:
			# getting the energy and mols after rotations
			if args.CSEARCH == "summ" and len(rotmatches) != 0:
				status = self.dihedral_filter_and_sdf(
					name,
					args,
					log,
					dup_data,
					dup_data_idx,
					coord_Map,
					alg_Map,
					mol_template,
				)

				# removes the rdkit file
				os.remove(name + "_" + "rdkit" + args.output)

		return status, update_to_rdkit

	def dihedral_filter_and_sdf(
		self, name, args, log, dup_data, dup_data_idx, coord_Map, alg_Map, mol_template
	):
		"""
		filtering after dihedral scan to sdf

		Parameters
		----------
		name : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]
		dup_data : [type]    [description]
		dup_data_idx : [type]    [description]
		coord_Map : [type]    [description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		int    status (job I guess?)
		"""
		rotated_energy = []

		rdmols = Chem.SDMolSupplier(name + "_" + "rdkit" + args.output, removeHs=False)

		if rdmols is None:
			log.write("Could not open " + name + args.output)
			sys.exit(-1)

		for i, rd_mol_i in enumerate(rdmols):
			if coord_Map is None and alg_Map is None and mol_template is None:
				energy = minimize_rdkit_energy(
					rd_mol_i, -1, log, args.ff, args.opt_steps_RDKit
				)
			else:
				rd_mol_i, energy = realign_mol(
					rd_mol_i, -1, coord_Map, alg_Map, mol_template, args.opt_steps_RDKit
				)
			rotated_energy.append(energy)

		rotated_cids = list(range(len(rdmols)))
		sorted_rotated_cids = sorted(rotated_cids, key=lambda cid: rotated_energy[cid])

		# filter based on energy window ewin_csearch
		sortedcids_rotated = ewin_filter(
			sorted_rotated_cids,
			rotated_energy,
			args,
			dup_data,
			dup_data_idx,
			log,
			"summ",
			args.ewin_csearch,
		)
		# pre-filter based on energy only
		selectedcids_initial_rotated = pre_E_filter(
			sortedcids_rotated,
			rotated_energy,
			dup_data,
			dup_data_idx,
			log,
			"summ",
			args.initial_energy_threshold,
			args.verbose,
		)
		# filter based on energy and RMSD
		selectedcids_rotated = RMSD_and_E_filter(
			rdmols,
			selectedcids_initial_rotated,
			rotated_energy,
			args,
			dup_data,
			dup_data_idx,
			log,
			"summ",
		)

		sdwriter_rd = Chem.SDWriter(name + "_" + "summ" + args.output)
		for i, cid in enumerate(selectedcids_rotated):
			mol_rd = Chem.RWMol(rdmols[cid])
			mol_rd.SetProp("_Name", rdmols[cid].GetProp("_Name") + " " + str(i))
			mol_rd.SetProp("Energy", str(rotated_energy[cid]))
			if args.metal_complex:
				set_metal_atomic_number(mol_rd, args.metal_idx, args.metal_sym)
			sdwriter_rd.write(mol_rd)
		sdwriter_rd.close()
		status = 1
		return status

	def clean_args(
		self, args, ori_ff, mol, ori_charge
	):  # RAUL: I hope this function does not survive the clean-ups - Shree : not needed anymore (just leaving it in here for now)
		"""
		returns the arguments to their original value after each calculation

		Parameters
		----------
		args : argparse.args    [description]
		ori_ff : [type]    original forcefield
		mol : rdkit.Chem.Mol    [description]
		ori_charge : [type]    original charge
		"""
		for atom in mol.GetAtoms():
			if atom.GetSymbol() in args.metal:
				args.metal_complex = True
				break
		else:
			args.metal_complex = False
		args.ff = ori_ff
		args.charge_default = ori_charge
		args.metal_idx = []
		args.complex_coord = []
		args.metal_sym = []

	def auto_sampling(self, mult_factor, mol, args, log):
		"""
		DETECTS INITIAL NUMBER OF SAMPLES AUTOMATICALLY

		Parameters
		----------
		mult_factor : [type]    [description]
		mol : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]

		Returns
		-------
		int    auto_samples
		"""
		if args.metal_complex:
			if len(args.metal_idx) > 0:
				mult_factor = (
					mult_factor * 3 * len(args.metal_idx)
				)  # this accounts for possible trans/cis isomers in metal complexes
		auto_samples = 0
		auto_samples += 3 * (Lipinski.NumRotatableBonds(mol))  # x3, for C3 rotations
		auto_samples += 3 * (Lipinski.NHOHCount(mol))  # x3, for OH/NH rotations
		auto_samples += 3 * (
			Lipinski.NumSaturatedRings(mol)
		)  # x3, for boat/chair/envelope confs
		if auto_samples == 0:
			auto_samples = mult_factor
		else:
			auto_samples = mult_factor * auto_samples
		return auto_samples

	def genConformer_r(
		self,
		mol,
		conf,
		i,
		matches,
		degree,
		sdwriter,
		args,
		name,
		log,
		update_to_rdkit,
		coord_Map,
		alg_Map,
		mol_template,
	):
		"""
		IF NOT USING DIHEDRALS, THIS REPLACES I BACK TO THE METAL WHEN METAL = TRUE
		AND WRITES THE RDKIT SDF FILES. WITH DIHEDRALS, IT OPTIMIZES THE ROTAMERS

		Parameters
		----------
		mol : rdkit.Chem.Mol    [description]
		conf : [type]    [description]
		i : [type]    [description]
		matches : [type]    [description]
		degree : [type]    [description]
		sdwriter : [type]    [description]
		args : [type]    [description]
		name : [type]    [description]
		log : [type]    [description]
		update_to_rdkit : [type]    [description]
		coord_Map : [type]    [description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		int    total number of conformers generated
		"""
		if i >= len(matches):  # base case, torsions should be set in conf
			# setting the metal back instead of I
			if args.metal_complex and (args.CSEARCH == "rdkit" or update_to_rdkit):
				if coord_Map is None and alg_Map is None and mol_template is None:
					energy = minimize_rdkit_energy(
						mol, conf, log, args.ff, args.opt_steps_RDKit
					)
				else:
					mol, energy = realign_mol(
						mol,
						conf,
						coord_Map,
						alg_Map,
						mol_template,
						args.opt_steps_RDKit,
					)
				mol.SetProp("Energy", str(energy))
				set_metal_atomic_number(mol, args.metal_idx, args.metal_sym)
			sdwriter.write(mol, conf)
			return 1
		else:
			total = 0
			deg = 0
			while deg < 360.0:
				rad = math.pi * deg / 180.0
				rdMolTransforms.SetDihedralRad(
					mol.GetConformer(conf), *matches[i], value=rad
				)
				mol.SetProp("_Name", name)
				total += self.genConformer_r(
					mol,
					conf,
					i + 1,
					matches,
					degree,
					sdwriter,
					args,
					name,
					log,
					update_to_rdkit,
					coord_Map,
					alg_Map,
					mol_template,
				)
				deg += degree
			return total

	def embed_conf(
		self, mol, initial_confs, args, log, coord_Map, alg_Map, mol_template
	):
		"""
		function to embed conformers

		Parameters
		----------
		mol : [type]    [description]
		initial_confs : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]
		coord_Map : [type]    [description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		list    cids
		"""
		is_sdf_mol_or_mol2 = os.path.splitext(args.input)[1] in [
			".sdf",
			".mol",
			".mol2",
		]

		if is_sdf_mol_or_mol2:
			Chem.AssignStereochemistryFrom3D(mol)

		embed_kwargs = dict()
		embed_kwargs["ignoreSmoothingFailures"] = True
		embed_kwargs["randomSeed"] = args.seed
		embed_kwargs["numThreads"] = 0

		if (coord_Map, alg_Map, mol_template) != (None, None, None):
			embed_kwargs["coordMap"] = coord_Map
		cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

		if len(cids) <= 1 and initial_confs != 1:
			log.write(
				"o  Normal RDKit embeding process failed, trying to "
				"generate conformers with random coordinates "
				f"(with {str(initial_confs)} possibilities)"
			)
			embed_kwargs["useRandomCoords"] = True
			embed_kwargs["boxSizeMult"] = 10.0
			embed_kwargs["numZeroFail"] = 1000
			embed_kwargs["numThreads"] = 1
			cids = rdDistGeom.EmbedMultipleConfs(mol, initial_confs, **embed_kwargs)

		if is_sdf_mol_or_mol2:
			# preserving AssignStereochemistryFrom3D
			for cid in cids:
				Chem.AssignAtomChiralTagsFromStructure(mol, confId=cid)

		return cids

	def min_and_E_calc(self, mol, cids, args, log, coord_Map, alg_Map, mol_template):
		"""
		minimization and E calculation with RDKit after embeding

		Parameters
		----------
		mol : [type] [description]
		cids : [type][description]
		args : [type][description]
		log : [type][description]
		coord_Map : [type][description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		outmols,cenergy    outmols is gonna be a list containing "initial_confs" mol objects    with "initial_confs" conformers. We do this to SetProp    (Name and Energy) to the different conformers    and log.write in the SDF file. At the end, since all the mol    objects has the same conformers, but the energies are different,    we can log.write conformers to SDF files with the energies of the    parent mol objects. We measured the computing time and it's the    same as using only 1 parent mol object with 10 conformers, but we    couldn'temp SetProp correctly.
		"""

		cenergy, outmols = [], []
		# bar = IncrementalBar('o  Minimizing', max = len(cids))
		for _, conf in enumerate(cids):
			if coord_Map is None and alg_Map is None and mol_template is None:
				energy = minimize_rdkit_energy(
					mol, conf, log, args.ff, args.opt_steps_RDKit
				)
			else:  # id template realign before doing calculations
				mol, energy = realign_mol(
					mol, conf, coord_Map, alg_Map, mol_template, args.opt_steps_RDKit
				)
			cenergy.append(energy)
			pmol = PropertyMol.PropertyMol(mol)
			outmols.append(pmol)
			# bar.next()
		# bar.finish()
		return outmols, cenergy

	def min_after_embed(
		self,
		mol,
		cids,
		name,
		initial_confs,
		rotmatches,
		dup_data,
		dup_data_idx,
		sdwriter,
		args,
		log,
		update_to_rdkit,
		coord_Map,
		alg_Map,
		mol_template,
	):
		"""
		minimizes, gets the energy and filters RDKit conformers after embeding

		Parameters
		----------
		mol : [type]    [description]
		cids : [type]    [description]
		name : [type]    [description]
		initial_confs : [type]    [description]
		rotmatches : [type]    [description]
		dup_data : [type]    [description]
		dup_data_idx : [type]    [description]
		sdwriter : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]
		update_to_rdkit : [type]    [description]
		coord_Map : [type]    [description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		int    status
		"""

		# gets optimized mol objects and energies
		outmols, cenergy = self.min_and_E_calc(
			mol, cids, args, log, coord_Map, alg_Map, mol_template
		)

		# writing charges after RDKit
		if (
			os.path.splitext(args.input)[1] == ".cdx"
			or os.path.splitext(args.input)[1] == ".smi"
			or os.path.splitext(args.input)[1] == ".csv"
		):
			dup_data.at[dup_data_idx, "Overall charge"] = np.sum(args.charge)
		else:
			dup_data.at[dup_data_idx, "Overall charge"] = args.charge_default

		for i, cid in enumerate(cids):
			outmols[cid].SetProp("_Name", name + " " + str(i + 1))
			outmols[cid].SetProp("Energy", str(cenergy[cid]))
			outmols[cid].SetProp(
				"Real charge", str(dup_data.at[dup_data_idx, "Overall charge"])
			)

		# sorts the energies
		cids = list(range(len(outmols)))
		sorted_all_cids = sorted(cids, key=lambda cid: cenergy[cid])

		log.write("\n\no  Applying filters to intial conformers")

		# filter based on energy window ewin_csearch
		sortedcids_rdkit = ewin_filter(
			sorted_all_cids,
			cenergy,
			args,
			dup_data,
			dup_data_idx,
			log,
			"rdkit",
			args.ewin_csearch,
		)

		# pre-filter based on energy only
		selectedcids_initial_rdkit = pre_E_filter(
			sortedcids_rdkit,
			cenergy,
			dup_data,
			dup_data_idx,
			log,
			"rdkit",
			args.initial_energy_threshold,
			args.verbose,
		)

		# filter based on energy and RMSD
		selectedcids_rdkit = RMSD_and_E_filter(
			outmols,
			selectedcids_initial_rdkit,
			cenergy,
			args,
			dup_data,
			dup_data_idx,
			log,
			"rdkit",
		)

		if args.CSEARCH == "summ" or args.CSEARCH == "rdkit":
			# now exhaustively drive torsions of selected conformers
			n_confs = int(
				len(selectedcids_rdkit) * (360 / args.degree) ** len(rotmatches)
			)
			if args.verbose and len(rotmatches) != 0:
				log.write(
					"\n\no  Systematic generation of " + str(n_confs) + " confomers"
				)
				# bar = IncrementalBar('o  Generating conformations based on dihedral rotation', max = len(selectedcids_rdkit))
			# else:
			#     bar = IncrementalBar('o  Writing unique conformers into an sdf file', max = len(selectedcids_rdkit))

			total = 0
			for conf in selectedcids_rdkit:
				if args.CSEARCH == "summ" and not update_to_rdkit:
					sdwriter.write(outmols[conf], conf)
					for m in rotmatches:
						rdMolTransforms.SetDihedralDeg(
							outmols[conf].GetConformer(conf), *m, 180.0
						)
				total += self.genConformer_r(
					outmols[conf],
					conf,
					0,
					rotmatches,
					args.degree,
					sdwriter,
					args,
					outmols[conf].GetProp("_Name"),
					log,
					update_to_rdkit,
					coord_Map,
					alg_Map,
					mol_template,
				)
				# bar.next()
			# bar.finish()
			if args.verbose and len(rotmatches) != 0:
				log.write("o  %d total conformations generated" % total)
			status = 1

		if args.CSEARCH == "summ":
			dup_data.at[dup_data_idx, "summ-conformers"] = total

		if args.CSEARCH == "fullmonte":
			status = generating_conformations_fullmonte(
				name,
				args,
				rotmatches,
				log,
				selectedcids_rdkit,
				outmols,
				sdwriter,
				dup_data,
				dup_data_idx,
				coord_Map,
				alg_Map,
				mol_template,
			)
			# removes the rdkit file
			os.remove(name + "_" + "rdkit" + args.output)

		return status

	def rdkit_to_sdf(
		self,
		mol,
		name,
		args,
		log,
		dup_data,
		dup_data_idx,
		sdwriter,
		coord_Map,
		alg_Map,
		mol_template
	):

		"""
		conversion from rdkit to sdf

		Parameters
		----------
		mol : [type]    [description]
		name : [type]    [description]
		args : [type]    [description]
		log : [type]    [description]
		dup_data : [type]    [description]
		dup_data_idx : [type]    [description]
		coord_Map : [type]    [description]
		alg_Map : [type]    [description]
		mol_template : [type]    [description]

		Returns
		-------
		tuple    status,rotmatches,update_to_rdkit
		"""
		Chem.SanitizeMol(mol)

		mol = Chem.AddHs(mol)

		mol.SetProp("_Name", name)

		# detects and applies auto-detection of initial number of conformers
		if args.sample == "auto":
			initial_confs = int(self.auto_sampling(args.auto_sample, mol, args, log))
		else:
			initial_confs = int(args.sample)

		dup_data.at[dup_data_idx, "Molecule"] = name
		update_to_rdkit = False

		rotmatches = getDihedralMatches(mol, args.heavyonly, log)

		if len(rotmatches) > args.max_torsions:
			log.write(
				"x  Too many torsions (%d). Skipping %s"
				% (len(rotmatches), (name + args.output))
			)
			status = -1
		elif args.CSEARCH == "summ" and len(rotmatches) == 0:
			update_to_rdkit = True
			log.write(
				"\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to SUMM SDF"
			)
		elif args.CSEARCH == "fullmonte" and len(rotmatches) == 0:
			update_to_rdkit = True
			log.write(
				"\nx  No rotatable dihedral found. Updating to CSEARCH to RDKit, writing to FULLMONTE SDF"
			)

		# csearch_folder = Path(self.w_dir_initial).joinpath(f"CSEARCH/{args.CSEARCH}")
		# csearch_folder.mkdir(exist_ok=True)
		# csearch_file = csearch_folder.joinpath(name + "_" + args.CSEARCH + args.output)
		# sdwriter = Chem.SDWriter(str(csearch_file))

		# if update_to_rdkit and args.CSEARCH == "summ":
		#     sdwriter = Chem.SDWriter(name + "_" + "summ" + args.output)
		# elif update_to_rdkit and args.CSEARCH == "fullmonte":
		#     sdwriter = Chem.SDWriter(name + "_" + "fullmonte" + args.output)
		# elif args.CSEARCH == "fullmonte":
		#     sdwriter = Chem.SDWriter(name + "_" + "fullmonte" + args.output)
		# elif args.CSEARCH == "crest":
		#     sdwriter = Chem.SDWriter(name + "_" + "crest" + args.output)
		# else:
		#     sdwriter = Chem.SDWriter(name + "_" + "rdkit" + args.output)

		if args.CSEARCH != "crest":
			dup_data.at[dup_data_idx, "RDKit-Initial-samples"] = initial_confs
			if args.CSEARCH == "rdkit":
				rotmatches = []
			cids = self.embed_conf(
				mol, initial_confs, args, log, coord_Map, alg_Map, mol_template
			)

			# energy minimize all to get more realistic results
			# identify the atoms and decide Force Field
			for atom in mol.GetAtoms():
				if (
					atom.GetAtomicNum() > 36
				):  # up to Kr for MMFF, if not the code will use UFF
					log.write(
						"x  "
						+ args.ff
						+ " is not compatible with the molecule, changing to UFF"
					)
					args.ff = "UFF"
			if args.verbose:
				log.write(
					"o  Optimizing "
					+ str(len(cids))
					+ " initial conformers with "
					+ args.ff
				)
				if args.CSEARCH == "summ":
					log.write(
						"o  Found " + str(len(rotmatches)) + " rotatable torsions"
					)
				elif args.CSEARCH == "fullmonte":
					log.write(
						"o  Found " + str(len(rotmatches)) + " rotatable torsions"
					)
				else:
					log.write("o  Systematic torsion rotation is set to OFF")

			status = self.min_after_embed(
				mol,
				cids,
				name,
				initial_confs,
				rotmatches,
				dup_data,
				dup_data_idx,
				sdwriter,
				args,
				log,
				update_to_rdkit,
				coord_Map,
				alg_Map,
				mol_template,
			)
		else:
			args.charge = rules_get_charge(mol, args)
			dup_data.at[dup_data_idx, "Overall charge"] = np.sum(args.charge)
			status = crest_opt(mol, name, dup_data, dup_data_idx, sdwriter, args, log)

		sdwriter.close()

		return status, rotmatches, update_to_rdkit

def smi_to_mol(smi,args,name, constraints_dist, constraints_angle, constraints_dihedral):
	smi = check_for_pieces(smi)
	if len(smi)>1:
		mol, constraints_dist, constraints_angle, constraints_dihedral = nci_ts_mol(smi, args, constraints_dist, constraints_angle, constraints_dihedral, name)
	else:
		mol = Chem.MolFromSmiles(smi[0])
	return mol, constraints_dist, constraints_angle, constraints_dihedral

# MAIN FUNCTION

# main function to generate conformers

def prepare_smiles_files(args, w_dir_initial):
	with open(args.input) as smifile:
		lines = [line for line in smifile if line.strip()]
	job_inputs = []
	for i, line in enumerate(lines):
		smi, name, args, constraints_dist, constraints_angle, constraints_dihedral = prepare_smiles_from_line(line, i, args)
		obj = smi, name, w_dir_initial, args.varfile, args.charge_default, constraints_dist, constraints_angle, constraints_dihedral
		job_inputs.append(obj)
	return job_inputs


def prepare_smiles_from_line(line, i, args):
	toks = line.split()
	# editing part
	smiles = toks[0]
	smi = check_for_pieces(smiles)
	if args.prefix == "None":  # I assume no AttributeError
		name = "".join(toks[1])  # I assume no AttributeError
	else:  # I assume no AttributeError
		name = f"{args.prefix}_{i}_{''.join(toks[1])}"  # I assume no AttributeError
	if len(smi)>1:
		constraints_angle, constraints_dist, constraints_dihedral = None, None, None
		if len(toks) > 2:
			constraints_dist = toks[2]
			constraints_dist = constraints_dist.split('/')
			for i,c in enumerate(constraints_dist):
				constraints_dist[i] = c.split(',')
			if len(toks) > 3:
				constraints_angle = toks[3]
				constraints_angle = constraints_angle.split('/')
				for i,c in enumerate(constraints_angle):
					constraints_angle[i] = c.split(',')
				if len(toks) > 4:
					constraints_dihedral = toks[4]
					constraints_dihedral = constraints_dihedral.split('/')
					for i,c in enumerate(constraints_dihedral):
						constraints_dihedral[i] = c.split(',')
	if args.charge_default == "auto":  # I assume no AttributeError
		if not args.metal_complex:  # I assume no AttributeError
			args.charge_default = check_charge_smi(smi, args.ts_complex)  # I assume no AttributeError
	return smiles, name, args, constraints_dist, constraints_angle, constraints_dihedral

def prepare_csv_files(args, w_dir_initial):
	csv_smiles = pd.read_csv(args.input)
	job_inputs = []
	for i in range(len(csv_smiles)):
		obj = generate_mol_from_csv(args, w_dir_initial, csv_smiles, i)
		job_inputs.append(obj)
	return job_inputs


def generate_mol_from_csv(args, w_dir_initial, csv_smiles, index):
	# assigning names and smi i  each loop
	smi = csv_smiles.loc[index, "SMILES"]
	pruned_smi = check_for_pieces(smi)
	mol = Chem.MolFromSmiles(pruned_smi)
	if args.charge_default == "auto":
		if not args.metal_complex:
			args.charge_default = check_charge_smi(pruned_smi)
	if args.prefix == "None":
		name = csv_smiles.loc[index, "code_name"]
	else:
		name = "comp_" + str(index) + "_" + csv_smiles.loc[index, "code_name"]
	constraints = []
	if 'constraints' in csv_smiles.columns:
		constraints = csv_smiles.loc[index, "constraints"]
	obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
	return obj


def prepare_cdx_files(args, w_dir_initial):
	# converting to smiles from chemdraw
	molecules = generate_mol_from_cdx(args)
	job_inputs = []
	for i, (smi, mol) in enumerate(molecules):
		name = f"{args.input.split('.')[0]}_{str(i)}"
		if args.charge_default == "auto":
			if not args.metal_complex:
				args.charge_default = check_charge_smi(smi)
		constraints = []
		obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
		job_inputs.append(obj)
	return job_inputs


def generate_mol_from_cdx(args):
	cmd_cdx = ["obabel", "-icdx", args.input, "-osmi", "-Ocdx.smi"]
	subprocess.call(cmd_cdx)
	with open("cdx.smi", "r") as smifile:
		smi_lines = [line for line in smifile]
	os.remove("cdx.smi")
	molecules = []
	for smi in smi_lines:
		pruned_smi = check_for_pieces(smi)
		molecule = Chem.MolFromSmiles(pruned_smi)
		molecules.append((pruned_smi, molecule))
	return molecules


def prepare_gaussian_files(args, w_dir_initial):
	job_inputs = []
	charge_com = com_2_xyz_2_sdf(args.input, args.default_charge)
	name = os.path.splitext(args.input)[0]
	sdffile = f"{name}.sdf"
	suppl, _, _ = mol_from_sdf_or_mol_or_mol2(sdffile)

	for i, mol in enumerate(suppl):
		if args.charge_default == "auto":
			args.charge_default = charge_com
		constraints = []
		obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
		job_inputs.append(obj)
	return job_inputs


def prepare_xyz_files(args, w_dir_initial):
	job_inputs = []
	name = os.path.splitext(args.input)[0]
	sdffile = f"{name}.sdf"
	suppl, _, charge_com_list = mol_from_sdf_or_mol_or_mol2(sdffile)
	charge_com = charge_com_list[0]

	for i, mol in enumerate(suppl):
		if args.charge_default == "auto":
			args.charge_default = charge_com
		constraints = []
		obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
		job_inputs.append(obj)
	return job_inputs


def prepare_sdf_file(charge_sdf, w_dir_initial, mol, name, args, i):
	if args.charge_default == "auto":
		args.charge_default = charge_sdf
	constraints = []
	obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
	return obj


def prepare_sdf_files(args, w_dir_initial):
	suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
	job_inputs = []
	for i, (mol, name, charge_sdf) in enumerate(zip(suppl, IDs, charges)):
		obj = prepare_sdf_file(charge_sdf, w_dir_initial, mol, name, args, i)
		job_inputs.append(obj)
	return job_inputs


def prepare_mol_file(suppl, name, charge, args, w_dir_initial):
	if args.charge_default == "auto":
		args.charge_default = charge
	mol = suppl
	constraints= []
	obj = mol, name, w_dir_initial, args.varfile, args.charge_default, constraints
	return obj


def prepare_mol_files(
	args, w_dir_initial
):  # The extra variables are for API Consistency.
	suppl, IDs, charges = mol_from_sdf_or_mol_or_mol2(args.input)
	job_inputs = []
	obj = prepare_mol_file(suppl, IDs[0], charges[0], args, w_dir_initial)
	job_inputs.append(obj)
	return job_inputs
