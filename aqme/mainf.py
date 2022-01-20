#####################################################.
#       This file stores all the main functions     #
#####################################################.

import glob
import sys
import os
import time
import concurrent.futures as futures
import multiprocessing as mp
from pathlib import Path

import pandas as pd
from progress.bar import IncrementalBar
from rdkit.Chem import AllChem as Chem

# from aqme.csearch import (check_charge_smi, clean_args,
#                                compute_confs,
#                                mol_from_sdf_or_mol_or_mol2, creation_of_dup_csv)
from aqme.csearch import *
from aqme.filter import geom_rules_output

from aqme.qprep import qprep, get_molecule_list, load_charge_data

# from aqme.grapher import graph
# from aqme.descp import calculate_parameters
# from aqme.nmr import calculate_boltz_and_nmr
# from aqme.energy import calculate_boltz_and_energy,calculate_avg_and_energy
# from aqme.dbstep_conf import calculate_db_parameters,calculate_boltz_and_dbstep
# from aqme.nics_conf import calculate_nics_parameters,calculate_boltz_for_nics,calculate_avg_nics
# from aqme.cclib_conf import calculate_cclib,get_avg_cclib_param,calculate_boltz_for_cclib
from aqme.cmin import mult_min
from aqme.utils import (
	Logger,
	move_file,
	get_filenames,
	creation_of_dup_csv_csearch,
	creation_of_dup_csv_cmin,
	read_energies,
	get_name_and_charge,
)

# need to and in energy


def csearch_main(w_dir_initial, args, log_overall):

	file_format = os.path.splitext(args.input)[1]

	# Checks
	if file_format not in SUPPORTED_INPUTS:
		log_overall.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!")
		sys.exit()
	if not os.path.exists(args.input):
		log_overall.write("\nx  INPUT FILE NOT FOUND!")
		sys.exit()

	# if large system increase stack size
	# if args.STACKSIZE != '1G':
	#     os.environ['OMP_STACKSIZE'] = args.STACKSIZE
	smi_derivatives = [".smi", ".txt", ".yaml", ".yml", ".rtf"]
	Extension2inputgen = dict()
	for key in smi_derivatives:
		Extension2inputgen[key] = prepare_smiles_files
	Extension2inputgen[".csv"] = prepare_csv_files
	Extension2inputgen[".cdx"] = prepare_cdx_files
	Extension2inputgen[".gjf"] = prepare_gaussian_files
	Extension2inputgen[".com"] = prepare_gaussian_files
	Extension2inputgen[".xyz"] = prepare_gaussian_files
	Extension2inputgen[".sdf"] = prepare_sdf_files
	Extension2inputgen[".mol"] = prepare_mol_files
	Extension2inputgen[".mol2"] = prepare_mol_files

	with futures.ProcessPoolExecutor(
		max_workers=args.cpus, mp_context=mp.get_context("fork")
	) as executor:
		# Submit a set of asynchronous jobs
		jobs = []
		count_mol = 0

		# Prepare the Jobs
		prepare_function = Extension2inputgen[file_format]
		job_inputs = prepare_function(args, w_dir_initial)

		# Submit the Jobs
		for job_input in job_inputs:
			smi_, name_, dir_, varfile_, charge_default_, constraints_dist_, constraints_angle_, constraints_dihedral_ = job_input
			job = executor.submit(
				process_csearch, smi_, name_, dir_, varfile_, charge_default_, constraints_dist_, constraints_angle_, constraints_dihedral_
			)
			jobs.append(job)
			count_mol += 1

		final_dup_data = creation_of_dup_csv_csearch(args.CSEARCH)
		bar = IncrementalBar("o  Number of finished jobs from CSEARCH", max=count_mol)
		# Process the job results (in submission order) and save the conformers.
		for i, job in enumerate(jobs):
			total_data = job.result()
			frames = [final_dup_data, total_data]
			final_dup_data = pd.concat(frames, ignore_index=True, sort=True)
			bar.next()
		bar.finish()

		# removing temporary files
		temp_files = [
			"gfn2.out",
			"xTB_opt.traj",
			"ANI1_opt.traj",
			"wbo",
			"xtbrestart",
			"ase.opt",
			"xtb.opt",
			"gfnff_topo",
		]
		for file in temp_files:
			if os.path.exists(file):
				os.remove(file)

	return final_dup_data


def cmin_main(w_dir_initial, args, log_overall, dup_data):
	bar = IncrementalBar("o  Number of finished jobs from CMIN", max=len(dup_data))
	final_dup_data = creation_of_dup_csv_cmin(args.CMIN)
	for dup_data_idx in range(len(dup_data)):
		# update_to_rdkit = dup_data.at[dup_data_idx,'update_to_rdkit']
		name = dup_data.at[dup_data_idx, "Molecule"]
		charge = dup_data.at[dup_data_idx, "Overall charge"]
		if dup_data.at[dup_data_idx, "status"] != -1:
			if args.CMIN == "ani":
				min_suffix = "ani"
			elif args.CMIN == "xtb":
				min_suffix = "xtb"
			if args.CSEARCH in ["rdkit", "summ", "fullmonte"]:

				csearch_folder = Path(w_dir_initial).joinpath(f"CSEARCH/{args.CSEARCH}")
				fullname = str(csearch_folder.joinpath(name + "_" + args.CSEARCH))

				# fullname = f'{name}_{args.CSEARCH}'
				# try:
				total_data = mult_min(
					fullname, args, min_suffix, charge, log_overall, w_dir_initial
				)
				# except:
				#     pass
				frames = [final_dup_data, total_data]
				final_dup_data = pd.concat(frames, ignore_index=True, sort=True)
		bar.next()
	bar.finish()

	# removing temporary files
	temp_files = [
		"gfn2.out",
		"xTB_opt.traj",
		"ANI1_opt.traj",
		"wbo",
		"xtbrestart",
		"ase.opt",
		"xtb.opt",
		"gfnff_topo",
	]
	for file in temp_files:
		if Path(file).exists():
			os.remove(file)

	return final_dup_data


# MAIN QPREP FUNCTION
def qprep_main(w_dir_initial, args, log):

	if len(args.geom_rules) >= 1:
		conf_files = glob.glob("*_rules.sdf")
	# define the SDF files to convert to COM Gaussian files
	elif args.CMIN == "xtb":
		conf_files = glob.glob(w_dir_initial + "/CMIN/" + args.CMIN + "/*_xtb.sdf")
	elif args.CMIN == "ani":
		conf_files = glob.glob(w_dir_initial + "/CMIN/" + args.CMIN + "/*_ani.sdf")
	elif args.CSEARCH == "rdkit":
		conf_files = glob.glob(
			w_dir_initial + "/CSEARCH/" + args.CSEARCH + "/*_rdkit.sdf"
		)
	elif args.CSEARCH == "summ":
		conf_files = glob.glob(
			w_dir_initial + "/CSEARCH/" + args.CSEARCH + "/*_summ.sdf"
		)
	elif args.CSEARCH == "fullmonte":
		conf_files = glob.glob(
			w_dir_initial + "/CSEARCH/" + args.CSEARCH + "/*_fullmonte.sdf"
		)
	elif args.CSEARCH == "crest":
		conf_files = glob.glob(
			w_dir_initial + "/CSEARCH/" + args.CSEARCH + "/*_crest.sdf"
		)
	else:
		conf_files = glob.glob("*.sdf")

	# # NEED TO UPDATE THIS PART TO START FROM JSON!
	# if args.com_from_xyz:
	#     xyz_files = glob.glob("*.xyz")
	#     for file in xyz_files:
	#         mol = next(pybel.readfile("xyz", file))
	#         stem = Path(file).stem
	#         mol.write("sdf", f"{stem}.sdf")
	#     conf_files = glob.glob("*.sdf")

	if not conf_files:
		log.write("\nx  No SDF files detected to convert to gaussian COM files")
		return

	csv_name = args.input.split(".")[0]
	csv_file = f"{w_dir_initial}/CSEARCH/csv_files/{csv_name}-CSEARCH-Data.csv"
	charge_data, invalid_files = load_charge_data(csv_file, conf_files)

	# remove the invalid files and non-existing files
	accept_file = lambda x: x not in invalid_files and Path(x).exists()
	conf_files = [file for file in conf_files if accept_file(file)]

	# Prepare the list of molecules that are to be written
	mols = []

	for file in conf_files:
		filepath = f"{file}"
		new_mols = get_molecule_list(
			filepath,
			lowest_only=args.lowest_only,
			lowest_n=args.lowest_n,
			energy_threshold=args.energy_threshold_for_gaussian,
		)
		mols.extend(new_mols)

		name = os.path.basename(filepath).split(".")[0].split('_')[0]
		charge = charge_data[charge_data.Molecule == name]['Overall charge'][0]

		# writing the com files
		for i, mol in enumerate(mols):
			qprep(mol=mol, molecule=name+'_conf_'+str(i+1), charge=charge, atom_types = [], yaml_file=args.varfile)

      
# moving files after compute and/or write_gauss
def move_sdf_main(args):
	src = Path(os.getcwd())
	if len(args.geom_rules) >= 1:
		geom_rules_files = glob.glob("*_filter_geom_rules.sdf")
	if args.CMIN == "xtb":
		all_xtb_conf_files = glob.glob("*_xtb.sdf")
		destination_xtb = src.joinpath("CMIN/xtb/")
		for file in all_xtb_conf_files:
			move_file(destination_xtb, src, file)
		all_xtb_conf_files_all = glob.glob("*_xtb_all_confs.sdf")
		destination_xtb_all = src.joinpath("CMIN/xtb_all_confs/")
		for file in all_xtb_conf_files_all:
			move_file(destination_xtb_all, src, file)
		if len(args.geom_rules) >= 1:
			destination_geom_rules = src.joinpath("CMIN/xtb/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)
	if args.CMIN == "ani":
		all_ani_conf_files = glob.glob("*_ani.sdf")
		destination_ani = src.joinpath("CMIN/ani/")
		for file in all_ani_conf_files:
			move_file(destination_ani, src, file)
		all_ani_conf_files_all = glob.glob("*_ani_all_confs.sdf")
		destination_ani_all = src.joinpath("CMIN/ani_all_confs/")
		for file in all_ani_conf_files_all:
			move_file(destination_ani_all, src, file)
		if len(args.geom_rules) >= 1:
			destination_geom_rules = src.joinpath("CMIN/ani/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)
	if args.CSEARCH == "rdkit":
		all_name_conf_files = glob.glob("*_rdkit.sdf")
		destination_rdkit = src.joinpath("CSEARCH/rdkit/")
		for file in all_name_conf_files:
			move_file(destination_rdkit, src, file)
		if len(args.geom_rules) >= 1:
			destination_geom_rules = src.joinpath("CSEARCH/rdkit/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)

	if args.CSEARCH == "summ":
		all_name_conf_files = glob.glob("*_summ.sdf")
		destination_rdkit = src.joinpath("CSEARCH/summ/")
		for file in all_name_conf_files:
			move_file(destination_rdkit, src, file)
		if len(args.geom_rules) >= 1 and args.CMIN is None:
			destination_geom_rules = src.joinpath("CSEARCH/summ/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)

	if args.CSEARCH == "fullmonte":
		all_name_conf_files = glob.glob("*_fullmonte.sdf")
		destination_rdkit = src.joinpath("CSEARCH/fullmonte/")
		for file in all_name_conf_files:
			move_file(destination_rdkit, src, file)
		if len(args.geom_rules) >= 1 and args.CMIN is None:
			destination_geom_rules = src.joinpath("CSEARCH/fullmonte/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)

	if args.CSEARCH == "crest":
		all_name_conf_files = glob.glob("*_crest.sdf")
		destination_rdkit = src.joinpath("CSEARCH/crest/")
		for file in all_name_conf_files:
			move_file(destination_rdkit, src, file)
		if len(args.geom_rules) >= 1 and args.CMIN is None:
			destination_geom_rules = src.joinpath("CSEARCH/crest/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)

	if args.CSEARCH is None:
		if len(args.geom_rules) >= 1:
			destination_geom_rules = src.joinpath("QCALC/SDF_input/filter_geom_rules/")
			for file in geom_rules_files:
				move_file(destination_geom_rules, src, file)

		all_conf_files = glob.glob("*.sdf")
		destination = src.joinpath("QCALC/SDF_input/")
		for file in all_conf_files:
			move_file(destination, src, file)

	if args.com_from_xyz:
		all_xyz_conf_files = glob.glob("*.xyz") + glob.glob("*.sdf")
		destination_xyz = src.joinpath("QCALC/xyz_and_sdf/")
		for file in all_xyz_conf_files:
			move_file(destination_xyz, src, file)


# getting descriptors
def geom_par_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no  Calculating paramters for molecule : {0} ".format(name))

		sdf_ani, sdf_xtb = None, None
		if os.path.exists(w_dir_initial + "/CSEARCH/rdkit/" + name + "_rdkit.sdf"):
			sdf_rdkit = w_dir_initial + "/CSEARCH/rdkit/" + name + "_rdkit.sdf"
		elif os.path.exists(w_dir_initial + "/CSEARCH/summ/" + name + "_summ.sdf"):
			sdf_rdkit = w_dir_initial + "/CSEARCH/summ/" + name + "_summ.sdf"
		elif os.path.exists(
			w_dir_initial + "/CSEARCH/fullmonte/" + name + "_fullmonte.sdf"
		):
			sdf_rdkit = w_dir_initial + "/CSEARCH/fullmonte/" + name + "_fullmonte.sdf"
		if os.path.exists(w_dir_initial + "/CMIN/xtb/" + name + "_xtb.sdf"):
			sdf_xtb = w_dir_initial + "/CMIN/xtb/" + name + "_xtb.sdf"
		if os.path.exists(w_dir_initial + "/CMIN/ani/" + name + "_ani.sdf"):
			sdf_ani = w_dir_initial + "/CMIN/ani/" + name + "_ani.sdf"
		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
			# Sets the folder and find the log files to analyze
			for lot, bs, bs_gcp in zip(
				args.level_of_theory, args.basis_set, args.genecp_bs
			):
				# assign the path to the finished directory.
				if str(bs).find("/") > -1:
					w_dir = (
						args.path
						+ str(lot)
						+ "-"
						+ str(bs).split("/")[0]
						+ "/success/output_files"
					)
				else:
					w_dir = (
						args.path + str(lot) + "-" + str(bs) + "/success/output_files"
					)
				os.chdir(w_dir)
				qm_files = get_filenames("output", name)
				calculate_parameters(
					sdf_rdkit,
					sdf_ani,
					sdf_xtb,
					qm_files,
					args,
					log,
					w_dir_initial,
					name,
					lot,
					bs,
				)
		else:
			calculate_parameters(
				sdf_rdkit,
				sdf_ani,
				sdf_xtb,
				None,
				args,
				log,
				w_dir_initial,
				name,
				None,
				None,
			)

		os.chdir(w_dir_initial)


# function to plot graphs
def graph_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no  Plotting graphs for molecule : {0} ".format(name))

		sdf_ani, sdf_xtb = None, None
		if os.path.exists(w_dir_initial + "/CSEARCH/rdkit/" + name + "_rdkit.sdf"):
			sdf_rdkit = w_dir_initial + "/CSEARCH/rdkit/" + name + "_rdkit.sdf"
		elif os.path.exists(w_dir_initial + "/CSEARCH/summ/" + name + "_summ.sdf"):
			sdf_rdkit = w_dir_initial + "/CSEARCH/summ/" + name + "_summ.sdf"
		elif os.path.exists(
			w_dir_initial + "/CSEARCH/fullmonte/" + name + "_fullmonte.sdf"
		):
			sdf_rdkit = w_dir_initial + "/CSEARCH/fullmonte/" + name + "_fullmonte.sdf"
		if os.path.exists(
			w_dir_initial + "/CMIN/xtb_all_confs/" + name + "_xtb_all_confs.sdf"
		):
			sdf_xtb = (
				w_dir_initial + "/CMIN/xtb_all_confs/" + name + "_xtb_all_confs.sdf"
			)
		if os.path.exists(
			w_dir_initial + "/CMIN/ani_all_confs/" + name + "_ani_all_confs.sdf"
		):
			sdf_ani = (
				w_dir_initial + "/CMIN/ani_all_confs/" + name + "_ani_all_confs.sdf"
			)
		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
			# Sets the folder and find the log files to analyze
			for lot, bs, bs_gcp in zip(
				args.level_of_theory, args.basis_set, args.genecp_bs
			):
				# assign the path to the finished directory.
				w_dir = args.path + str(lot) + "-" + str(bs) + "/success/output_files"
				os.chdir(w_dir)
				qm_files = get_filenames("output", name)
				if os.path.exists(
					args.path + str(lot) + "-" + str(bs) + "/success/G16-SP_input_files"
				):
					for lot_sp, bs_sp, bs_gcp_sp in zip(
						args.level_of_theory_sp, args.basis_set_sp, args.gen_bs_sp
					):
						w_dir_sp = (
							args.path
							+ str(lot)
							+ "-"
							+ str(bs)
							+ "/success/G16-SP_input_files"
							+ "/"
							+ str(lot_sp)
							+ "-"
							+ str(bs_sp)
						)
						sp_files = get_filenames("output", name)
						graph(
							sdf_rdkit,
							sdf_xtb,
							sdf_ani,
							qm_files,
							sp_files,
							args,
							log,
							lot,
							bs,
							lot_sp,
							bs_sp,
							name,
							w_dir_initial,
							w_dir_sp,
							w_dir,
							"g16",
						)
				if os.path.exists(
					f"{args.path}{lot}-{bs}/success/TURBOMOLE-SP_input_files"
				):
					for lot_sp, bs_sp, bs_gcp_sp in zip(
						args.level_of_theory_sp, args.basis_set_sp, args.gen_bs_sp
					):
						w_dir_sp = f"{path_turbomole}/{lot_sp}-{bs_sp.split('/')[0]}"
						os.chdir(w_dir_sp)
						sp_files = []
						for path in Path(w_dir_sp).iterdir():
							if path.is_dir() and "_SP" in path.stem:
								sp_files.append(path)
						os.chdir(w_dir)
						graph(
							sdf_rdkit,
							sdf_xtb,
							sdf_ani,
							qm_files,
							sp_files,
							args,
							log,
							lot,
							bs,
							lot_sp,
							bs_sp.split("/")[0],
							name,
							w_dir_initial,
							w_dir_sp,
							w_dir,
							"turbomole",
						)

				if os.path.exists(
					args.path
					+ str(lot)
					+ "-"
					+ str(bs)
					+ "/success/ORCA-SP_input_files"
				):
					for lot_sp, bs_sp, bs_gcp_sp in zip(
						args.level_of_theory_sp, args.basis_set_sp, args.gen_bs_sp
					):
						w_dir_sp = (
							args.path
							+ str(lot)
							+ "-"
							+ str(bs)
							+ "/success/ORCA-SP_input_files"
							+ "/"
							+ str(lot_sp)
							+ "-"
							+ str(bs_sp.split("/")[0])
						)
						os.chdir(w_dir_sp)
						sp_files = get_filenames("output", name)
						os.chdir(w_dir)
						graph(
							sdf_rdkit,
							sdf_xtb,
							sdf_ani,
							qm_files,
							sp_files,
							args,
							log,
							lot,
							bs,
							lot_sp,
							bs_sp.split("/")[0],
							name,
							w_dir_initial,
							w_dir_sp,
							w_dir,
							"orca",
						)
				else:
					graph(
						sdf_rdkit,
						sdf_xtb,
						sdf_ani,
						qm_files,
						None,
						args,
						log,
						lot,
						bs,
						None,
						None,
						name,
						w_dir_initial,
						None,
						w_dir,
						None,
					)

		else:
			graph(
				sdf_rdkit,
				sdf_xtb,
				sdf_ani,
				None,
				None,
				args,
				log,
				None,
				None,
				None,
				None,
				name,
				w_dir_initial,
				None,
				None,
				None,
			)

	os.chdir(w_dir_initial)


# function for comparison of nmr
def nmr_main(args, log, w_dir_initial):

	if os.path.exists(w_dir_initial + "/QCALC/G16"):
		args.path = w_dir_initial + "/QCALC/G16/"
	else:
		log.write(
			"\nx  The path for NMR analysis was not set up properly! (check the tutorials for more information)"
		)

	# get sdf FILES from csv
	try:
		pd_name = pd.read_csv(
			w_dir_initial
			+ "/CSEARCH/csv_files/"
			+ args.input.split(".")[0]
			+ "-CSEARCH-Data.csv"
		)

	except FileNotFoundError:
		# detects all the unique molecules from the success folder
		if str(args.basis_set[0]).find("/") > -1:
			w_dir_fin = (
				args.path
				+ str(args.level_of_theory[0])
				+ "-"
				+ str(args.basis_set[0]).split("/")[0]
				+ "/success/output_files"
			)
		else:
			w_dir_fin = (
				args.path
				+ str(args.level_of_theory[0])
				+ "-"
				+ str(args.basis_set[0])
				+ "/success/output_files"
			)
		os.chdir(w_dir_fin)

		nmr_list = []
		standard_suffixes = ["xtb", "ani"]
		for name_nmr in glob.glob("*.*"):
			# discard_charact keeps track of the extra characters after the name of the molecule
			discard_charact = len(name_nmr.split(".")[1]) + 1
			potential_unique = name_nmr.split(".")[0].split("_")
			for i in reversed(range(len(potential_unique))):
				try:
					if potential_unique[i] not in standard_suffixes:
						int(potential_unique[i])
					discard_charact += len(potential_unique[i])
					discard_charact += 1
				except ValueError:
					original_name = name_nmr[:-discard_charact]
					if original_name not in nmr_list:
						nmr_list.append(original_name)
		pd_name = pd.DataFrame(data=nmr_list, columns=["Molecule"])

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no NMR analysis for molecule : {0} ".format(name))

		# Sets the folder and find the log files to analyze
		for lot, bs, bs_gcp in zip(
			args.level_of_theory, args.basis_set, args.genecp_bs
		):
			# assign the path to the finished directory.
			if str(bs).find("/") > -1:
				w_dir_fin = (
					args.path
					+ str(lot)
					+ "-"
					+ str(bs).split("/")[0]
					+ "/success/output_files"
				)
			else:
				w_dir_fin = (
					args.path + str(lot) + "-" + str(bs) + "/success/output_files"
				)
			os.chdir(w_dir_fin)

			qm_files = get_filenames("output", name)
			if len(qm_files) != 0:
				calculate_boltz_and_nmr(
					qm_files, args, log, name, w_dir_fin, w_dir_initial, lot, bs
				)
	os.chdir(w_dir_initial)


def energy_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write(
			"\no Boltzmann average energy analysis for molecule : {0} ".format(name)
		)
		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
		# Sets the folder and find the log files to analyze
		for lot, bs, bs_gcp in zip(
			args.level_of_theory, args.basis_set, args.genecp_bs
		):
			# assign the path to the finished directory.
			if str(bs).find("/") > -1:
				w_dir_fin = (
					args.path
					+ str(lot)
					+ "-"
					+ str(bs).split("/")[0]
					+ "/success/output_files/"
				)
			else:
				w_dir_fin = (
					args.path + str(lot) + "-" + str(bs) + "/success/output_files/"
				)
			os.chdir(w_dir_fin)
			qm_files = get_filenames("output", name)
			if len(qm_files) != 0:
				calculate_boltz_and_energy(
					qm_files, args, log, name, w_dir_fin, w_dir_initial, lot, bs
				)

	# combining the combining all files in different folders
	w_dir_boltz = w_dir_initial + "/QPRED/energy/boltz/"

	for lot, bs, bs_gcp in zip(args.level_of_theory, args.basis_set, args.genecp_bs):
		# assign the path to the finished directory.
		if str(bs).find("/") > -1:
			w_dir_fin = w_dir_boltz + str(lot) + "-" + str(bs).split("/")[0]
		else:
			w_dir_fin = w_dir_boltz + str(lot) + "-" + str(bs)
		os.chdir(w_dir_fin)
		dat_files = glob.glob("*.dat")
		if len(dat_files) != 0:
			calculate_avg_and_energy(
				dat_files,
				args,
				log,
				name,
				w_dir_fin,
				w_dir_initial,
				w_dir_boltz,
				lot,
				bs,
			)

	os.chdir(w_dir_initial)


def dbstep_par_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no  Calculating paramters for molecule : {0} ".format(name))

		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
			# Sets the folder and find the log files to analyze
			for lot, bs, bs_gcp in zip(
				args.level_of_theory, args.basis_set, args.genecp_bs
			):
				# assign the path to the finished directory.
				if str(bs).find("/") > -1:
					w_dir = (
						args.path
						+ str(lot)
						+ "-"
						+ str(bs).split("/")[0]
						+ "/success/output_files"
					)
				else:
					w_dir = (
						args.path + str(lot) + "-" + str(bs) + "/success/output_files"
					)
				os.chdir(w_dir)
				qm_files = get_filenames("output", name)
				calculate_db_parameters(
					qm_files, args, log, w_dir_initial, name, lot, bs
				)
				calculate_boltz_and_dbstep(
					qm_files, args, log, name, w_dir, w_dir_initial, lot, bs
				)
		os.chdir(w_dir_initial)


def nics_par_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no  Calculating nics for molecule : {0} ".format(name))

		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
			# Sets the folder and find the log files to analyze
			for lot, bs, bs_gcp in zip(
				args.level_of_theory, args.basis_set, args.genecp_bs
			):
				# assign the path to the finished directory.
				if str(bs).find("/") > -1:
					w_dir = (
						args.path
						+ str(lot)
						+ "-"
						+ str(bs).split("/")[0]
						+ "/success/output_files"
					)
				else:
					w_dir = (
						args.path + str(lot) + "-" + str(bs) + "/success/output_files"
					)
				os.chdir(w_dir)
				qm_files = get_filenames("output", name)
				# do boltz firsst
				calculate_boltz_for_nics(
					qm_files, args, log, name, w_dir, w_dir_initial, lot, bs
				)
				for lot_sp, bs_sp, bs_gcp_sp in zip(
					args.level_of_theory_sp, args.basis_set_sp, args.gen_bs_sp
				):
					if str(bs).find("/") > -1:
						w_dir_sp = (
							args.path
							+ str(lot)
							+ "-"
							+ str(bs).split("/")[0]
							+ "/success/G16-NICS_input_files/"
							+ str(lot_sp)
							+ "-"
							+ str(bs_sp)
						)
					else:
						w_dir_sp = (
							args.path
							+ str(lot)
							+ "-"
							+ str(bs)
							+ "/success/G16-NICS_input_files/"
							+ str(lot_sp)
							+ "-"
							+ str(bs_sp)
						)
					os.chdir(w_dir_sp)
					qm_files_sp = get_filenames("output", name)
					calculate_nics_parameters(
						qm_files_sp, args, log, w_dir_initial, name, lot_sp, bs_sp
					)
					calculate_avg_nics(
						qm_files_sp,
						args,
						log,
						name,
						w_dir_sp,
						w_dir_initial,
						lot_sp,
						bs_sp,
					)
		os.chdir(w_dir_initial)


def cclib_main(args, log, w_dir_initial):
	# get sdf FILES from csv
	pd_name = pd.read_csv(
		w_dir_initial
		+ "/CSEARCH/csv_files/"
		+ args.input.split(".")[0]
		+ "-CSEARCH-Data.csv"
	)

	for i in range(len(pd_name)):
		name = pd_name.loc[i, "Molecule"]

		log.write("\no  Calculating cclib paramters for molecule : {0} ".format(name))
		if os.path.exists(w_dir_initial + "/QCALC/G16"):
			args.path = w_dir_initial + "/QCALC/G16/"
			# Sets the folder and find the log files to analyze
			for lot, bs, bs_gcp in zip(
				args.level_of_theory, args.basis_set, args.genecp_bs
			):
				# assign the path to the finished directory.
				if str(bs).find("/") > -1:
					w_dir = (
						args.path
						+ str(lot)
						+ "-"
						+ str(bs).split("/")[0]
						+ "/success/output_files"
					)
				else:
					w_dir = (
						args.path + str(lot) + "-" + str(bs) + "/success/output_files"
					)
				os.chdir(w_dir)
				qm_files = get_filenames("output", name)
				# do boltz firsst
				calculate_cclib(qm_files, w_dir_initial, lot, bs)
				calculate_boltz_for_cclib(qm_files, name, w_dir_initial, lot, bs)
				if str(bs).find("/") > -1:
					os.chdir(
						w_dir_initial
						+ "/QPRED/cclib-json/all_confs_cclib/"
						+ str(lot)
						+ "-"
						+ str(bs).split("/")[0]
					)
				else:
					os.chdir(
						w_dir_initial
						+ "/QPRED/cclib-json/all_confs_cclib/"
						+ str(lot)
						+ "-"
						+ str(bs)
					)
				json_files = get_filenames("output", name)
				get_avg_cclib_param(json_files, name, w_dir_initial, lot, bs)


# MAIN OPTION FOR DISCARDING MOLECULES BASED ON USER INPUT DATA (REFERRED AS EXPERIMENTAL RULES)
def geom_rules_main(args, log, geom_rules_active):
	if geom_rules_active:
		if args.verbose:
			log.write(
				"\n   ----- Applying experimental rules to write the new confs file -----"
			)
		# do 2 cases, for RDKit only and RDKIt+xTB
		if args.CMIN == "xtb":
			conf_files = glob.glob("*_xtb.sdf")
		if args.CMIN == "ani":
			conf_files = glob.glob("*_ani.sdf")
		if args.CMIN is None:
			if args.CSEARCH == "rdkit":
				conf_files = glob.glob("*_rdkit.sdf")
			elif args.CSEARCH == "summ":
				conf_files = glob.glob("*_summ.sdf")
			elif args.CSEARCH == "fullmonte":
				conf_files = glob.glob("*_fullmonte.sdf")
			else:
				conf_files = glob.glob("*.sdf")

		for file in conf_files:
			try:
				allmols = Chem.SDMolSupplier(file, removeHs=False)
			except OSError:
				pass
			if allmols:
				sdwriter = Chem.SDWriter(file.split(".")[0] + "_filter_geom_rules.sdf")
				print_error_geom_rules = False
				for mol in allmols:
					check_mol = True
					ob_compat = True
					rdkit_compat = True
					check_mol = geom_rules_output(
						mol, args, log, file, print_error_geom_rules
					)
					print_error_geom_rules += 1
					if check_mol:
						sdwriter.write(mol)

				sdwriter.close()
