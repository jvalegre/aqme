#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	       used for genrating details for nmr       #
#####################################################.

from rdkit.Chem import AllChem as Chem
import numpy as np
import os,subprocess

from pyconfort.cheshire_lookup import cheshire

def get_exp_data(args,name,w_dir_initial):
	if args.nmr_exp == 'nmrsdf':
		sdf_file = w_dir_initial+'/'+name+'.sdf'

		sdf_lines = open(sdf_file,'r').readlines()

		found_nmr = 0
		for i,line in enumerate(outlines):
			if line.find('> <NMREDATA_ASSIGNMENT>') > -1:
				start = i+1
				found_nmr = 1
			if len(line.strip()) == 0 and found_nmr == 1 :
				stop = i
				break

		atom_num,atom_sheilding = [],[]
		for i in range(start,stop):
			if len(sdf_lines[idea].split()) == 3:
				atom_num.append(sdf_lines[idea].split()[2])
				atom_sheilding.append(sdf_lines.split()[1])
			if len(sdf_lines[idea].split()) > 3:
				for i in range(2,len(sdf_lines[idea].split()) - 2):
					atom_num.append(sdf_lines[idea].split()[i])
					atom_sheilding.append(sdf_lines.split()[1])

		return atom_num,atom_sheilding

def calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,w_dir_initial,lot_sp,bs_sp,lot,bs):

	#define the lot etc.,
	[opt_method, opt_basis, opt_solv] = [lot,bs,args.solvent_model]
	[nmr_method, nmr_basis, nmr_solv] = [lot_sp,bs_sp,args.solvent_model_sp]

	## Matching against the CHESHIRE database of scaling factors
	if args.nmr_online:
		for nuc in args.nmr_nucleus:
			scale = cheshire(args.online, nuc, opt_method, opt_basis, opt_solv, nmr_method, nmr_basis, nmr_solv, args.nmr_aos)
			try:
				slope.append(float(scale.split()[1]))
				intercept.append(float(scale.split()[3]))
			except AttributeError:
				log.write("   No scaling factors found for this level of theory!"); exit()

			log.write("   The slope for nucleus {0} = {1}".format(nuc, slope))
			log.write("   The intercept for nucleus {0} = {1}".format(nuc, intercept))
		tms_ref = args.nmr_tms_ref

	else:
		slope = args.nmr_slope
		intercept = args.nmr_intercept
		tms_ref = args.nmr_tms_ref

	for file in nmr_log_files:

		# list of H and C shieldings for each individual conformer
		conf_H_shieldings, conf_H_idx, conf_C_shieldings, conf_C_idx, nmr_file_idx = [],[],[],[],[]

		# read the file and detects start and stop point enclosing the NMR section of the LOG file
		start, stop = 0,0

		outfile = open(file,"r")
		outlines = outfile.readlines()

		for i in range(len(outlines)):
			if outlines[i].find('SCF GIAO Magnetic shielding tensor (ppm):') > -1:
				start = i
			if outlines[i].find('Population analysis using the SCF Density') > -1:
				stop = i
				break

		# stores H and C NMR shieldings within the NMR section of the LOG file
		for i in range(start,stop):
			try:
				if '1H' in args.nmr_nucleus:
					if outlines[i].split()[1] == 'H':
						#assigning values from arrays
						index = args.nmr_nucleus.index('1H')
						TMS_H_ref = tms_ref[index]
						slope_H = slope[index]
						intercept_H = intercept[index]

						#conf_H_idx.append('H'+str(outlines[i].split()[0]))
						conf_H_idx.append(str(i-start))
						if TMS_H_ref is not None:
							conf_H_shieldings.append(TMS_H_ref-float(outlines[i].split()[4]))
						else:
							conf_H_shieldings.append((intercept_H-float(outlines[i].split()[4])/(-slope_H)))
				if '13C' in args.nmr_nucleus:
					if outlines[i].split()[1] == 'C':
						#assigning values from arrays
						index = args.nmr_nucleus.index('13C')
						TMS_H_ref = tms_ref[index]
						slope_H = slope[index]
						intercept_H = intercept[index]

						#conf_C_idx.append('C'+str(outlines[i].split()[0]))
						conf_H_idx.append(str(i-start))
						if TMS_C_ref is not None:
							conf_C_shieldings.append(TMS_C_ref-float(outlines[i].split()[4]))
						else:
							conf_C_shieldings.append((intercept_C-float(outlines[i].split()[4])/(-slope_C)))
			except: pass

		outfile.close()

		# get the Boltzmann probabilities of each conformer from a GoodVibes file
		boltz_file = w_dir_initial+'/QPRED/nmr/boltz/'+str(lot)+'-'+str(bs) + '/Goodvibes_'+name+'.dat'

		boltz_outfile = open(boltz_file,"r")
		boltz_outlines = boltz_outfile.readlines()

		for i in range(len(boltz_outlines)):
			# I remove the NMR from the file names using [0:-4]
			if boltz_outlines[i].find(file.split('._NMR')[0]) > -1:
				boltz_factor = float(boltz_outlines[i].split()[-1])
				break

		# Multiply the shieldings by the corresponding Boltzmann factor
		if '1H' in args.nmr_nucleus:
			conf_H_shieldings = np.asarray(conf_H_shieldings, dtype=np.float64)*boltz_factor
			final_H_shieldings.append(conf_H_shieldings)
			# convert list of lists into numpy arrays
			final_H_shieldings = np.array(final_H_shieldings)
			# Sum all the individual H shieldings * Boltzmann prob
			final_H_shieldings = np.sum(final_H_shieldings, axis=0)
			# convert all the idx lists into numpy arrays to concatenate below
			conf_H_idx = np.array(conf_H_idx)

		if '13C' in args.nmr_nucleus:
			conf_C_shieldings = np.asarray(conf_C_shieldings, dtype=np.float64)*boltz_factor
			final_C_shieldings.append(conf_C_shieldings)
			final_C_shieldings = np.array(final_C_shieldings)
			final_C_shieldings = np.sum(final_C_shieldings,axis=0)
			conf_C_idx = np.array(conf_C_idx)

	# concatenate H and C data one after the other (for printing in the CSV)
	if final_C_shieldings is not None:
		all_shieldings = np.concatenate((final_H_shieldings, final_C_shieldings))
		all_idx = np.concatenate((conf_H_idx, conf_C_idx))


	#getting experimental shiftd for get_atom
	atom_num_exp,atom_sheilding_exp = get_exp_data(args,name,w_dir_initial)

	#make final array with atom_num, exp, dft
	df = pd.DataFrame(list(zip(atom_num_exp,atom_sheilding_exp)),
			   columns =['Atom', 'Shielding-Exp'])

	for i, atom_num in enumerate(atom_num_exp):
		for j,atom_cal in enumerate(all_idx):
			if atom_num == atom_cal:
				df.at[i,'Shielding-Cal'] = all_shieldings[j]

	#calculation of MAE and SD


	# create a CSV file containing all the data for the dr
	df = pd.DataFrame(list(zip(all_idx,all_shieldings)),
			   columns =['Atom', 'Shielding'])

	w_dir = os.getcwd()
	folder = w_dir_initial+'/QPRED/nmr/average-nmr/'+str(lot_sp)+'-'+str(bs_sp)
	log.write("\no  Preparing final NMR results in {}".format(folder))
	try:
		os.makedirs(folder)
		os.chdir(folder)
	except OSError:
		if os.path.isdir(folder):
			os.chdir(folder)
		else:
			pass
	export_param_excel = df.to_csv(name+'_'+lot_sp+'-'+bs_sp+'_predicted_shifts.csv', index = None, header=False)

def calculate_boltz_and_nmr(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	destination = w_dir_initial+'/QPRED/nmr/boltz/'+str(lot)+'-'+str(bs)
	moving_files('Goodvibes_'+name+'.dat', destination)


	for lot_sp in args.level_of_theory_sp:
		for bs_sp in args.basis_set_sp:
			for bs_gcp_sp in args.basis_set_genecp_atoms_sp:
				#performing the nmr calculations
				dir_sp_nmr =  w_dir_fin+'/output-files/G16-SP_input_files/'+str(lot_sp)+'-'+str(bs_sp)
				os.chdir(dir_sp_nmr)
				#grabbing the respective NMR files for a given molecules
				nmr_log_files = glob.glob(name+'*_NMR.log') + glob.glob(name+'*_NMR_*.log')
				calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,w_dir_initial,lot_sp,bs_sp,lot,bs)
