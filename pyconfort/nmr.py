#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	       used for genrating details for nmr       #
#####################################################.

from rdkit.Chem import AllChem as Chem
import numpy as np
import os,subprocess

def calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,lot_sp,bs_sp):

	#some constants which can be changed
	TMS_H_ref = 0
	TMS_C_ref = 0

	slope_H = -1.0565
	intercept_H = 31.9340
	slope_C = -1.0427
	intercept_C = 181.7173

	for file in nmr_log_files:

		# list of H and C shieldings for each individual conformer
		conf_H_shieldings, conf_H_idx, conf_C_shieldings, conf_C_idx = [],[],[],[]

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
				if outlines[i].split()[1] == 'H':
					conf_H_idx.append('H'+str(outlines[i].split()[0]))
					if TMS_H_ref is not None:
						conf_H_shieldings.append(TMS_H_ref-float(outlines[i].split()[4]))
					else:
						conf_H_shieldings.append((intercept_H-float(outlines[i].split()[4])/(-slope_H)))

				if outlines[i].split()[1] == 'C':
					conf_C_idx.append('C'+str(outlines[i].split()[0]))
					if TMS_C_ref is not None:
						conf_C_shieldings.append(TMS_C_ref-float(outlines[i].split()[4]))
					else:
						conf_C_shieldings.append((intercept_C-float(outlines[i].split()[4])/(-slope_C)))

			except: pass

		outfile.close()

		# get the Boltzmann probabilities of each conformer from a GoodVibes file
		boltz_file = w_dir_fin + '/Goodvibes_'+name+'.dat'

		boltz_outfile = open(boltz_file,"r")
		boltz_outlines = boltz_outfile.readlines()

		for i in range(len(boltz_outlines)):
			# I remove the NMR from the file names using [0:-4]
			if boltz_outlines[i].find(file.split('._NMR')[0]) > -1:
				boltz_factor = float(boltz_outlines[i].split()[-1])
				break

		# Multiply the shieldings by the corresponding Boltzmann factor
		conf_H_shieldings = np.asarray(conf_H_shieldings, dtype=np.float64)*boltz_factor
		conf_C_shieldings = np.asarray(conf_C_shieldings, dtype=np.float64)*boltz_factor

		final_H_shieldings.append(conf_H_shieldings)
		final_C_shieldings.append(conf_C_shieldings)

	# convert list of lists into numpy arrays
	final_H_shieldings = np.array(final_H_shieldings)
	final_C_shieldings = np.array(final_C_shieldings)

	# Sum all the individual H shieldings * Boltzmann prob
	final_H_shieldings = np.sum(final_H_shieldings, axis=0)
	final_C_shieldings = np.sum(final_C_shieldings,axis=0)

	# convert all the idx lists into numpy arrays to concatenate below
	conf_H_idx = np.array(conf_H_idx)
	conf_C_idx = np.array(conf_C_idx)

#     # this part is just in case you want all the info together
#     H_shieldings.append(final_H_shieldings)
#     C_shieldings.append(final_C_shieldings)
#     H_idx.append(conf_H_idx)
#     C_idx.append(conf_C_idx)

	# concatenate H and C data one after the other (for printing in the CSV)
	H_C_shieldings = np.concatenate((final_H_shieldings, final_C_shieldings))
	H_C_idx = np.concatenate((conf_H_idx, conf_C_idx))

	# create a CSV file containing all the data for the dr
	df = pd.DataFrame(list(zip(H_C_idx,H_C_shieldings)),
			   columns =['Atom', 'Shielding'])

	w_dir = os.getcwd()
	folder = w_dir + '/NMR_averaged'
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

def calculate_boltz_and_nmr(val,args,log,name,w_dir_fin):
	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name, val]
	subprocess.call(cmd_boltz)
	for lot_sp in args.level_of_theory_sp:
		for bs_sp in args.basis_set_sp:
			for bs_gcp_sp in args.basis_set_genecp_atoms_sp:
				#performing the nmr calculations
				dir_sp_nmr =  w_dir_fin+'/single_point_input_files/'+str(lot_sp)+'-'+str(bs_sp)
				os.chdir(dir_sp_nmr)
				#grabbing the respective NMR files for a given molecules
				nmr_log_files = glob.glob(name+'*_NMR.log') + glob.glob(name+'*_NMR_*.log')
				calculate_nmr(nmr_log_files,args,log,name,w_dir_fin,lot_sp,bs_sp)
