#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	     used for genrating details for energy      #
#####################################################.

from pyconfort.qcorr_gaussian import moving_files

from rdkit.Chem import AllChem as Chem
import numpy as np
import pandas as pd
import os
import subprocess

def calculate_boltz_and_energy(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):
	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	if str(bs).find('/') > -1:
		destination = w_dir_initial+'/QPRED/energy/boltz/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		destination = w_dir_initial+'/QPRED/energy/boltz/'+str(lot)+'-'+str(bs)
	moving_files('Goodvibes_'+name+'.dat', destination)

def calculate_avg_and_energy(val,args,log,name,w_dir_fin,w_dir_initial,w_dir_boltz,lot,bs):

	all_file_data = []
	for file in val:
		outlines = open(file,"r").readlines()

		for i in range(len(outlines)):
			# I remove the NMR from the file names using [0:-4]
			if outlines[i].find('   ***************************************************************************************************************************************\n') > -1 and outlines[i-1].find('   Structure') > -1:
				start_line = i+1
			elif outlines[i].find('   ***************************************************************************************************************************************\n') > -1:
				end_line = i

		one_file_data,values_avg= [],[]
		for i in range(start_line,end_line):
			values = outlines[i].split()
			for j in range(2,9):
				values_avg.append(float(values[j])*float(values[-1]))
			one_file_data.append(values_avg)
			values_avg = []

		res = []
		res.append(name)
		for j in range(0, len(one_file_data[0])):
			tmp = 0
			for i in range(0, len(one_file_data)):
				tmp = tmp + one_file_data[i][j]
			res.append(tmp)
		all_file_data.append(res)
	if str(bs).find('/') > -1:
		folder = w_dir_initial+'/QPRED/energy/average-energy/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial+'/QPRED/energy/average-energy/'+str(lot)+'-'+str(bs)
	try:
		os.makedirs(folder)
		os.chdir(folder)
	except OSError:
		if os.path.isdir(folder):
			os.chdir(folder)
		else:
			raise
	all_file_data_df = pd.DataFrame(all_file_data,columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)'])
	if str(bs).find('/') > -1:
		all_file_data_df.to_csv(w_dir_initial+'/QPRED/energy/average-energy/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+args.input.split('.')[0]+'.csv', index=False)
	else:
		all_file_data_df.to_csv(w_dir_initial+'/QPRED/energy/average-energy/'+str(lot)+'-'+str(bs)+'/'+args.input.split('.')[0]+'.csv', index=False)
