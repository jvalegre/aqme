#!/usr/bin/env python

#####################################################.
#            This file stores all the functions         #
#      used for genrating all paramets from cclib       #
#####################################################.

import subprocess, sys, os, math
from pyconfort.qprep_gaussian import moving_files
import pandas as pd
from pyconfort.argument_parser import possible_atoms
import json
import numpy as np

possible_atoms = possible_atoms()

def calculate_boltz_for_cclib(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	if str(bs).find('/') > -1:
		destination = w_dir_initial+'/QPRED/cclib-json/boltz/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		destination = w_dir_initial+'/QPRED/cclib-json/boltz/'+str(lot)+'-'+str(bs)

	moving_files(destination, os.getcwd(),'Goodvibes_'+name+'.dat')

def log_json(log_file,w_dir_initial,args,lot,bs):
	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/cclib-json/all_confs_cclib/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/cclib-json/all_confs_cclib/'+str(lot)+'-'+str(bs)
	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	file_out = folder +'/'+log_file.split('.log')[0]+'.json'
	cmd = ['ccframe',log_file,'-O',file_out]
	cmd_sub = " ".join(cmd)
	os.system(cmd_sub)

def calculate_cclib(log_files,args,log,name,w_dir,w_dir_initial,lot,bs):
	for file in log_files:
		log_json(file,w_dir_initial,args,lot,bs)

def calcualte_average_cclib_parameter(json_files,args,log,name,w_dir,w_dir_initial,lot,bs):
	charge,dipole = [],[]
	for file in json_files:
		mol = open(file, 'r')
		mol_data = json.load(mol)

		# get the Boltzmann probabilities of each conformer from a GoodVibes file
		if str(bs).find('/') > -1:
			boltz_file = w_dir_initial+'/QPRED/cclib-json/boltz/'+str(lot)+'-'+str(bs).split('/')[0] + '/Goodvibes_'+name+'.dat'
		else:
			boltz_file = w_dir_initial+'/QPRED/cclib-json/boltz/'+str(lot)+'-'+str(bs) + '/Goodvibes_'+name+'.dat'

		boltz_outfile = open(boltz_file,"r")
		boltz_outlines = boltz_outfile.readlines()
		for i in range(len(boltz_outlines)):
			# I remove the NMR from the file names using [0:-4]
			if boltz_outlines[i].find(file.split('.json')[0]) > -1:
				boltz_factor = float(boltz_outlines[i].split()[-1])
				break
		#charge
		charge.append(np.array(mol_data['atomcharges']['0']['mulliken']).astype(float)*boltz_factor)
		#dipole
		dipole.append(np.linalg.norm(np.array(mol_data['moments']['0'][1]))*boltz_factor)

	charge = np.sum(charge, axis=0).tolist()
	dipole = np.sum(dipole)

	dict_param = {'name': name,'atomcharges': charge, 'dipole': dipole}

	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/cclib-json/average_cclib/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/cclib-json/average_cclib/'+str(lot)+'-'+str(bs)

	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	with open(folder+'/'+name+'.json', "w") as outfile:
		json.dump(dict_param, outfile)
