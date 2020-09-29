#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	 used for genrating all parameters from DBSTEP  #
#####################################################.

import os
import numpy as np
import pandas as pd
import subprocess
from pyconfort.qprep_gaussian import moving_files

def calculate_db_parameters(log_files,args,log,w_dir_initial,name_mol,lot,bs):

	try:
		from dbstep.Dbstep import dbstep
	except (ModuleNotFoundError,AttributeError):
		log.write('\nx  DBSTEP is not installed correctly - DBSTEP is not available')
		sys.exit()

	total_data = pd.DataFrame()

	#find the ring atoms in the File
	filelines =  open(w_dir_initial+'/'+args.dbstep_cen_lig_file,'r').readlines()

	for counter,log in enumerate(log_files):
		for line in (filelines):
			split_line = line.strip().split(',')
			if split_line[0] == name_mol:
				C = split_line[1]
				L =  split_line[2]
				break
		sterics = dbstep(log, atom1=str(C),atom2=str(L), volume=True, sterimol=True, commandline=True)
		total_data.at[counter,'Name'] = name_mol
		total_data.at[counter,'log'] = log.split('.log')[0]
		total_data.at[counter,'bv'] = sterics.bur_vol
		total_data.at[counter,'bmax'] = sterics.Bmax
		total_data.at[counter,'bmin'] = sterics.Bmin
		total_data.at[counter,'L'] = sterics.L

	#creating folder for all molecules to write geom parameter
	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/dbstep_parameters/all_confs_sterics/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/dbstep_parameters/all_confs_sterics/'+str(lot)+'-'+str(bs)

	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	total_data.to_csv(folder+'/'+name_mol+'-all-steric-data.csv',index=False)

def calculate_boltz_and_dbstep(val,args,log,name,w_dir,w_dir_initial,lot,bs):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	if str(bs).find('/') > -1:
		destination = w_dir_initial+'/QPRED/dbstep_parameters/boltz/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		destination = w_dir_initial+'/QPRED/dbstep_parameters/boltz/'+str(lot)+'-'+str(bs)

	moving_files(destination, os.getcwd(),'Goodvibes_'+name+'.dat')

	if str(bs).find('/') > -1:
		dbstep_parm_file = w_dir_initial + '/QPRED/dbstep_parameters/all_confs_sterics/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+name+'-all-steric-data.csv'
	else:
		dbstep_parm_file = w_dir_initial + '/QPRED/dbstep_parameters/all_confs_sterics/'+str(lot)+'-'+str(bs)+'/'+name+'-all-steric-data.csv'


	df_dbstep =  pd.read_csv(dbstep_parm_file)

	if str(bs).find('/') > -1:
		file = w_dir_initial+'/QPRED/dbstep_parameters/boltz/'+str(lot)+'-'+str(bs).split('/')[0]+'/Goodvibes_'+name+'.dat'
	else:
		file = w_dir_initial+'/QPRED/dbstep_parameters/boltz/'+str(lot)+'-'+str(bs)+'/Goodvibes_'+name+'.dat'

	outlines = open(file,"r").readlines()
	#reading the data from boltz fileyiu
	for i in range(len(outlines)):
		# I remove the NMR from the file names using [0:-4]
		if outlines[i].find('   ***************************************************************************************************************************************\n') > -1 and outlines[i-1].find('   Structure') > -1:
			start_line = i+1
		elif outlines[i].find('   ***************************************************************************************************************************************\n') > -1:
			end_line = i

	boltz_values = pd.DataFrame()
	counter = 0
	for i in range(start_line,end_line):
		values = outlines[i].split()
		boltz_values.at[counter,'Name'] = values[1]
		boltz_values.at[counter,'boltz'] = values[-1]
		counter +=1

	# multiply actual file and sum and add write to new FILE
	df_dbstep['bv'] = df_dbstep['bv'].astype(float)*boltz_values['boltz'].astype(float)
	df_dbstep['bmin'] = df_dbstep['bmin'].astype(float)*boltz_values['boltz'].astype(float)
	df_dbstep['bmax'] = df_dbstep['bmax'].astype(float)*boltz_values['boltz'].astype(float)
	df_dbstep['L'] = df_dbstep['L'].astype(float)*boltz_values['boltz'].astype(float)

	avg_data = pd.DataFrame()
	avg_data.at[0,'Name'] = name
	avg_data.at[0,'bv'] = df_dbstep.sum().bv
	avg_data.at[0,'bmin'] = df_dbstep.sum().bmin
	avg_data.at[0,'bmax'] = df_dbstep.sum().bmax
	avg_data.at[0,'L'] = df_dbstep.sum().L

	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/dbstep_parameters/average_sterics/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/dbstep_parameters/average_sterics/'+str(lot)+'-'+str(bs)

	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	avg_data.to_csv(folder+'/'+name+'-steric-data.csv',index=False)
