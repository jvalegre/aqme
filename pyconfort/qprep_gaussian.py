#!/usr/bin/env python

#########################################################.
# 		   This file stores all the functions 		    #
# 	       used in writing SDF and COM files,           #
#              as well as the logger and                #
#                 yaml file importer		            #
#########################################################.

import os
import sys
import subprocess
import glob
import shutil
import pandas as pd
from rdkit.Chem import AllChem as Chem
from pyconfort.argument_parser import possible_atoms


possible_atoms = possible_atoms()

def convert_xyz_to_sdf(xyz_files,args,log):
	for file in xyz_files:
		name=file.split('.xyz')[0]
		command_xyz = ['obabel', '-ixyz', file, '-osdf', '-O'+name+'.sdf']
		subprocess.call(command_xyz)

def header_com(name,lot,bs,bs_gcp, args, log, input_route, genecp):
	if genecp != 'None':
		input_route_to_write = input_route
		genecp_or_bs_to_write = genecp
	else:
		input_route_to_write = input_route
		genecp_or_bs_to_write = bs
	#chk option
	if args.chk:
		header = [
			'%chk={}.chk'.format(name),
			'%mem={}'.format(args.mem),
			'%nprocshared={}'.format(args.nprocs),
			'# {0}'.format(lot)+ '/'+ genecp_or_bs_to_write + ' '+ input_route_to_write ]
	else:
		header = [
			'%mem={}'.format(args.mem),
			'%nprocshared={}'.format(args.nprocs),
			'# {0}'.format(lot)+ '/'+ genecp_or_bs_to_write + ' '+ input_route_to_write]

	return header

def convert_sdf_to_com(path_for_file,file,com,com_low,energies,header,args,log):

	if args.lowest_only and args.lowest_n:
		log.write('x  Both lowest \'n\' and lowest are turned on. Writing lowest \'n\'')
		args.lowest_only = False

	if args.lowest_only:
		command_lowest = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com_low,'-l' , '1', '-xk', '\n'.join(header)]
		subprocess.call(command_lowest) #takes the lowest conformer which is the first in the file

	elif args.lowest_n:
		no_to_write = 0
		if len(energies) != 1:
			for i,_ in enumerate(energies):
				energy_diff = energies[i] - energies[0]
				if energy_diff < args.energy_threshold_for_gaussian: # thershold is in kcal/mol and energies are in kcal/mol as well
					no_to_write +=1
			command_n = ['obabel', '-isdf', path_for_file+file, '-f', '1', '-l' , str(no_to_write), '-osdf', '-Otemp.sdf']
			subprocess.call(command_n)
			command_n_2 =  ['obabel', '-isdf', 'temp.sdf', '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.call(command_n_2)
			os.remove('temp.sdf')
		else:
			command_n_3 = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.call(command_n_3)
	else:
		command_no_lowest = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
		subprocess.call(command_no_lowest)

def input_route_line(args):
	#definition of input_route lines
	if args.input_for_gauss == 'None':
		input_route = ''
		if args.frequencies:
			input_route += 'freq=noraman'
		if args.empirical_dispersion != 'None':
			input_route += ' empiricaldispersion={0}'.format(args.empirical_dispersion)
		if not args.QCORR:
			input_route += ' opt=(maxcycles={0})'.format(args.max_cycle_opt)
		else:
			input_route += ' opt=(calcfc,maxcycles={0})'.format(args.max_cycle_opt)
		if args.solvent_model != 'gas_phase':
			input_route += ' scrf=({0},solvent={1})'.format(args.solvent_model,args.solvent_name)
	else:
		input_route = args.input_for_gauss

	return input_route

def rename_file_and_charge_chk_change(read_lines,file,args,charge_com):
	#changing the name of the files to the way they are in xTB Sdfs
	#getting the title line
	if not args.com_from_xyz:
		for i,line in enumerate(read_lines):
			if len(line.strip()) == 0:
				title_line = read_lines[i+1]
				title_line = title_line.lstrip()
				rename_file_name = title_line.replace(" ", "_")
				break

		rename_file_name = rename_file_name.strip()+'.com'
	else:
		rename_file_name = file

	#change charge and multiplicity for all molecules
	for i,_ in enumerate(read_lines):
		if args.chk:
			if read_lines[i].find('%chk') > -1:
				read_lines[i] = '%chk='+rename_file_name.split('.com')[0]+'.chk\n'
		if len(read_lines[i].strip()) == 0:
			if  args.com_from_xyz:
				read_lines[i+1] = rename_file_name.split('.com')[0]+'\n'
			if charge_com is not None:
				read_lines[i+3] = str(charge_com)+' '+ str(args.complex_spin)+'\n'
			break
	out = open(file, 'w')
	out.writelines(read_lines)
	out.close()
	return rename_file_name

def get_name_and_charge(name,charge_data):

	if charge_data is not None:
		name_list = name.split('_')

		if 'rules' in name_list:
			name_molecule = name[:-23]
		elif 'xtb' in name_list or 'ani' in name_list:
			name_molecule = name[:-4]
		elif 'rotated' in name_list:
			name_molecule = name[:-14]
		elif 'rdkit' in name_list:
			name_molecule = name[:-6]
		for i in range(len(charge_data)):
			if charge_data.loc[i,'Molecule'] == name_molecule:
				charge_com = charge_data.loc[i,'Overall charge']
	else:
		charge_com = None

	return charge_com

def get_genecp(file,args):
	genecp = 'None'

	try:
		#reading the sdf to check for I atom_symbol
		suppl = Chem.SDMolSupplier(file)
		for atom in suppl[0].GetAtoms():
			if atom.GetSymbol() in args.genecp_atoms:
				genecp = 'genecp'
				break
			elif atom.GetSymbol() in args.gen_atoms:
				genecp = 'gen'
				break
	except OSError:
		outfile = open(file,"r")
		outlines = outfile.readlines()
		for _,line in enumerate(outlines):
			if args.genecp_atoms:
				for atom in args.genecp_atoms:
					if line.find(atom)>-1:
						genecp = 'genecp'
						break
			elif args.gen_atoms:
				for atom in args.gen_atoms:
					if line.find(atom)>-1:
						genecp = 'gen'
						break
		outfile.close()

	return genecp

# MAIN FUNCTION TO CREATE GAUSSIAN JOBS
def write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args, log, charge_data=None):

	# get the names of the SDF files to read from depending on the optimizer and their suffixes. Also, get molecular charge
	charge_com = get_name_and_charge(name,charge_data)


	input_route = input_route_line(args)

	#defining genecp
	genecp = get_genecp(file,args)

	# defining path to place the new COM files
	if args.single_point:
		path_write_gjf_files = 'QMCALC/G16-SP/' + str(lot) + '-' + str(bs)
	else:
		path_write_gjf_files = 'QMCALC/G16/' + str(lot) + '-' + str(bs)

	os.chdir(path_write_gjf_files)

	path_for_file = '../../../'

	com = '{0}_.com'.format(name)
	com_low = '{0}_low.com'.format(name)

	header = header_com(name,lot, bs, bs_gcp,args,log, input_route, genecp)

	convert_sdf_to_com(path_for_file,file,com,com_low,energies,header,args,log)

	com_files = glob.glob('{0}_*.com'.format(name))

	for file in com_files:
		if genecp =='genecp' or genecp == 'gen':
			ecp_list,ecp_genecp_atoms,ecp_gen_atoms = [],False,False
			read_lines = open(file,"r").readlines()

			rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

			read_lines = open(file,"r").readlines()

			fileout = open(file, "a")
			# Detect if there are atoms to use genecp or not (to use gen)
			for i in range(4,len(read_lines)):
				if read_lines[i].split(' ')[0] not in ecp_list and read_lines[i].split(' ')[0] in possible_atoms:
					ecp_list.append(read_lines[i].split(' ')[0])
				if read_lines[i].split(' ')[0] in args.genecp_atoms:
				   ecp_genecp_atoms = True
				if read_lines[i].split(' ')[0] in args.gen_atoms:
				   ecp_gen_atoms = True

			#error if both genecp and gen are
			if ecp_genecp_atoms and ecp_gen_atoms:
				sys.exit("ERROR: Can't use Gen and GenECP at the same time")

			for _,element_ecp in enumerate(ecp_list):
				if element_ecp not in (args.genecp_atoms or args.gen_atoms):
					fileout.write(element_ecp+' ')
			fileout.write('0\n')
			fileout.write(bs+'\n')
			fileout.write('****\n')
			if not ecp_genecp_atoms and not ecp_gen_atoms:
				fileout.write('\n')
			else:
				format_ext_genecp = ['.txt','.yaml','.yml','.rtf']
				if len(bs_gcp.split('.')) > 1:
					if bs_gcp.split('.')[1] in format_ext_genecp:
						os.chdir(path_for_file)
						read_lines = open(bs_gcp,"r").readlines()
						os.chdir(path_write_gjf_files)
						#chaanging the name of the files to the way they are in xTB Sdfs
						#getting the title line
						for line in read_lines:
							fileout.write(line)
						fileout.write('\n\n')
				else:
					for _,element_ecp in enumerate(ecp_list):
						if element_ecp in args.genecp_atoms :
							fileout.write(element_ecp+' ')
						elif element_ecp in args.gen_atoms :
							fileout.write(element_ecp+' ')
					fileout.write('0\n')
					fileout.write(bs_gcp+'\n')
					fileout.write('****\n\n')
					if ecp_genecp_atoms:
						for _,element_ecp in enumerate(ecp_list):
							if element_ecp in args.genecp_atoms:
								fileout.write(element_ecp+' ')
						fileout.write('0\n')
						fileout.write(bs_gcp+'\n\n')
			fileout.close()

		else:
			read_lines = open(file,"r").readlines()

			rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

		#change file by moving to new file
		os.rename(file,rename_file_name)

		# #submitting the gaussian file on summit
		if args.qsub:
			cmd_qsub = [args.submission_command, rename_file_name]
			subprocess.call(cmd_qsub)

	os.chdir(path_for_file)

# MOVES SDF FILES TO THEIR CORRESPONDING FOLDERS
def moving_files(destination,src,file):
	try:
		os.makedirs(destination)
		shutil.move(os.path.join(src, file), os.path.join(destination, file))
	except OSError:
		if  os.path.isdir(destination):
			shutil.move(os.path.join(src, file), os.path.join(destination, file))
		else:
			raise

# WRITE SDF FILES FOR xTB AND ANI1
def write_confs(conformers, energies,selectedcids, name, args, program,log):
	if len(conformers) > 0:
		name = name.split('_rdkit')[0]# a bit hacky
		sdwriter = Chem.SDWriter(name+'_'+program+args.output)

		write_confs = 0
		for cid in selectedcids:
			sdwriter.write(conformers[cid])
			write_confs += 1

		if args.verbose:
			log.write("o  Writing "+str(write_confs)+ " conformers to file " + name+'_'+program+args.output)
		sdwriter.close()
	else:
		log.write("x  No conformers found!")

# PARSES THE ENERGIES FROM SDF FILES
def read_energies(file,log): # parses the energies from sdf files - then used to filter conformers
	energies = []
	f = open(file,"r")
	readlines = f.readlines()
	for i,_ in enumerate(readlines):
		if readlines[i].find('>  <Energy>') > -1:
			energies.append(float(readlines[i+1].split()[0]))
	f.close()
	return energies