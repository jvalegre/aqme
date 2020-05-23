#!/usr/bin/env python

"""#####################################################.
# 		   This file stores all the functions 		    #
# 	       used in writing SDF and COM files,           #
#              as well as the logger and                #
#                 yaml file importer		            #
######################################################"""

import os
import sys
import subprocess
import glob
import shutil
import yaml
import pandas as pd
from rdkit.Chem import AllChem as Chem
from pyconfort.argument_parser import possible_atoms

possible_atoms = possible_atoms()

# CLASS FOR LOGGING
class Logger:
	"""
	 Class Logger to write the output to a file
	"""
	def __init__(self, filein, append):
		"""
		Logger to write the output to a file
		"""
		suffix = 'dat'
		self.log = open('{0}_{1}.{2}'.format(filein, append, suffix), 'w')

	def write(self, message):
		print(message, end='\n')
		self.log.write(message+ "\n")

	def fatal(self, message):
		print(message, end='\n')
		self.log.write(message + "\n")
		self.finalize()
		sys.exit(1)

	def finalize(self):
		self.log.close()

# LOAD PARAMETERS FROM A YAML FILE
def load_from_yaml(args,log):
	# Variables will be updated from YAML file
	if args.varfile is not None:
		if os.path.exists(args.varfile):
			if os.path.splitext(args.varfile)[1] == '.yaml':
				log.write("\no  IMPORTING VARIABLES FROM " + args.varfile)
				with open(args.varfile, 'r') as file:
					param_list = yaml.load(file, Loader=yaml.FullLoader)
	else:
		log.write("\no  No yaml file containing parameters was found (the program might crush unless you specify the input files manually!).\n")
	for param in param_list:
		if hasattr(args, param):
			if getattr(args, param) != param_list[param]:
				log.write("o  RESET " + param + " from " + str(getattr(args, param)) + " to " + str(param_list[param]))
				setattr(args, param, param_list[param])
			else:
				log.write("o  DEFAULT " + param + " : " + str(getattr(args, param)))

def creation_of_dup_csv(args):
	# writing the list of DUPLICATES
	if args.nodihedrals:
		if not args.xtb and not args.ANI1ccx:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window', 'RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','time (seconds)','Overall charge'])
		elif args.xtb:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
		elif args.ANI1ccx:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])
	else:
		if not args.xtb and not args.ANI1ccx:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','time (seconds)','Overall charge'])
		elif args.xtb:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-window','RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','xTB-Initial-samples','xTB-energy-window','xTB-initial_energy_threshold','xTB-RMSD-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
		elif args.ANI1ccx:
			dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples','RDKit-energy-window', 'RDKit-energy-duplicates','RDKit-RMSD-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-energy-window','ANI1ccx-initial_energy_threshold','ANI1ccx-RMSD-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])
	return dup_data

def header_com(name,lot,bs,bs_gcp, args, log, input_sp, input, genecp):
	#chk option
	if args.chk:
		if args.single_point:
			if genecp != 'None':
				header = [
					'%chk={}.chk'.format(name),
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input_sp ]
			else:
				header = [
					'%chk={}.chk'.format(name),
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input_sp ]
		else:
			if genecp != 'None':
				header = [
						'%chk={}.chk'.format(name),
						'%mem={}'.format(args.mem),
						'%nprocshared={}'.format(args.nprocs),
						'# {0}'.format(lot)+ '/'+ genecp + ' '+ input ]
			else:
				header = [
						'%chk={}.chk'.format(name),
						'%mem={}'.format(args.mem),
						'%nprocshared={}'.format(args.nprocs),
						'# {0}'.format(lot)+ '/'+ bs + ' '+ input ]
	else:
		if args.single_point:
			if genecp != 'None':
				header = [
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input_sp ]
			else:
				header = [
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input_sp ]

		else:
			if genecp != 'None':
				header = [
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input ]
			else:
				header = [
					'%mem={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input ]
	return header

def convert_sdf_to_com(path_for_file,file,com,com_low,energies,header,args):

	if args.lowest_only:
		command_lowest = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com_low,'-l' , '1', '-xk', '\n'.join(header)]
		subprocess.run(command_lowest) #takes the lowest conformer which is the first in the file

	elif args.lowest_n:
		no_to_write = 0
		if len(energies) != 1:
			for i,_ in enumerate(energies):
				energy_diff = energies[i] - energies[0]
				if energy_diff < args.energy_threshold_for_gaussian: # thershold is in kcal/mol and energies are in kcal/mol as well
					no_to_write +=1
			command_n = ['obabel', '-isdf', path_for_file+file, '-f', '1', '-l' , str(no_to_write), '-osdf', '-Otemp.sdf']
			subprocess.run(command_n)
			command_n_2 =  ['obabel', '-isdf', 'temp.sdf', '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.run(command_n_2)
		else:
			command_n_3 = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.run(command_n_3)
	else:
		command_no_lowest = ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
		subprocess.run(command_no_lowest)

def input_line(args):
	#definition of input lines
	if args.frequencies:
		if args.dispersion_correction:
			if args.solvent_model == 'gas_phase':
				input = 'opt=(maxcycles={0}) freq=noraman empiricaldispersion={1}'.format(args.max_cycle_opt,args.empirical_dispersion)
				input_sp = 'nmr=giao empiricaldispersion={0}'.format(args.empirical_dispersion)  #input for single point nmr
			else :
				input = 'opt=(maxcycles={0}) freq=noraman scrf=({1},solvent={2}) empiricaldispersion={3}'.format(args.max_cycle_opt, args.solvent_model, args.solvent_name,args.empirical_dispersion ) #add solvent if needed
				input_sp = 'scrf=({0},solvent={1}) nmr=giao empiricaldispersion={2}'.format(args.solvent_model, args.solvent_name, args.empirical_dispersion)  ##add solvent if needed
		else:
			if args.solvent_model == 'gas_phase':
				input = 'opt=(maxcycles={0}) freq=noraman'.format(args.max_cycle_opt)
				input_sp = 'nmr=giao ' #input for single point nmr
			else :
				input = 'opt=(maxcycles={0}) freq=noraman scrf=({1},solvent={2})'.format(args.max_cycle_opt,args.solvent_model, args.solvent_name) #add solvent if needed
				input_sp = 'scrf=({0},solvent={1}) nmr=giao'.format(args.solvent_model, args.solvent_name)  ##add solvent if needed
	else:
		if args.dispersion_correction:
			if args.solvent_model == 'gas_phase':
				input = 'opt=(maxcycles={0}) empiricaldispersion={1}'.format(args.max_cycle_opt,args.empirical_dispersion)
				input_sp = 'nmr=giao empiricaldispersion={0}'.format(args.empirical_dispersion)  #input for single point nmr
			else :
				input = 'opt=(maxcycles={0}) scrf=({1},solvent={2}) empiricaldispersion={3}'.format(args.max_cycle_opt,args.solvent_model, args.solvent_name,args.empirical_dispersion ) #add solvent if needed
				input_sp = 'scrf=({0},solvent={1}) nmr=giao empiricaldispersion={2}'.format(args.solvent_model, args.solvent_name, args.empirical_dispersion)  ##add solvent if needed
		else:
			if args.solvent_model == 'gas_phase':
				input = 'opt=(maxcycles={0})'.format(args.max_cycle_opt)
				input_sp = 'nmr=giao ' #input for single point nmr
			else :
				input = 'opt=(maxcycles={0}) scrf=({1},solvent={2})'.format(args.max_cycle_opt,args.solvent_model, args.solvent_name) #add solvent if needed
				input_sp = 'scrf=({0},solvent={1}) nmr=giao'.format(args.solvent_model, args.solvent_name)  ##add solvent if needed
	return input, input_sp

# MAIN FUNCTION TO CREATE GAUSSIAN JOBS
def write_gaussian_input_file(file, name,lot, bs, bs_gcp, energies, args,log,charge_data):

	#find location of molecule and respective scharges
	name_list = name.split('_')
	if 'xtb' or 'ani' in name_list:
		name_molecule = name[:-4]
	if 'rdkit' in name_list:
		name_molecule = name[:-6]
	if 'rotated' in name_list:
		name_molecule = name[:-14]

	for i in range(len(charge_data)):
		if charge_data.loc[i,'Molecule'] == name_molecule:
			charge_com = charge_data.loc[i,'Overall charge']

	input, input_sp = input_line(args)

	#defining genecp
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
	except:
		read_lines = open(file,"r").readlines()
		for line,_ in enumerate(read_lines):
			for atom in args.genecp_atoms:
				if read_lines[line].find(atom)>-1:
					genecp = 'genecp'
					break
			for atom in args.gen_atoms:
				if read_lines[line].find(atom)>-1:
					genecp = 'gen'
					break

	if args.single_point:
		#pathto change to
		path_write_gjf_files = 'generated_sp_files/' + str(lot) + '-' + str(bs)
		#log.write(path_write_gjf_files)
		os.chdir(path_write_gjf_files)
	else:
		#path to change to
		path_write_gjf_files = 'generated_gaussian_files/' + str(lot) + '-' + str(bs)
		os.chdir(path_write_gjf_files)

	path_for_file = '../../'

	com = '{0}_.com'.format(name)
	com_low = '{0}_low.com'.format(name)

	header = header_com(name,lot, bs, bs_gcp,args,log,input_sp, input, genecp)

	convert_sdf_to_com(path_for_file,file,com,com_low,energies,header,args)

	com_files = glob.glob('{0}_*.com'.format(name))

	for file in com_files:
		if genecp =='genecp' or genecp == 'gen':
			ecp_list,ecp_genecp_atoms,ecp_gen_atoms = [],False,False
			read_lines = open(file,"r").readlines()

			#chaanging the name of the files to the way they are in xTB Sdfs
			#getting the title line
			for i,line in enumerate(read_lines):
				if len(line.strip()) == 0:
					title_line = read_lines[i+1]
					title_line = title_line.lstrip()
					rename_file_name = title_line.replace(" ", "_")
					break

			rename_file_name = rename_file_name.strip()+'.com'

			#change charge and multiplicity for Octahydrasl
			if args.metal_complex:
				for i,line in enumerate(read_lines):
					if len(line.strip()) == 0:
						read_lines[i+3] = str(charge_com)+' '+ str(args.complex_spin)+'\n'
						break
				out = open(file, 'w')
				out.writelines(read_lines)
				out.close()
				read_lines = open(file,"r").readlines()

			fileout = open(file, "a")
			# Detect if there are I atoms to use genecp or not (to use gen)
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

			for i,element_ecp in enumerate(ecp_list):
				if element_ecp not in (args.genecp_atoms or args.gen_atoms):
					fileout.write(element_ecp+' ')
			fileout.write('0\n')
			fileout.write(bs+'\n')
			fileout.write('****\n')
			if not ecp_genecp_atoms and not ecp_gen_atoms:
				fileout.write('\n')
			else:
				if len(bs_gcp.split('.')) > 1:
					if bs_gcp.split('.')[1] == 'txt' or bs_gcp.split('.')[1] == 'yaml':
						os.chdir(path_for_file)
						read_lines = open(bs_gcp,"r").readlines()
						os.chdir(path_write_gjf_files)
						#chaanging the name of the files to the way they are in xTB Sdfs
						#getting the title line
						for line in read_lines:
							fileout.write(line)
						fileout.write('\n\n')
				else:
					for i,element_ecp in enumerate(ecp_list):
						if element_ecp in args.genecp_atoms :
							fileout.write(element_ecp+' ')
						elif element_ecp in args.gen_atoms :
							fileout.write(element_ecp+' ')
					fileout.write('0\n')
					fileout.write(bs_gcp+'\n')
					fileout.write('****\n\n')
					if ecp_genecp_atoms:
						for i,element_ecp in enumerate(ecp_list):
							if element_ecp in args.genecp_atoms:
								fileout.write(element_ecp+' ')
						fileout.write('0\n')
						fileout.write(bs_gcp+'\n\n')
			fileout.close()

		else:
			read_lines = open(file,"r").readlines()

			#changing the name of the files to the way they are in xTB Sdfs
			#getting the title line
			for i,line in enumerate(read_lines):
				if len(line.strip()) == 0:
					title_line = read_lines[i+1]
					title_line = title_line.lstrip()
					rename_file_name = title_line.replace(" ", "_")
					break

			rename_file_name = rename_file_name.strip()+'.com'

			#change charge and multiplicity for Octahydrasl
			if args.metal_complex:
				for i,_ in enumerate(read_lines):
					if len(read_lines[i].strip()) == 0:
						read_lines[i+3] = str(charge_com)+' '+ str(args.complex_spin)+'\n'
						break
				out = open(file, 'w')
				out.writelines(read_lines)
				out.close()

		#change file by moving to new file
		os.rename(file,rename_file_name)

		#submitting the gaussian file on summit
		if args.qsub:
			cmd_qsub = [args.submission_command, rename_file_name]
			subprocess.run(cmd_qsub)

	os.chdir(path_for_file)

# MOVES SDF FILES TO THEIR CORRESPONDING FOLDERS
def moving_sdf_files(destination,src,file):
	try:
		os.makedirs(destination)
		shutil.move(os.path.join(src, file), os.path.join(destination, file))
	except OSError:
		if  os.path.isdir(destination):
			shutil.move(os.path.join(src, file), os.path.join(destination, file))
		else:
			raise

# WRITE SDF FILES FOR xTB AND ANI1
def write_confs(conformers, energies, name, args, program,log):
	if len(conformers) > 0:
		# list in energy order
		cids = list(range(len(conformers)))
		sortedcids = sorted(cids, key = lambda cid: energies[cid])

		name = name.split('_rdkit')[0]# a bit hacky
		sdwriter = Chem.SDWriter(name+'_'+program+args.output)

		write_confs = 0
		for cid in sortedcids:
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

def write_gauss_main(args,log):
	if args.exp_rules:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif not args.xtb and not args.ANI1ccx and args.nodihedrals:
			conf_files =  glob.glob('*_rdkit.sdf')
	elif not args.xtb and not args.ANI1ccx and not args.nodihedrals:
		conf_files =  glob.glob('*_rdkit_rotated.sdf')
	elif args.xtb:
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.ANI1ccx:
		conf_files =  glob.glob('*_ani.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	# names for directories created
	sp_dir = 'generated_sp_files'
	g_dir = 'generated_gaussian_files'

	#read in dup_data to get the overall charge of MOLECULES
	charge_data = pd.read_csv(args.input.split('.')[0]+'-Duplicates Data.csv', usecols=['Molecule','Overall charge'])

	for lot in args.level_of_theory:
		for bs in args.basis_set:
			for bs_gcp in args.basis_set_genecp_atoms:
				# only create this directory if single point calculation is requested
				if args.single_point:
					folder = sp_dir + '/' + str(lot) + '-' + str(bs)
					log.write("\no  PREPARING SINGLE POINT INPUTS in {}".format(folder))
				else:
					folder = g_dir + '/' + str(lot) + '-' + str(bs)
					log.write("\no  Preparing Gaussian COM files in {}".format(folder))
				try:
					os.makedirs(folder)
				except OSError:
					if os.path.isdir(folder):
						pass
					else:
						raise
				# writing the com files
				# check conf_file exists, parse energies and then write dft input
				for file in conf_files:
					if os.path.exists(file):
						if args.verbose:
							log.write("   -> Converting from {}".format(file))
						energies = read_energies(file,log)
						name = os.path.splitext(file)[0]

						write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args,log,charge_data)
