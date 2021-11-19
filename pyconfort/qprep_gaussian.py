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
import shlex
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
		genecp_or_bs_to_write = genecp
	else:
		genecp_or_bs_to_write = bs
	#chk option
	if args.chk:
		header = [
			'%chk={}.chk'.format(name),
			'%mem={}'.format(args.mem),
			'%nprocshared={}'.format(args.nprocs),
			'# {0}'.format(lot)+ '/'+ genecp_or_bs_to_write + ' '+ input_route]
	else:
		header = [
			'%mem={}'.format(args.mem),
			'%nprocshared={}'.format(args.nprocs),
			'# {0}'.format(lot)+ '/'+ genecp_or_bs_to_write + ' '+ input_route]

	return header

def convert_sdf_to_com(w_dir_initial,file,com,com_low,energies,header,args,log):

	if args.lowest_only and args.lowest_n:
		log.write('x  The lowest_n and lowest_only options are both True, lowest_n will be used')
		args.lowest_only = False

	if args.lowest_only:
		command_lowest = ['obabel', '-isdf', w_dir_initial+'/'+file, '-ocom', '-O'+com_low,'-l' , '1', '-xk', '\n'.join(header)]
		subprocess.call(command_lowest) #takes the lowest conformer which is the first in the file
		log.write('o  The lowest_only option is activated (only using the lowest energy conformer)')

	elif args.lowest_n:
		log.write('o  The lowest_n option is True (only using conformers within the specified E window)')
		no_to_write = 0
		if len(energies) != 1:
			for i,_ in enumerate(energies):
				energy_diff = energies[i] - energies[0]
				if energy_diff < args.energy_threshold_for_gaussian: # thershold is in kcal/mol and energies are in kcal/mol as well
					no_to_write +=1
			command_n = ['obabel', '-isdf', w_dir_initial+'/'+file, '-f', '1', '-l' , str(no_to_write), '-osdf', '-Otemp.sdf']
			subprocess.call(command_n)
			command_n_2 =  ['obabel', '-isdf', 'temp.sdf', '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.call(command_n_2)
			os.remove('temp.sdf')
		else:
			command_n_3 = ['obabel', '-isdf', w_dir_initial+'/'+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
			subprocess.call(command_n_3)

	else:
		command_no_lowest = ['obabel', '-isdf', w_dir_initial+'/'+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)]
		subprocess.call(command_no_lowest)

def input_route_line(args):
	#definition of input_route lines
	if args.set_input_line == 'None':
		input_route = ''
		if args.QPREP == 'gaussian' or args.QCORR == 'gaussian':
			if args.frequencies:
				input_route += 'freq=noraman'
			if args.empirical_dispersion != 'None':
				input_route += ' empiricaldispersion={0}'.format(args.empirical_dispersion)
			if not args.calcfc:
				input_route += ' opt=(maxcycles={0})'.format(args.max_cycle_opt)
			else:
				input_route += ' opt=(calcfc,maxcycles={0})'.format(args.max_cycle_opt)
			if args.solvent_model != 'gas_phase':
				input_route += ' scrf=({0},solvent={1})'.format(args.solvent_model,args.solvent_name)
	else:
		input_route = args.set_input_line

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

		if args.QPREP == 'gaussian':
			rename_file_name = rename_file_name.strip()+'.com'
		elif args.QPREP == 'orca':
			rename_file_name = rename_file_name.strip()+'.inp'
	else:
		rename_file_name = file

	if args.QPREP == 'gaussian':
		#change charge and multiplicity for all molecules
		for i,_ in enumerate(read_lines):
			if args.chk:
				if read_lines[i].find('%chk') > -1:
					read_lines[i] = '%chk='+rename_file_name.split('.com')[0]+'.chk\n'
			if len(read_lines[i].strip()) == 0:
				if  args.com_from_xyz:
					read_lines[i+1] = rename_file_name.split('.com')[0]+'\n'
				if charge_com is not None:
					read_lines[i+3] = str(charge_com)+' '+ str(args.mult)+'\n'
				break
		out = open(file, 'w')
		out.writelines(read_lines)
		out.close()

	return rename_file_name

def get_name_and_charge(name,charge_data):

	name_list = name.split('_')

	if 'xtb' in name_list or 'ani' in name_list:
		if 'filter' in name_list:
			name_molecule = name[:-21]
		else:
			name_molecule = name[:-4]
	elif 'summ' in name_list:
		if 'filter' in name_list:
			name_molecule = name[:-22]
		else:
			name_molecule = name[:-5]
	elif 'rdkit' in name_list:
		if 'filter' in name_list:
			name_molecule = name[:-23]
		else:
			name_molecule = name[:-6]
	elif 'fullmonte' in name_list:
		if 'filter' in name_list:
			name_molecule = name[:-27]
		else:
			name_molecule = name[:-10]
	if charge_data is not None:
		for i in range(len(charge_data)):
			if charge_data.loc[i,'Molecule'] == name_molecule:
				charge_com = charge_data.loc[i,'Overall charge']
			else:
				try:
					suppl = Chem.SDMolSupplier(name+'.sdf', removeHs=False)
				except OSError:
					suppl = False
				if suppl:
					mol = suppl[0]
					charge_com = mol.GetProp('Real charge')
				else:
					charge_com = 'Invalid'

		return charge_com

	else:
		return name_molecule

# DETECTION AND LISTING OF GEN/GENECP FROM COM FILES
def check_for_gen_or_genecp(ATOMTYPES,args,type_of_check,program_gen):
	# Options for genecp
	ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = [],False,False,None,False
	if type_of_check == 'analysis':
		genecp_atoms_include = args.genecp_atoms
		gen_atoms_include = args.gen_atoms
		aux_atoms_include = args.aux_atoms_orca

	elif type_of_check == 'sp':
		genecp_atoms_include = args.genecp_atoms_sp
		gen_atoms_include = args.gen_atoms_sp
		aux_atoms_include = args.aux_atoms_orca_sp

	for _,atomtype in enumerate(ATOMTYPES):
		if program_gen == 'orca':
			if atomtype in aux_atoms_include:
				orca_aux_section = True

		if program_gen == 'gaussian':
			if atomtype not in ecp_list and atomtype in possible_atoms:
				ecp_list.append(atomtype)
			if atomtype in genecp_atoms_include:
				ecp_genecp_atoms = True
			if atomtype in gen_atoms_include:
				ecp_gen_atoms = True

	if ecp_gen_atoms:
		genecp = 'gen'
	if ecp_genecp_atoms:
		genecp = 'genecp'

	return ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section

# DETECTION OF GEN/GENECP FROM SDF FILES
def get_genecp(args):
	genecp = 'None'

	if len(args.genecp_atoms) > 0:
		genecp = 'genecp'

	if len(args.gen_atoms) > 0:
		genecp = 'gen'

	return genecp

# write genecp/gen part
def write_genecp(ATOMTYPES,type_gen,fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,new_gaussian_input_files):

	ecp_not_used = 0
	for _,element_ecp in enumerate(ecp_list):
		if element_ecp in ATOMTYPES and element_ecp != '':
			if type_gen == 'sp':
				if element_ecp not in (args.genecp_atoms_sp or args.gen_atoms_sp):
					fileout.write(element_ecp+' ')
					ecp_not_used += 1
			else:
				if element_ecp not in (args.genecp_atoms or args.gen_atoms):
					fileout.write(element_ecp+' ')
					ecp_not_used += 1

	if ecp_not_used > 0:
		fileout.write('0\n')
		fileout.write(bs_com+'\n')
		fileout.write('****\n')

	if len(bs_gcp_com.split('.')) > 1:
		if bs_gcp_com.split('.')[1] == ('txt' or 'yaml' or 'yml' or 'rtf'):
			os.chdir(w_dir_initial)
			read_lines = open(bs_gcp_com,"r").readlines()
			os.chdir(new_gaussian_input_files)
			#getting the title line
			for line in read_lines:
				fileout.write(line)
			fileout.write('\n\n')
	else:
		ecp_used = 0
		for _,element_ecp in enumerate(ecp_list):
			if type_gen == 'sp':
				if element_ecp in (args.genecp_atoms_sp or args.gen_atoms_sp):
					fileout.write(element_ecp+' ')
					ecp_used += 1
			else:
				if element_ecp in (args.genecp_atoms or args.gen_atoms):
					fileout.write(element_ecp+' ')
					ecp_used += 1

		if ecp_used > 0:
			fileout.write('0\n')
			fileout.write(bs_gcp_com+'\n')
			fileout.write('****\n\n')

		else:
			fileout.write('\n')

		if ecp_genecp_atoms:
			for _,element_ecp in enumerate(ecp_list):
				if type_gen == 'sp':
					if element_ecp in args.genecp_atoms_sp:
						fileout.write(element_ecp+' ')
				else:
					if element_ecp in args.genecp_atoms:
						fileout.write(element_ecp+' ')

			fileout.write('0\n')
			fileout.write(bs_gcp_com+'\n\n')

def orca_file_gen(read_lines,rename_file_name,bs,lot,genecp,ecp_list,bs_gcp,bs_gcp_fit,charge_com,mult_com,orca_aux_section,args,extra_input,solvation,solvent_orca,cpcm_input_orca,scf_iters_orca,orca_mdci,print_mini):
	try:
		write_orca_lines = open(rename_file_name,"w")

	except FileExistsError:
		os.remove(rename_file_name)
		write_orca_lines = open(rename_file_name,"w")

	write_orca_lines.write("# "+rename_file_name.split('.')[0]+"\n")
	write_orca_lines.write("# Memory per core\n")
	# calculate memory for ORCA input
	if args.mem.find('GB') > -1 or args.mem.find('gb') > -1:
		mem_orca = int(args.mem[:-2])*1000
	elif args.mem.find('MB') > -1 or args.mem.find('mb') > -1:
		mem_orca = int(args.mem[:-2])

	write_orca_lines.write("%maxcore "+ str(mem_orca) + "\n")
	write_orca_lines.write("# Number of processors\n")
	write_orca_lines.write("%pal nprocs " + str(args.nprocs) +" end\n")

	if extra_input != 'None':
		write_orca_lines.write("! " + bs + " " + lot + " " + extra_input +"\n")
	else:
		write_orca_lines.write("! " + bs + " " + lot + "\n")
	if orca_aux_section:
		write_orca_lines.write("%basis\n")
		gen_orca_line_1 = 'NewGTO '
		gen_orca_line_2 = 'NewAuxCGTO '
		for gen_orca_atom in ecp_list:
			gen_orca_line_1 += str(possible_atoms.index(gen_orca_atom))+' '
			gen_orca_line_2 += str(possible_atoms.index(gen_orca_atom))+' '
		if len(bs_gcp) != 0:
			gen_orca_line_1 += '"' + bs_gcp[0] + '" end\n'
		if len(bs_gcp_fit) != 0:
			gen_orca_line_2 += '"' + bs_gcp_fit[0] + '" end\n'
		write_orca_lines.write(gen_orca_line_1)
		write_orca_lines.write(gen_orca_line_2)
		write_orca_lines.write("end\n")
	if solvation != 'gas_phase':
		if solvation.lower() == 'cpcm':
			write_orca_lines.write("! CPCM("+solvent_orca+")\n")
		if cpcm_input_orca != 'None' or solvation.lower() == 'smd':
			write_orca_lines.write("%cpcm\n")
			if cpcm_input_orca != 'None':
				for cpcm_line in cpcm_input_orca:
					write_orca_lines.write(cpcm_line+'\n')
			if solvation.lower() == 'smd':
				write_orca_lines.write("smd true\n")
				write_orca_lines.write('SMDsolvent "'+solvent_orca+'"\n')
			write_orca_lines.write("end\n")
	write_orca_lines.write("%scf maxiter "+str(scf_iters_orca)+"\n")
	write_orca_lines.write("end\n")
	if orca_mdci != 'None':
		write_orca_lines.write("% mdci\n")
		for mdci_line in orca_mdci:
			write_orca_lines.write(mdci_line+'\n')
		write_orca_lines.write("end\n")
	if print_mini:
		write_orca_lines.write("% output\n")
		write_orca_lines.write("printlevel mini\n")
		write_orca_lines.write("print[ P_SCFInfo ] 1\n")
		write_orca_lines.write("print[ P_SCFIterInfo ] 1\n")
		write_orca_lines.write("print[ P_OrbEn ] 0\n")
		write_orca_lines.write("print[ P_Cartesian ] 0\n")
		write_orca_lines.write("end\n")
		write_orca_lines.write("% elprop\n")
		write_orca_lines.write("Dipole False\n")
		write_orca_lines.write("end\n")
	write_orca_lines.write("* xyz "+str(charge_com)+" "+str(mult_com)+"\n")
	# this part gets the coordinates from above
	start_coords = 10000
	for i,line in enumerate(read_lines):
		if len(line.split()) > 0:
			# this if statement gets where the coordinates start
			if line.split()[0] == '#':
				start_coords = i+5
		if i >= start_coords:
			if len(line.split()) == 0:
				break
			write_orca_lines.write(line)
	write_orca_lines.write("*")
	write_orca_lines.close()

# MAIN FUNCTION TO CREATE GAUSSIAN JOBS
def write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args, log, charge_data, w_dir_initial):

	# get the names of the SDF files to read from depending on the optimizer and their suffixes. Also, get molecular charge
	charge_com = get_name_and_charge(name,charge_data)

	if charge_com != 'Invalid':
		input_route = input_route_line(args)

		#defining genecp
		genecp = get_genecp(file,args)

		# defining path to place the new COM files
		if args.QPREP == 'gaussian':
			if str(bs).find('/') > -1:
				path_write_input_files = '/QMCALC/G16/' + str(lot) + '-' + str(bs).split('/')[0]
			else:
				path_write_input_files = '/QMCALC/G16/' + str(lot) + '-' + str(bs)
		elif args.QPREP == 'orca':
			if str(bs).find('/') > -1:
				path_write_input_files = '/QMCALC/ORCA/' + str(lot) + '-' + str(bs).split('/')[0]
			else:
				path_write_input_files = '/QMCALC/ORCA/' + str(lot) + '-' + str(bs)

		os.chdir(w_dir_initial+path_write_input_files)

		try:
			com = '{0}_.com'.format(name)
			com_low = '{0}_low.com'.format(name)

			header = header_com(name,lot, bs, bs_gcp,args,log, input_route, genecp)

			convert_sdf_to_com(w_dir_initial,file,com,com_low,energies,header,args,log)

			com_files = glob.glob('{0}_*.com'.format(name))

			for file in com_files:
				#patch for Isotopes
				cmd = shlex.split(r"sed -i 's/\([a-zA-Z]\{1,3\}\)[^(]*\([(]Iso\)/\1\2/g' "+file)
				subprocess.call(cmd)
				ecp_list,ecp_genecp_atoms,ecp_gen_atoms = [],False,False
				read_lines = open(file,"r").readlines()
				rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)
				read_lines = open(file,"r").readlines()

				# Detect if there are atoms to use genecp or not (to use gen)
				ATOMTYPES = []
				for i in range(4,len(read_lines)):
					if read_lines[i].split(' ')[0] not in ATOMTYPES and read_lines[i].split(' ')[0] in possible_atoms:
						ATOMTYPES.append(read_lines[i].split(' ')[0])

				if genecp =='genecp' or genecp == 'gen' or args.last_line_for_input != 'None':
					# write genecp/gen part
					type_gen = 'qprep'
					if args.QPREP == 'gaussian':
						# define genecp/gen atoms
						ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')

						#error if both genecp and gen are
						if ecp_genecp_atoms and ecp_gen_atoms:
							sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")
						fileout = open(file, "a")

						write_genecp(type_gen,fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs,lot,bs_gcp,args,w_dir_initial,path_write_input_files)

						if args.last_line_for_input != 'None':
							fileout.write(args.last_line_for_input+'\n\n')

						fileout.close()

						read_lines = open(file,"r").readlines()
						rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

				else:
					read_lines = open(file,"r").readlines()
					rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

				#change file by moving to new file
				try:
					os.rename(file,rename_file_name)

				except FileExistsError:
					os.remove(rename_file_name)
					os.rename(file,rename_file_name)

				if args.QPREP == 'orca':

					# define auxiliary atoms
					ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','orca')

					rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

					#create input file
					orca_file_gen(read_lines,rename_file_name,bs,lot,genecp,args.aux_atoms_orca,args.aux_basis_set_genecp_atoms,args.aux_fit_genecp_atoms,charge_com,args.mult,orca_aux_section,args,args.set_input_line,args.solvent_model,args.solvent_name,args.cpcm_input,args.orca_scf_iters,args.mdci_orca,args.print_mini_orca)

				# submitting the input file on a HPC
				if args.qsub:
					cmd_qsub = [args.submission_command, rename_file_name]
					subprocess.call(cmd_qsub)

		except OSError:
			pass

	os.chdir(w_dir_initial)

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
		# name = name.split('_'+args.CSEARCH)[0]# a bit hacky
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
