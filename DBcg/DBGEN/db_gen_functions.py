"""

* In this file, the paths to helper programs are collected.
* You must make sure that all the variables are correct before launching db_gen.py.

* OTHER functions USED THROUGHOUT THE PROGRAM ARE ALSO SET HERE.

"""

# Since there are many imports, we could save some CPU time if we moved
# some of these import inside the functions that need them
import glob, os, shutil, sys, time
import numpy as np
import pandas as pd

from rdkit import Chem,DataStructs
from rdkit.Chem import PropertyMol, rdChemReactions,AllChem,Lipinski,Descriptors
import openbabel as ob
import subprocess

from confgen import *

from db_gen import columns,possible_atoms


" FUCNTION WORKING WITH MOL OBJECT TO CREATE CONFORMERS"
def conformer_generation(mol,name,args):
	valid_structure = filters(mol, args)
	#print(valid_structure)
	#print(valid_structure)
	if valid_structure:
		if args.verbose: print("\n   ----- {} -----".format(name))

		try:
			# the conformational search
			summ_search(mol, name,args)

			# the multiple minimization returns a list of separate mol objects
			conformers, energies = mult_min(mol, name, args)
			#print(energies)

			#only print if within predefined energy window
			if len(conformers) > 0:
				# list in energy order
				cids = list(range(len(conformers)))
				sortedcids = sorted(cids, key = lambda cid: energies[cid])

				np.set_printoptions(precision=3)
				#print(name, *sorted(energies), sep = ", ")
				print("o  Final energies:", np.array(sorted(energies)))
				#print(["{0:0.2f}".format(i) for i in sorted(energies)])

				sdwriter = Chem.SDWriter(name+final_output)
				glob_min = min(energies)

				write_confs = 0
				for cid in sortedcids:
					if args.ANI1ccx == True:
						if (energies[cid] - glob_min) < args.ewin / 2625.5:
							sdwriter.write(conformers[cid])
							write_confs += 1

					elif args.xtb == True:
						if (energies[cid] - glob_min) < args.ewin / 2625.5:
							sdwriter.write(conformers[cid])
							write_confs += 1

					else:
						if (energies[cid] - glob_min) < args.ewin:
							sdwriter.write(conformers[cid])
							write_confs += 1


				if args.verbose == True: print("   ----- {} conformers written to {} -----".format(write_confs, name+final_output))
				sdwriter.close()
			else: print("x  No conformers found!\n")

		except (KeyboaddInterrupt, SystemExit):
			raise
		except Exception as e: print(traceback.print_exc())
	else: print("ERROR: The structure is not valid")

	if args.time: print("Execution time: %s seconds" % (round(time.time() - start_time,2)))

" FILTER TO BE APPLIED FOR SMILES"
def filters(mol,args):
	valid_structure = True
	# First filter: number of rotatable bonds
	if Lipinski.NumRotatableBonds(mol) < args.max_torsions:
		# Second filter: molecular weight
		if Descriptors.MolWt(mol) < args.max_MolWt:
			# Third filter: this filters salts off (2 separated components)
			#if len(Chem.MolToSmiles(mol).split('.')) == 1:
			for atom in mol.GetAtoms():
				#Fourth filter: atoms outside the scope chosen in 'possible_atoms'
				if atom.GetSymbol() not in possible_atoms:
					valid_structure = False
					if args.verbose == True: print(" Exiting as atom isn't in atoms in the periodic table")
			#else: valid_structure = False
		else:
			valid_structure = False
			if args.verbose == True: print(" Exiting as total molar mass > 1000")
	else:
		valid_structure = False
		if args.verbose == True: print(" Exiting as number of rotatable bonds > 10")
	return valid_structure

def read_energies(file): # parses the energies from sdf files - then used to filter conformers
	energies = []
	f = open(file,"r")
	readlines = f.readlines()
	for i in range(len(readlines)):
		if readlines[i].find('>  <Energy>') > -1:
			energies.append(float(readlines[i+1].split()[0]))
	f.close()
	return energies

" MAIN FUNCTION TO CREATE GAUSSIAN JOBS"
def write_gaussian_input_file(file, name,lot, bs, bs_gcp, energies, args):

	#definition of input lines
	if args.frequencies == True:
		if args.dispersion_correction == True:
		    if args.solvent_model == 'gas_phase':
		        input = 'opt freq=noraman EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)
		        input_sp = 'nmr=giao EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)  #input for single point nmr
		    else :
		        input = 'opt freq=noraman SCRF=({0},Solvent={1}) EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name,args.empirical_dispersion ) #add solvent if needed
		        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name, args.empirical_dispersion)  ##add solvent if needed
		else:
		    if args.solvent_model == 'gas_phase':
		        input = 'opt freq=noraman '
		        input_sp = 'nmr=giao ' #input for single point nmr
		    else :
		        input = 'opt freq=noraman SCRF=({0},Solvent={1})'.format(args.solvent_model, args.solvent_name) #add solvent if needed
		        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao'.format(args.solvent_model, args.solvent_name)  ##add solvent if needed
	else:
		if args.dispersion_correction == True:
		    if args.solvent_model == 'gas_phase':
		        input = 'opt EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)
		        input_sp = 'nmr=giao EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)  #input for single point nmr
		    else :
		        input = 'opt SCRF=({0},Solvent={1}) EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name,args.empirical_dispersion ) #add solvent if needed
		        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name, args.empirical_dispersion)  ##add solvent if needed
		else:
		    if args.solvent_model == 'gas_phase':
		        input = 'opt '
		        input_sp = 'nmr=giao ' #input for single point nmr
		    else :
		        input = 'opt SCRF=({0},Solvent={1})'.format(args.solvent_model, args.solvent_name) #add solvent if needed
		        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao'.format(args.solvent_model, args.solvent_name)  ##add solvent if needed

	#defining genecp
	genecp = 'gen'

	#reading the sdf to check for I atom_symbol
	suppl = Chem.SDMolSupplier(file)
	for mol in suppl:
		for atom in mol.GetAtoms():
			if atom.GetSymbol() in args.genecp_atoms:
				genecp = 'genecp'

	if args.single_point == True:
		#pathto change to
		path_write_gjf_files = 'sp/' + str(lot) + '-' + str(bs)
		#print(path_write_gjf_files)
		os.chdir(path_write_gjf_files)
	else:
		#pathto change to
		path_write_gjf_files = 'gaussian/' + str(lot) + '-' + str(bs)
		#print(path_write_gjf_files)
		os.chdir(path_write_gjf_files)

	path_for_file = '../../'

	com = '{0}_.com'.format(name)
	com_low = '{0}_low.com'.format(name)

	if genecp =='genecp':
		#chk option
		if args.chk == True:
			if args.single_point == True:
				header = [
					'%chk={}.chk'.format(name),
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input_sp ]
			else:
				header = [
						'%chk={}.chk'.format(name),
						'%MEM={}'.format(args.mem),
						'%nprocshared={}'.format(args.nprocs),
						'# {0}'.format(lot)+ '/'+ genecp + ' '+ input ]

		else:
			if args.single_point == True:
				header = [
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input_sp ]
			else:
				header = [
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ genecp + ' '+ input ]

		if args.lowest_only == True:
			subprocess.run(
				  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com_low,'-l' , '1', '-xk', '\n'.join(header)]) #takes the lowest conformer which is the first in the file
		elif args.lowest_n == True:
			no_to_write = 0
			if len(energies) != 1:
				for i in range(len(energies)):
					energy_diff = energies[i] - energies[0]
					if energy_diff < args.energy_threshold_for_gaussian:
						no_to_write +=1
				subprocess.run(
					 ['obabel', '-isdf', path_for_file+file, '-f', '1', '-l' , str(no_to_write), '-osdf', '-Otemp.sdf'])
				subprocess.run(
					  ['obabel', '-isdf', 'temp.sdf', '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])
			else:
				subprocess.run(
					  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])
		else:
			subprocess.run(
				  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])

		#adding the basis set at the end of the FILES
		#grab all the com FILES
		com_files = glob.glob('{0}_*.com'.format(name))

		for file in com_files:
			ecp_list,ecp_genecp_atoms = [],False
			read_lines = open(file,"r").readlines()
			fileout = open(file, "a")
			# Detect if there are I atoms to use genecp or not (to use gen)
			for i in range(4,len(read_lines)):
				if read_lines[i].split(' ')[0] not in ecp_list and read_lines[i].split(' ')[0] in possible_atoms:
					ecp_list.append(read_lines[i].split(' ')[0])
				if read_lines[i].split(' ')[0] in args.genecp_atoms:
				   ecp_genecp_atoms = True

			for i in range(len(ecp_list)):
				if ecp_list[i] not in args.genecp_atoms:
					fileout.write(ecp_list[i]+' ')
			fileout.write('0\n')
			fileout.write(bs+'\n')
			fileout.write('****\n')
			if ecp_genecp_atoms == False:
				fileout.write('\n')
			else:
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n')
				fileout.write('****\n\n')
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n\n')
			fileout.close()

			#submitting the gaussian file on summit
			if args.qsub == True:
				os.system(submission_command + file)

		os.chdir(path_for_file)

	else:
		#chk option
		if args.chk == True:
			if args.single_point == True:
				header = [
					'%chk={}.chk'.format(name),
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input_sp ]
			else:
				header = [
						'%chk={}.chk'.format(name),
						'%MEM={}'.format(args.mem),
						'%nprocshared={}'.format(args.nprocs),
						'# {0}'.format(lot)+ '/'+ bs + ' '+ input ]

		else:
			if args.single_point == True:
				header = [
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input_sp ]
			else:
				header = [
					'%MEM={}'.format(args.mem),
					'%nprocshared={}'.format(args.nprocs),
					'# {0}'.format(lot)+ '/'+ bs + ' '+ input ]

		if args.lowest_only == True:
			subprocess.run(
				  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com_low,'-l' , '1', '-xk', '\n'.join(header)]) #takes the lowest conformer which is the first in the file
		elif args.lowest_n == True:
			no_to_write = 0
			if len(energies) != 1:
				for i in range(len(energies)):
					energy_diff = energies[i] - energies[0]
					if energy_diff < args.energy_threshold_for_gaussian:
						no_to_write +=1
				subprocess.run(
					 ['obabel', '-isdf', path_for_file+file, '-f', '1', '-l' , str(no_to_write), '-osdf', '-Otemp.sdf'])
				subprocess.run(
					  ['obabel', '-isdf', 'temp.sdf', '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])
			else:
				subprocess.run(
					  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])
		else:
			subprocess.run(
				  ['obabel', '-isdf', path_for_file+file, '-ocom', '-O'+com,'-m', '-xk', '\n'.join(header)])

" DEFINTION OF OUTPUT ANALYSER and NMR FILES CREATOR"
def output_analyzer(log_files, w_dir, lot, bs,bs_gcp, args):

	print(w_dir)

	#definition of input lines
	if args.dispersion_correction == True:
	    if args.solvent_model == 'gas_phase':
	        input = 'opt freq=noraman EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)
	        input_sp = 'nmr=giao EmpiricalDispersion=G{0}'.format(args.empirical_dispersion)  #input for single point nmr
	    else :
	        input = 'opt freq=noraman SCRF=({0},Solvent={1}) EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name,args.empirical_dispersion ) #add solvent if needed
	        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao EmpiricalDispersion=G{2}'.format(args.solvent_model, args.solvent_name, args.empirical_dispersion)  ##add solvent if needed
	else:
	    if args.solvent_model == 'gas_phase':
	        input = 'opt freq=noraman '
	        input_sp = 'nmr=giao ' #input for single point nmr
	    else :
	        input = 'opt freq=noraman SCRF=({0},Solvent={1})'.format(args.solvent_model, args.solvent_name) #add solvent if needed
	        input_sp = 'SCRF=({0},Solvent={1}) nmr=giao'.format(args.solvent_model, args.solvent_name)  ##add solvent if needed

	for file in log_files:

		#made it global for all functions
		rms = 10000
		#defined the variable stop_rms, standor
		stop_rms = 0
		standor = 0
		NATOMS = 0

		outfile = open(file,"r")
		outlines = outfile.readlines()
		ATOMTYPES, CARTESIANS = [],[]
		FREQS, REDMASS, FORCECONST, NORMALMODE = [],[],[],[]; IM_FREQS = 0
		freqs_so_far = 0
		TERMINATION = "unfinished"
		###stop=0
		### Change to reverse
		for i in range(0,len(outlines)):
			###if stop == 3: break
			# Get the name of the compound (specified in the title)
			if outlines[i].find('Symbolic Z-matrix:') > -1:
				name = outlines[i-2]
				###stop=stop+1
			# Determine the kind of job termination
			if outlines[i].find("Normal termination") > -1:
				TERMINATION = "normal"
				###stop=stop+1
			elif outlines[i].find("Error termination") > -1:
				TERMINATION = "error"
				###stop=stop+1
			# Determine charge and multiplicity
			if outlines[i].find("Charge = ") > -1:
				CHARGE = int(outlines[i].split()[2])
				MULT = int(outlines[i].split()[5].rstrip("\n"))
				###stop=stop+1


		###reverse
		for i in range(0,len(outlines)):
			if TERMINATION == "normal":
				# Sets where the final coordinates are inside the file
				###if outlines[i].find("Input orientation") > -1: standor = i
				if outlines[i].find("Standard orientation") > -1: standor = i
				if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1:
					if outlines[i-1].find("-------") > -1:
						NATOMS = i-standor-6
						###break
				# Get the frequencies and identifies negative frequencies
				if outlines[i].find(" Frequencies -- ") > -1:
					nfreqs = len(outlines[i].split())
					for j in range(2, nfreqs):
						FREQS.append(float(outlines[i].split()[j]))
						NORMALMODE.append([])
						if float(outlines[i].split()[j]) < 0.0: IM_FREQS += 1
					for j in range(3, nfreqs+1): REDMASS.append(float(outlines[i+1].split()[j]))
					for j in range(3, nfreqs+1): FORCECONST.append(float(outlines[i+2].split()[j]))
					for j in range(0,NATOMS):
						for k in range(0, nfreqs-2):
							NORMALMODE[(freqs_so_far + k)].append([float(outlines[i+5+j].split()[3*k+2]), float(outlines[i+5+j].split()[3*k+3]), float(outlines[i+5+j].split()[3*k+4])])
					freqs_so_far = freqs_so_far + nfreqs - 2
			if TERMINATION != "normal":
				if outlines[i].find('Cartesian Forces:  Max') > -1:
					if float(outlines[i].split()[5]) < rms:
						rms = float(outlines[i].split()[5])
						stop_rms = i

		if TERMINATION == "normal":
			# Get the coordinates for jobs that finished well with and without imag. freqs
			try: standor
			except NameError: pass
			else:
				for i in range (standor+5,standor+5+NATOMS):
					massno = int(outlines[i].split()[1])
					if massno < len(possible_atoms):
						atom_symbol = possible_atoms[massno]
					else: atom_symbol = "XX"
					ATOMTYPES.append(atom_symbol)
					CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

		if TERMINATION != "normal":
			# Get the coordinates for jobs that did not finished or finished with an error
			###reverse copy from above
			for i in range(0,stop_rms):
				# Sets where the final coordinates are inside the file
				if outlines[i].find("Input orientation") > -1: standor = i
				if outlines[i].find("Standard orientation") > -1: standor = i
				if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1:
					if outlines[i-1].find("-------") > -1:
						NATOMS = i-standor-6
			###no change after this
			for i in range (standor+5,standor+5+NATOMS):
				massno = int(outlines[i].split()[1])
				if massno < len(possible_atoms):
					atom_symbol = possible_atoms[massno]
				else: atom_symbol = "XX"
				ATOMTYPES.append(atom_symbol)
				CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

		# This part fixes jobs with imaginary freqs
		if IM_FREQS > 0:
			# Multiplies the imaginary normal mode vector by this amount (from -1 to 1).
			amplitude = 0.2 # default in pyQRC
			shift = []

			# Save the original Cartesian coordinates before they are altered
			orig_carts = []
			for atom in range(0,NATOMS):
				orig_carts.append([CARTESIANS[atom][0], CARTESIANS[atom][1], CARTESIANS[atom][2]])

			# could get rid of atomic units here, if zpe_rat definition is changed
			for mode, wn in enumerate(FREQS):
				# Either moves along any and all imaginary freqs, or a specific mode requested by the user
				if FREQS[mode] < 0.0:
					shift.append(amplitude)
				else: shift.append(0.0)

			# The starting geometry is displaced along the each normal mode according to the random shift
				for atom in range(0,NATOMS):
					for coord in range(0,3):
						CARTESIANS[atom][coord] = CARTESIANS[atom][coord] + NORMALMODE[mode][atom][coord] * shift[mode]
		outfile.close()

		# This part places the calculations in different folders depending on the type of
		# termination and number of imag. freqs
		source = w_dir+file

		if IM_FREQS == 0 and TERMINATION == "normal" and args.nmr != True:
			#only if normally terminated move to the finished folder of first run.
			if args.secondrun != True:
				destination = w_dir+'Finished/'
			else:
				destination = w_dir+'../Finished/'
			try:
				os.makedirs(destination)
				shutil.move(source, destination)
			except OSError:
				if  os.path.isdir(destination) and not os.path.exists(destination+file):
					shutil.move(source, destination)
				else:
					raise

		if IM_FREQS > 0:
			destination = w_dir+'Imaginary_frequencies/'
			try:
				os.makedirs(destination)
				shutil.move(source, destination)
			except OSError:
				if  os.path.isdir(destination) and not os.path.exists(destination+file):
					shutil.move(source, destination)
				else:
					raise

		if IM_FREQS == 0 and TERMINATION == "error":
			destination = w_dir+'Failed_Error/'
			try:
				os.makedirs(destination)
				shutil.move(source, destination)
			except OSError:
				if  os.path.isdir(destination) and not os.path.exists(destination+file):
					shutil.move(source, destination)
				else:
					raise

		if IM_FREQS == 0 and TERMINATION == "unfinished":
			destination = w_dir+'Failed_Unfinished/'
			try:
				os.makedirs(destination)
				shutil.move(source, destination)
			except OSError:
				if  os.path.isdir(destination) and not os.path.exists(destination+file):
					shutil.move(source, destination)
				else:
					raise

		if IM_FREQS > 0 or TERMINATION != "normal":
			# creating new folder with new input gaussian files
			new_gaussian_input_files = w_dir+'New Gaussian Input Files'

			try:
				os.makedirs(new_gaussian_input_files)
			except OSError:
				if  os.path.isdir(new_gaussian_input_files):
					os.chdir(new_gaussian_input_files)
				else:
					raise

			os.chdir(new_gaussian_input_files)
			print("vrbbtbtwrb")
			print('Creating new gaussian input files')

			# Options for genecp
			ecp_list,ecp_genecp_atoms = [],False

			for i in range(len(ATOMTYPES)):
				if ATOMTYPES[i] not in ecp_list and ATOMTYPES[i] in possible_atoms:
					ecp_list.append(ATOMTYPES[i])
				if ATOMTYPES[i] in args.genecp_atoms:
				   ecp_genecp_atoms = True
			if ecp_genecp_atoms == False:
				genecp = 'gen'
			if ecp_genecp_atoms == True:
				genecp = 'genecp'

			if args.single_point == True:
				keywords_opt = lot +'/'+ genecp+' '+ input_sp
			else:
				keywords_opt = lot +'/'+ genecp+' '+ input

			fileout = open(file.split(".")[0]+'.com', "w")
			fileout.write("%mem="+str(args.mem)+"\n")
			fileout.write("%nprocshared="+str(args.nprocs)+"\n")
			fileout.write("# "+keywords_opt+"\n")
			fileout.write("\n")
			fileout.write(name+"\n")
			fileout.write("\n")
			fileout.write(str(CHARGE)+' '+str(MULT)+'\n')
			for atom in range(0,NATOMS):
				fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(ATOMTYPES[atom], CARTESIANS[atom][0],  CARTESIANS[atom][1],  CARTESIANS[atom][2]))
				fileout.write("\n")
			fileout.write("\n")
			for i in range(len(ecp_list)):
				if ecp_list[i] not in args.genecp_atoms:
					fileout.write(ecp_list[i]+' ')
			fileout.write('0\n')
			fileout.write(bs+'\n')
			fileout.write('****\n')
			if ecp_genecp_atoms == False:
				fileout.write('\n')
			else:
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n')
				fileout.write('****\n\n')
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n\n')
			fileout.close()

		#changing directory back to where all files are from new files created.
		os.chdir(w_dir)

		#adding in the NMR componenet only to the finished files after reading from normally finished log files
		if args.nmr == True:

			# creating new folder with new input gaussian files
			nmr_gaussian_input_files = w_dir+'/NMR_Gaussian_input_files'

			try:
				os.makedirs(nmr_gaussian_input_files)
			except OSError:
				if  os.path.isdir(nmr_gaussian_input_files):
					os.chdir(nmr_gaussian_input_files)
				else:
					raise

			os.chdir(nmr_gaussian_input_files)
			print('Creating new NMR files')

			# Options for genecp
			ecp_list,ecp_genecp_atoms = [],False

			for i in range(len(ATOMTYPES)):
				if ATOMTYPES[i] not in ecp_list and ATOMTYPES[i] in possible_atoms:
					ecp_list.append(ATOMTYPES[i])
				if ATOMTYPES[i] in args.genecp_atoms:
				   ecp_genecp_atoms = True
			if ecp_genecp_atoms == False:
				genecp = 'gen'
			if ecp_genecp_atoms == True:
				genecp = 'genecp'

			keywords_opt = lot +'/'+ genecp+' '+ ' nmr=giao'

			fileout = open(file.split(".")[0]+'.com', "w")
			fileout.write("%mem="+str(args.mem)+"\n")
			fileout.write("%nprocshared="+str(args.nprocs)+"\n")
			fileout.write("# "+keywords_opt+"\n")
			fileout.write("\n")
			fileout.write(name+"\n")
			fileout.write("\n")
			fileout.write(str(CHARGE)+' '+str(MULT)+'\n')
			for atom in range(0,NATOMS):
				fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(ATOMTYPES[atom], CARTESIANS[atom][0],  CARTESIANS[atom][1],  CARTESIANS[atom][2]))
				fileout.write("\n")
			fileout.write("\n")
			for i in range(len(ecp_list)):
				if ecp_list[i] not in args.genecp_atoms:
					fileout.write(ecp_list[i]+' ')
			fileout.write('0\n')
			fileout.write(bs+'\n')
			fileout.write('****\n')
			if ecp_genecp_atoms == False:
				fileout.write('\n')
			else:
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n')
				fileout.write('****\n\n')
				for i in range(len(ecp_list)):
					if ecp_list[i] in args.genecp_atoms:
						fileout.write(ecp_list[i]+' ')
				fileout.write('0\n')
				fileout.write(bs_gcp+'\n\n')
			fileout.close()

		#changing directory back to where all files are from new files created.
		os.chdir(w_dir)

" CALCULATION OF BOLTZMANN FACTORS "
def boltz_calculation(val,i):
	#need to have good vibes
	cmd = 'python' +  ' -m' + ' goodvibes' + ' --csv' + ' --boltz ' +'--output ' + str(i) + ' ' + val
	os.system(cmd)

"COMBINING FILES FOR DIFFERENT MOLECULES"
def combine_files(csv_files, lot, bs, args):
	#final dataframe with only the boltzmann averaged values
	final_file_avg_thermo_data = pd.DataFrame(columns=columns)
	compare_G = pd.DataFrame(columns=['Structure_of_min_conf','min_qh-G(T)','boltz_avg_qh-G(T)'])

	files = []
	#combine all the csv_files

	for f in csv_files:

		df = pd.read_csv(f, skiprows = 14)
		# df['Structure']= df['Structure'].astype(str)
		df = df.rename(columns={"   Structure": "Structure"})

		#dropping the ************* line
		df = df.drop(df.index[0])
		df.iloc[-1] = np.nan

		for col in columns:
			if col == 'Structure':
				#identifyin the minmum energy if the conformers
				min_G = df['qh-G(T)'].min()
				#getting the name of the structure of the min G
				idx_name_of_min_conf = df['qh-G(T)'].idxmin() - 1
				name_of_min_conf = df.iloc[idx_name_of_min_conf]['Structure']
				#df.at[df.index[-1], col] = name_of_min_conf
			elif col != 'Structure':
				boltz_avg = np.sum(df[col] * df['Boltz'])
				df.at[df.index[-1], col] = boltz_avg
				if col == 'qh-G(T)':
					compare_G = compare_G.append({'Structure_of_min_conf': name_of_min_conf,'min_qh-G(T)': min_G,'boltz_avg_qh-G(T)': boltz_avg}, ignore_index=True)

		final_file_avg_thermo_data = final_file_avg_thermo_data.append({'Structure':name_of_min_conf , 'E': df.iloc[-1]['E'] , 'ZPE': df.iloc[-1]['ZPE'], 'H':df.iloc[-1]['H'] , 'T.S':df.iloc[-1]['T.S'] , 'T.qh-S':df.iloc[-1]['T.qh-S'] , 'G(T)': df.iloc[-1]['G(T)'], 'qh-G(T)':df.iloc[-1]['qh-G(T)'] },ignore_index=True)

		files.append(df)

	final_file_all_data = pd.concat(files, axis=0, ignore_index=True)

	#combined_csv = pd.concat([pd.read_csv(f, skiprows = 14, skipfooter = 1) for f in csv_files ])
	#change directory to write all files in one place
	destination = args.path+'All csv files/'+ str(lot)+ '-'+ str(bs)
	try:
		os.makedirs(destination)
	except OSError:
		if  os.path.isdir(destination):
			pass
		else:
			raise
	os.chdir(destination)

	#export to csv
	final_file_all_data.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_all data.csv', index=False, encoding='utf-8-sig')
	final_file_avg_thermo_data.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_avg_thermo_data.csv', index=False, encoding='utf-8-sig')
	compare_G.to_csv( str(lot) + '-' + str(bs) + '_all_molecules_compare_G(T).csv', index=False, encoding='utf-8-sig')
