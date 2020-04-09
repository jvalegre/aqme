'''
#ADD description
DBcg v1.0.0
Authors: Shree Sowndarya S. V., Juan V. Alegre Requena,
please report any bugs to svss@colostate.edu or juanvi89@hotmail.com.
'''

from DBGEN.db_gen_functions import *

def main():
	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")

	#INput details
	parser.add_argument("--varfile",dest="varfile",default=None,help="Parameters in python format")
	parser.add_argument("--input",help="Molecular structure")

	#metal complex
	parser.add_argument("--metal_complex", action="store_true", default=False, help="If studying Metal Complexes of coordination number 4,5 or 6")
	parser.add_argument("--metal",  help="If metal complex, specify the metal involed", default="Ir", dest="metal", type=str)
	parser.add_argument("--complex_spin",  help="If metal complex, specify the complex spin involed", default="1", dest="complex_spin", type=int)
	parser.add_argument("--complex_coord",  help="If metal complex, specify the complex coordination involed", default="6", dest="complex_coord", type=int)
	parser.add_argument("--complex_type",  help="If metal complex, specify the complex type involed", default="octahedral", dest="complex_type", type=str)
	parser.add_argument("--m_oxi",  help="If metal complex, specify the metal oxidation state involed", default="3", dest="m_oxi", type=int)
	parser.add_argument("--charge",  help="If metal complex,charge of metal complex. Will automatically update for metals", default="0", dest="charge", type=int)

	#NCI complex
	parser.add_argument("--nci_complex", action="store_true", default=False, help="If studying NCI Complexes")

	parser.add_argument("-m","--maxnumber", help="Number of compounds", type=int, metavar="maxnumber")
	parser.add_argument("--prefix", help="Prefix for naming files", default=None, metavar="prefix")

	#work the script has to do
	parser.add_argument("-w","--compute", action="store_true", default=False, help="Create conformers")
	parser.add_argument("--write_gauss", action="store_true", default=False, help="Create input files for Gaussian")
	parser.add_argument("-a","--analysis", action="store_true", default=False, help="Fix and analyze Gaussian output files")
	parser.add_argument("-r","--resubmit", action="store_true", default=False, help="Resubmit Gaussian input files")


	#Post analysis
	parser.add_argument("--dup",action="store_true",default=False, help="Remove Duplicates after DFT optimization")
	parser.add_argument("-n","--nmr",action="store_true",default=False, help="Create Files for Single Point which includes NMR calculation after DFT Optimization")
	parser.add_argument("-b","--boltz", action="store_true", default=False, help="Boltzmann factor for each conformers from Gaussian output files")
	parser.add_argument("-f","--combine", action="store_true", default=False, help="Combine files of differnt molecules including boltzmann weighted energies")

	#aaply exp rules
	parser.add_argument("--exp_rules", action="store_true", default=False, help="Experimental rules applied to make Gaussian input files")
	parser.add_argument("--angle_off", type=float,help="Any limit to set for check rules",default=30)

	#pass the argument for path for the gaussian folder.
	parser.add_argument("--path", help="Path for analysis/boltzmann factor/combining files where the gaussian folder created is present")
	parser.add_argument("-v","--verbose",action="store_true",default=False, help="verbose output")

	#argumets for conformer generation
	parser.add_argument("--ANI1ccx", "--ani", action="store_true",default=False, help="request ANI1ccx optimizations")
	parser.add_argument("--xtb", action="store_true",default=False, help="request xtb optimizations")
	parser.add_argument("--ewin", action="store",default=40.0, help="energy window to print conformers (kJ/mol)", type=float)
	parser.add_argument("--convergence","-c", action="store",default=1.0, help="Adjust convergence criteria of ANI and xtb optimizations (set at 0.005)", type=float)
	parser.add_argument("--time","-t",action='store_true',default=False,help="request program runtime")
	parser.add_argument("--heavyonly", help="only consider torsion angles involving heavy (non H) elements (default=True)", default=True, metavar="heavyonly")
	parser.add_argument("--constraints", action="store_true",default=None, help="distance constraint")
	parser.add_argument("--nodihedrals", action="store_true", default=False, help="turn off dihedral scan")
	parser.add_argument("-d","--degree", type=float,help="Amount, in degrees, to enumerate torsions by (default 30.0)",default=30.0)
	parser.add_argument("--etkdg", dest="etkdg",action="store_true",default=False, help="use new ETKDG knowledge-based method instead of distance geometry")
	parser.add_argument("--max_torsions",type=int,help="Skip any molecules with more than this many torsions (default 5)",default=5)
	parser.add_argument("--sample", help="number of conformers to sample to get non-torsional differences (default 100)", default=100, type=int, metavar="sample")
	parser.add_argument("--ff", help="force field (MMFF or UFF)", default="MMFF", metavar="ff")
	parser.add_argument("--seed", help="random seed (default 062609)", default="062609", type=int, metavar="s")
	parser.add_argument("--rms_threshold", help="cutoff for considering sampled conformers the same (default 0.25)", default=0.25, type=float, metavar="R")
	parser.add_argument("--energy_threshold", dest="energy_threshold",action="store",default=0.05, help="energy difference between unique conformers")
	parser.add_argument("--max_MolWt", help="Max. molecular weight of molecule", default=1000, type=int, metavar="max_MolWt")
	parser.add_argument("--large_sys", action="store_true",default=False, help="Large systems for xtb optimizations")
	parser.add_argument("--STACKSIZE", help="STACKSIZE for optimization of large systems", default="500m")

	#arguments for gaussian files Creation
	parser.add_argument("-l", "--level_of_theory",help="Level of Theory", default=['wB97xd'], dest="level_of_theory", type=str, nargs='*')
	parser.add_argument("--basis_set",  help="Basis Set", default=['6-31g*'], dest="basis_set", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms",default=['LANL2DZ'], help="Basis Set genecp: The length has to be the same as basis_set", dest="basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--genecp_atoms",  help="genecp atoms",default=[], dest="genecp_atoms",type=str, nargs='*')
	parser.add_argument("--max_cycle_opt", help="Number of cycles for DFT optimization", default="300", type=int, dest="max_cycle_opt")
	parser.add_argument("--frequencies",action="store_true", default=False, help="Request only optimization without any frequency calculation")
	parser.add_argument("--single_point",action="store_true", default=False, help="Request only single point calculation")
	parser.add_argument("--lowest_only", action="store_true", default=False, help="Lowest conformer to write for gaussian")
	parser.add_argument("--lowest_n", action="store_true", default=False, help="Lowest Number of conformers to write for gaussian")
	parser.add_argument("--energy_threshold_for_gaussian", help="cutoff for considering sampled conformers for gaussian input", default="4.0", type=float, dest="energy_threshold_for_gaussian")
	parser.add_argument("--dispersion_correction",action="store_true", default=False, help="Add Dispersion Correction")
	parser.add_argument("--empirical_dispersion",  help="Type of Dispersion ", default="D3BJ", dest="empirical_dispersion", type=str)
	parser.add_argument("--solvent_model",  help="Type of solvent model", default="gas_phase", dest="solvent_model", type=str)
	parser.add_argument("--solvent_name",  help="Name of Solvent", default="Acetonitrile", dest="solvent_name", type=str)
	parser.add_argument("--nprocs", help="Number of Processors", default="24", type=int, dest="nprocs")
	parser.add_argument("--mem", help="Memory", default="96GB", type=str, dest="mem")
	parser.add_argument("--chk", action="store_true", default=False, help="Create .chk files for Gaussian")

	#submmsion of gaussion files
	parser.add_argument("--qsub", action="store_true", default=False, help="Submit Gaussian files")
	parser.add_argument("--submission_command",  help="Queueing system that the submission is done on", default="qsub_summit", metavar="submission_command", type=str)

	args = parser.parse_args()

	var = False

	if args.varfile != None: var = True
	if args.time: start_time = time.time()

	### If the input file is python format then we will assume it contains variables ... ###
	if var == True:
		if os.path.splitext(args.varfile)[1] == '.py':
			db_gen_variables = os.path.splitext(args.varfile)[0]
			print("\no  IMPORTING VARIABLES FROM", args.varfile)
			args = __import__(db_gen_variables)

### check if it's working with *.smi (if it isn't, we need to change this a little bit)
### no it isnt working for *.smi, same or .sdf. We can add that feature in the end.

	if args.compute == True: # this will perform conformational analysis and create inputs for Gaussian

		# input file format specified
		[file_name, file_format] = os.path.splitext(args.input)

		if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com','.gjf']:
			print("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!"); sys.exit()
		else:
			mol_objects = [] # a list of mol objects that will be populated

		# first check if octahydral of not

		if file_format == '.smi': # SMILES input specified
			smifile = open(args.input)

			for i, line in enumerate(smifile):
				toks = line.split()

				if args.metal_complex == True and  args.complex_coord == 6:
					#find metal and replace with I+ for octahydral
					smi = toks[0].replace(args.metal,'I+')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I++','I+')

					#taking largest component for salts
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

				elif args.metal_complex == True and args.complex_coord == 5:
					#find metal and replace with I+ for octahydral
					smi = toks[0].replace(args.metal,'I')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I+','I')

				elif args.metal_complex == True and args.complex_coord == 4:
					#find metal and replace with I+ for octahydral
					smi = toks[0].replace(args.metal,'I-')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I-+','I-')

					#taking largest component for salts
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

					#taking largest component for salts
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

				# elif args.metal_complex == True and  args.complex_coord == 7:
				# 	#find metal and replace with I+ for octahydral
				# 	smi = toks[0].replace(args.metal,'I+2')
				# 	#Ir+ exixted then we need to change back to I+
				# 	smi = smi.replace('I+2+','I+2')
				#
				# elif args.metal_complex == True and  args.complex_coord == 8:
				# 	#find metal and replace with I+ for octahydral
				# 	smi = toks[0].replace(args.metal,'I+3')
				# 	#Ir+ exixted then we need to change back to I+
				# 	smi = smi.replace('I+3+','I+3')

				else: smi = toks[0]

				if args.prefix == None: name = ''.join(toks[1:])
				else: name = args.prefix+str(i)+'_'+''.join(toks[1:])

					# Converts each line to a rdkit mol object
				if args.verbose: print("   -> Input Molecule {} is {}".format(i, smi))
				mol = Chem.MolFromSmiles(smi)
				# get manually for square planar and SQUAREPYRIMIDAL
				if complex_type == 'squareplanar' or complex_type == 'squarepyrimidal':
					file_template = 'template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])

		elif os.path.splitext(args.input)[1] == '.csv': # CSV file with one columns SMILES and code_name
			csv_smiles = pd.read_csv(args.input)
			m = 0
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				if args.prefix == None: name = csv_smiles.loc[i, 'code_name']
				else: name = 'comp_'+str(m)+'_'+csv_smiles.loc[i, 'code_name']

				smi = csv_smiles.loc[i, 'SMILES']
				#checking for salts

				if args.metal_complex == True and args.complex_coord == 6:
					#find metal and replace with I+ for octahydral
					smi = smi.replace(args.metal,'I+')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I++','I+')

					#taking the largest piece
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

				elif args.metal_complex == True and args.complex_coord == 4:
					#find metal and replace with I+ for octahydral
					smi = toks[0].replace(args.metal,'I-')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I-+','I-')

					#taking largest component for salts
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

				elif args.metal_complex == True and args.complex_coord == 5:
					#find metal and replace with I+ for octahydral
					smi = toks[0].replace(args.metal,'I')
					#Ir+ exixted then we need to change back to I+
					smi = smi.replace('I+','I')

					#taking largest component for salts
					pieces = smi.split('.')
					if len(pieces) > 1:
						smi = max(pieces, key=len) #take largest component by length

				# if args.metal_complex == True and  args.complex_coord == 7:
				# 	#find metal and replace with I+ for octahydral
				# 	smi = toks[0].replace(args.metal,'I+2')
				# 	#Ir+ exixted then we need to change back to I+
				# 	smi = smi.replace('I+2+','I+2')
				#
				# if args.metal_complex == True and  args.complex_coord == 8:
				# 	#find metal and replace with I+ for octahydral
				# 	smi = toks[0].replace(args.metal,'I+3')
				# 	#Ir+ exixted then we need to change back to I+
				# 	smi = smi.replace('I+3+','I+3')

				if args.verbose: print("   -> Input Molecule {} is {}".format(i, smi))
				mol = Chem.MolFromSmiles(smi)
				if complex_type == 'squareplanar' or complex_type == 'squarepyrimidal':
					file_template = 'template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])
				m += 1

		elif os.path.splitext(args.input)[1] == '.cdx': # CDX file
				#converting to smiles from chemdraw
			subprocess.run(['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi'])
			smifile = open('cdx.smi',"r")

			for i, line in enumerate(smifile):
				name = 'comp' + str(i)+'_'
				if args.metal_complex == True and args.complex_coord == 6 :
					#find metal and replace with I+ for octahydral
					line = line.replace(args.metal,'I+')
					#Ir+ exixted then we need to change back to I+
					line = line.replace('I++','I+')

				elif args.metal_complex == True and args.complex_coord == 4:
					#find metal and replace with I+ for octahydral
					line = line.replace(args.metal,'I-')
					#Ir+ exixted then we need to change back to I+
					line = line.replace('I-+','I+')

				elif args.metal_complex == True and args.complex_coord == 5:
					#find metal and replace with I+ for octahydral
					line = line.replace(args.metal,'I')
					#Ir+ exixted then we need to change back to I+
					line = line.replace('I+','I')

				mol = Chem.MolFromSmiles(line)
				if complex_type == 'squareplanar' or complex_type == 'squarepyrimidal':
					file_template = 'template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])

		elif os.path.splitext(args.input)[1] == '.com' or os.path.splitext(args.input)[1] == '.gjf': # COM file
				#converting to sdf from comfile to preserve geometry

			comfile = open(args.input,"r")
			comlines = comfile.readlines()

			emptylines=[]

			for i in range(0,len(comlines)):
				if len(comlines[i].strip()) == 0: emptylines.append(i)

			#assigning the charges
			args.charge = comlines[(emptylines[1]+1)].split(' ')[0]

			xyzfile = open(os.path.splitext(args.input)[0]+'.xyz',"w")
			xyzfile.write(str(emptylines[2]- (emptylines[1]+2)))
			xyzfile.write('\n')
			xyzfile.write(os.path.splitext(args.input)[0])
			xyzfile.write('\n')
			for i in range((emptylines[1]+2), emptylines[2]):
				xyzfile.write(comlines[i])

			xyzfile.close()
			comfile.close()

			subprocess.run(['obabel', '-ixyz', os.path.splitext(args.input)[0]+'.xyz', '-osdf', '-O', os.path.splitext(args.input)[0]+'.sdf'])

			# obConversion = ob.OBConversion()
			# obConversion.SetInAndOutFormats("com", "sdf")
			#
			# mol = ob.OBMol()
			# obConversion.ReadFile(mol,os.path.splitext(args.input)[0]+'.com')
			# obConversion.WriteFile(mol, os.path.splitext(args.input)[0]+'.sdf')

			sdffile = os.path.splitext(args.input)[0]+'.sdf'

			IDs = []
			f = open(sdffile,"r")
			readlines = f.readlines()

			for i, line in enumerate(readlines):
				if line.find('>  <ID>') > -1:
					ID = readlines[i+1].split()[0]
					IDs.append(ID)
				else: IDs.append(os.path.splitext(args.input)[0])

			suppl = Chem.SDMolSupplier(sdffile)

			for mol in suppl:
				#FInding the metal and replacing it for RDkit embedding
				if args.metal_complex == True:
					for el in elementspt:
						if el.symbol == 'I':
							atomic_number = el.number
					#changing name for Metal
					for atom in mol.GetAtoms():
						if atom.GetSymbol() in args.metal:
							atom.SetAtomicNum(atomic_number)
							if len(atom.GetNeighbors()) == 2:
								atom.SetFormalCharge(-3)
							if len(atom.GetNeighbors()) == 3:
								atom.SetFormalCharge(-2)
							if len(atom.GetNeighbors()) == 4:
								atom.SetFormalCharge(-1)
							if len(atom.GetNeighbors()) == 5:
								atom.SetFormalCharge(0)
							if len(atom.GetNeighbors()) == 6:
								atom.SetFormalCharge(1)
							if len(atom.GetNeighbors()) == 7:
								atom.SetFormalCharge(2)
							if len(atom.GetNeighbors()) == 8:
								atom.SetFormalCharge(3)

				if args.prefix == None: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				if complex_type == 'squareplanar' or complex_type == 'squarepyrimidal':
					file_template = 'template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])

#------------------ Check for metals ----------------------------------------------

		elif file_format == '.sdf': # SDF input specified
			sdffile = args.input
			IDs = []
			f = open(sdffile,"r")
			readlines = f.readlines()

			for i, line in enumerate(readlines):
				if line.find('>  <ID>') > -1:
					ID = readlines[i+1].split()[0]
					IDs.append(ID)
				else: IDs.append(str(i))

			suppl = Chem.SDMolSupplier(sdffile)

			for i, mol in enumerate(suppl):
				#FInding the metal and replacing it for RDkit embedding
				if args.metal_complex == True:
					#changing name for Metal
					for atom in mol.GetAtoms():
						if atom.GetSymbol() == args.metal:
							for el in elementspt:
								if el.symbol == 'I':
									atomic_number = el.number
							atom.SetAtomicNum(atomic_number)
							if len(atom.GetNeighbors()) == 4:
								atom.SetFormalCharge(-1)
							if len(atom.GetNeighbors()) == 6:
								atom.SetFormalCharge(1)

				if args.prefix == None: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				if complex_type == 'squareplanar' or complex_type == 'squarepyrimidal':
					file_template = 'template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])

#------------------------------------------------------------------------------------------

		for [mol, name] in mol_objects: # Run confomer generation for each mol object
			conformer_generation(mol,name,args)

		conf_files = [name+'_confs.sdf' for mol,name in mol_objects]

		# names for directories created
		sp_dir = 'sp'
		g_dir = 'gaussian'

		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:

					if args.single_point ==  True: # only create this directory if single point calculation is requested
						folder = sp_dir + '/' + str(lot) + '-' + str(bs)
						print("\no  PREPARING SINGLE POINT INPUTS in {}".format(folder))
						try: os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder): pass
							else: raise

						#writing the com files
						for file in conf_files: # check conf_file exists, parse energies and then write dft input
							if os.path.exists(file):
								if args.verbose: print("   -> Converting from {}".format(file))
								energies = read_energies(file)
								name = os.path.splitext(file)[0]
								write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args)

					else: # else create the directory for optimizations
						folder = g_dir + '/' + str(lot) + '-' + str(bs)
						print("\no  PREPARING GAUSSIAN INPUTS in {}".format(folder))
						try: os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder): pass
							else: raise

						#writing the com files
						for file in conf_files: # check conf_file exists, parse energies and then write dft input
							if os.path.exists(file):
								if args.verbose: print("   -> Converting from {}".format(file))
								energies = read_energies(file)
								name = os.path.splitext(file)[0]
								write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args)

	if args.exp_rules == True:

		if args.verbose == True: print("   ----- Applying experimental rules to write the new confs file -----")

		conf_files = glob.glob('*.sdf')

		for file in conf_files:
			print('vfodbobnoigbogibgoi')
			allmols = Chem.SDMolSupplier(file, removeHs=False)
			if allmols is None:
				print("Could not open ", file)
				sys.exit(-1)

			sdwriter = Chem.SDWriter(file.split('.')[0]+exp_rules_output_ext)

			for mol in allmols:
				check_mol = True
				check_mol = exp_rules_output(mol,args)
				if check_mol == True:
					sdwriter.write(mol)
			sdwriter.close()

	if args.analysis == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/'
					#check if New_Gaussian_Input_Files folder exists
					w_dir = check_for_final_folder(w_dir)
					#assign the path to the finished directory.
					w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/Finished'
					#print(w_dir)
					os.chdir(w_dir)
					print(w_dir)
					log_files = glob.glob('*.log')
					output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin)

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.resubmit == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'
				#check if New_Gaussian_Input_Files folder exists
				w_dir = check_for_final_folder(w_dir)
				os.chdir(w_dir)
				cmd = args.submission_command + ' *.com'
				if args.qsub == True:
					os.system(cmd)

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.nmr == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					w_dir = args.path + str(lot) + '-' + str(bs) +'/Finished'
					if os.path.isdir(w_dir):
						os.chdir(w_dir)
						log_files = glob.glob('*.log')
						output_analyzer(log_files, w_dir, lot, bs,bs_gcp, args)
					else:
						pass

	#once all files are finished are in the Finished folder
	if args.dup == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'Finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				try:
					log_files = glob.glob('*.log')
					if len(log_files) != 0:
						val = ' '.join(log_files)
						dup_calculation(val,w_dir)
					else:
						print(' Files for are not there!')
					#print(log_files)
				except:
					pass

	#once all files are finished are in the Finished folder
	if args.boltz == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'Finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				for i in range(args.maxnumber):
					#grab all the corresponding files make sure to renamme prefix when working with differnet files
					try:
						log_files = glob.glob('RE' + '_' + str(i)+'_'+'confs_low.log')
						if len(log_files) != 0:
							val = ' '.join(log_files)
							boltz_calculation(val,i)
						else:
							print(' Files for {} are not there!'.format(i))
						#print(log_files)
					except:
						pass

	if args.combine == True:
		#combines the files and gives the boltzmann weighted energies
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) + '/Finished'
				os.chdir(w_dir)
				#read the csv log_files
				csv_files = glob.glob('Goodvibes*.csv')
				combine_files(csv_files, lot, bs, args)

if __name__ == "__main__":
	main()
