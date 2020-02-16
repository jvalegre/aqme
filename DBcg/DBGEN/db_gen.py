'''
#ADD description
DBcg v1.0.0
Authors: Shree Sowndarya S. V., Juan V. Alegre Requena,
please report any bugs to svss@colostate.edu or juanvi89@hotmail.com.
'''

from db_gen_functions import *

possible_atoms = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
                 "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                 "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                 "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
                 "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
                 "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
                 "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                 "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"]
columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']

def main():
	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")

	parser.add_argument("--varfile",dest="varfile",default=None,help="Parameters in python format")
	parser.add_argument("--input",help="Molecular structure")
	parser.add_argument("-m","--maxnumber", help="Number of compounds", type=int, metavar="maxnumber")
	parser.add_argument("--prefix", help="Prefix for naming files", default=None, metavar="prefix")

	parser.add_argument("-w","--compute", action="store_true", default=False, help="Create input files for Gaussian")
	parser.add_argument("-a","--analysis", action="store_true", default=False, help="Fix and analyze Gaussian output files")
	parser.add_argument("-r","--resubmit", action="store_true", default=False, help="Resubmit Gaussian input files")
	parser.add_argument("-s","--secondrun", action="store_true", default=False, help="Set true for second run analysis")
	parser.add_argument("-n","--nmr",action="store_true",default=False, help="Create Files for Single Point which includes NMR calculation after DFT Optimization")
	parser.add_argument("-b","--boltz", action="store_true", default=False, help="Boltzmann factor for each conformers from Gaussian output files")
	parser.add_argument("-f","--combine", action="store_true", default=False, help="Combine files of differnt molecules including boltzmann weighted energies")
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

	#arguments for gaussian files Creation
	parser.add_argument("-l", "--level_of_theory",help="Level of Theory", default=['wB97xd'], dest="level_of_theory", type=str, nargs='*')
	parser.add_argument("--basis_set",  help="Basis Set", default=['6-31g*'], dest="basis_set", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms",default=['LANL2DZ'], help="Basis Set genecp: The length has to be the same as basis_set", dest="basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--genecp_atoms",  help="genecp atoms",default=[], dest="genecp_atoms",type=str, nargs='*')
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


	if args.varfile != None: args.var = True
	if args.time: start_time = time.time()

	### If the input file is python format then we will assume it contains variables ... ###
	if args.var == True:
		if os.path.splitext(args.varfile)[1] == '.py':
			db_gen_variables = os.path.splitext(args.varfile)[0]
			print("\no  IMPORTING VARIABLES FROM", args.varfile)
			args = __import__(db_gen_variables)

### check if it's working with *.smi (if it isn't, we need to change this a little bit)
### no it isnt working for *.smi, same or .sdf. We can add that feature in the end.

	if args.compute == True: # this will perform conformational analysis and create inputs for Gaussian

		# input file format specified
		[file_name, file_format] = os.path.splitext(args.input)

		if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com']:
			print("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!"); sys.exit()
		else:
			mol_objects = [] # a list of mol objects that will be populated

		if file_format == '.smi': # SMILES input specified
			smifile = open(args.input)

			for i, line in enumerate(smifile):
				toks = line.split()
				smi = toks[0]
				if args.prefix == None: name = ''.join(toks[1:])
				else: name = args.prefix+str(i)+'_'+''.join(toks[1:])

				# Converts each line to a rdkit mol object
				if args.verbose: print("   -> Input Molecule {} from {}".format(smi, args.input))
				mol = Chem.MolFromSmiles(smi)
				mol_objects.append([mol, name])

		elif file_format == '.sdf': # SDF input specified
			sdffile = open(args.input)
			IDs = []
			f = open(sdffile,"r")
			readlines = f.readlines()

			for i, line in enumerate(readlines):
				if line.find('>  <ID>') > -1:
					ID = readlines[i+1].split()[0]
					IDs.append(ID)
				else: IDs.append(i)

			suppl = Chem.SDMolSupplier(sdffile)

			for i, mol in enumerate(suppl):
				if args.prefix == None: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				mol_objects.append([mol, name])

		elif os.path.splitext(args.input)[1] == '.csv': # CSV file with one columns SMILES and code_name
			csv_smiles = pd.read_csv(args.input)
			m = 0
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				if args.prefix == None: name = csv_smiles.loc[i, 'code_name']
				else: name = 'comp_'+str(m)+'_'+csv_smiles.loc[i, 'code_name']
				smi = csv_smiles.loc[i, 'SMILES']
				mol = Chem.MolFromSmiles(smi)
				mol_objects.append([mol, name])
				m += 1

		elif os.path.splitext(args.input)[1] == '.cdx': # CDX file
			#converting to smiles from chemdraw
			subprocess.run(['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi'])
			smifile = open('cdx.smi',"r")

			for i, line in enumerate(smifile):
				name = 'comp' + str(i)+'_'
				mol = Chem.MolFromSmiles(line)
				mol_objects.append([mol, name])

		elif os.path.splitext(args.input)[1] == '.com': # COM file
			#converting to sdf from comfile to preserve geometry

			comfile = open(args.input,"r")
			comlines = comfile.readlines()

			for i in range(0,len(comlines)):
				if comlines[i][0].isdigit() and len(comlines[i-1].strip()) == 0:standor_initial = i+1
				if len(comlines[i].strip()) == 0:standor_final = i

			xyzfile = open(os.path.splitext(args.input)[0]+'.xyz',"w")

			xyzfile.write(str(standor_final - standor_initial))
			xyzfile.write('\n')
			xyzfile.write(os.path.splitext(args.input)[0])
			xyzfile.write('\n')
			for i in range(standor_initial, standor_final):
				xyzfile.write(comlines[i])

			xyzfile.close()
			comfile.close()

			subprocess.run(['obabel', '-ixyz', os.path.splitext(args.input)[0]+'.xyz', '-osdf', '-O', os.path.splitext(args.input)[0]+'.sdf'])

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

			for i, mol in enumerate(suppl):
				if args.prefix == None: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				mol_objects.append([mol, name])

		for [mol, name] in mol_objects: # Run confomer generation for each mol object
			conformer_generation(mol,name,args)

		# creating the com files from the sdf files created read *_confs.sdf
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

	if args.analysis == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					if args.secondrun != True:
						w_dir = args.path + str(lot) + '-' + str(bs) +'/'
						os.chdir(w_dir)
						log_files = glob.glob('*.log')
						output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args)

					else:
						w_dir = args.path + str(lot) + '-' + str(bs) +'/New_Gaussian_Input_Files/'
						if os.path.isdir(w_dir):
							os.chdir(w_dir)
							log_files = glob.glob('*.log')
							output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args)
						else:
							pass

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.resubmit == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/New_Gaussian_Input_Files'
				if  os.path.isdir(w_dir):
					os.chdir(w_dir)
					cmd = submission_command + ' *.com'
					os.system(cmd)
				else:
					pass

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
						log_files = glob.glob(comp + '_' + str(i)+'-'+'*.log')
					except:
						pass
					val = ' '.join(log_files)
					boltz_calculation(val,i)

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
