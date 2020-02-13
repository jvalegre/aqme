'''
#ADD description
DBcg v1.0.0
Authors: Shree Sowndarya S. V., Juan V. Alegre Requena,
please report any bugs to svss@colostate.edu or juanvi89@hotmail.com.
'''

from db_gen_functions import *

possible_atoms = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
	"Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
columns = ['Structure', 'E', 'ZPE', 'H', 'T.S', 'T.qh-S', 'G(T)', 'qh-G(T)']

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")

	parser.add_argument("--var",dest="var",action="store_true",default=False,help="If providing a Variable files")
	parser.add_argument("--varfile",dest="varfile",help="Variable file")
	parser.add_argument("--input",help="Input smi file pr SDF files")
	parser.add_argument("-m","--maxnumber", help="Number of compounds", type=int, metavar="maxnumber")

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
	parser.add_argument("--basis_set",  help="Basis Set", default=['6-31g**'], dest="basis_set", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms",default=['LANL2DZ'], help="Basis Set genecp: The length has to be the same as basis_set", dest="basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--genecp_atoms",  help="genecp atoms",default=[], dest="genecp_atoms",type=str, nargs='*')
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

	### If the input file is python format then we will assume it contains variables ... ###
	if args.var == True:
		if os.path.splitext(args.varfile)[1] == '.py':
		   print("o  READING VARIABLES FROM", args.varfile)
		   print(os.path.splitext(args.varfile)[0])
		   db_gen_variables = os.path.splitext(args.varfile)[0]
		   import db_gen_variables as args

	if args.time: start_time = time.time()

### check if it's working with *.smi (if it isn't, we need to change this a little bit)
### no it isnt working for *.smi, same or .sdf. We can add that feature in the end.
	if args.compute == True:
		#checking if .smi
		if os.path.splitext(args.input)[1] == '.smi':
			smifile = open(args.input)
			root = os.path.splitext(args.input)[0]

			m = 0
			for line in smifile:
				toks = line.split()
				smi = toks[0]
				name = 'comp_'+str(m)+'_'+''.join(toks[1:])

				"""
				if len(name) == 0: name = root
				pieces = smi.split('.')
				if len(pieces) > 1:
					smi = max(pieces, key=len) #take largest component by length
					print("Taking largest component: %s\t%s" % (smi,name))
					"""
					# finally converts each line to a rdkit mol object
				mol = Chem.MolFromSmiles(smi)
				print(Chem.MolToSmiles(mol))
					#doing confomer genertion for each mol object
				conformer_generation(mol,name,args)
				m += 1
					#checking if .sdf
		elif os.path.splitext(args.input)[1] == '.sdf':
			IDs = []
			sdffile = args.input
			root = os.path.splitext(args.input)[0]

			f = open(sdffile,"r")
			readlines = f.readlines()
			for i in range(len(readlines)):
				if readlines[i].find('>  <ID>') > -1:
					ID = readlines[i+1].split()[0]
					#if readlines[i].find('>  <NAME>') > -1:
					#name = readlines[i+1].split()[0]
					IDs.append(ID)
			suppl = Chem.SDMolSupplier(sdffile)
			i=0
			m=0
			for mol in suppl:
				name = 'comp_'+str(m)+'_'+IDs[i]
				#doing confomer genertion for each mol object
				conformer_generation(mol,name,args)
				i += 1
				m += 1

		#provide a .csv file with one columns SMILES and code_name
		elif os.path.splitext(args.input)[1] == '.csv':
			csv_smiles = pd.read_csv(args.input)
			m = 0
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				name = 'comp_'+str(m)+'_'+csv_smiles.loc[i, 'Ir_cat']
				smi = csv_smiles.loc[i, 'SMILES']
				mol = Chem.MolFromSmiles(smi)
				conformer_generation(mol,name,args)
				m += 1

		#checking if .cdx
		elif os.path.splitext(args.input)[1] == '.cdx':
			#converting to smiles from chemdraw
			subprocess.run(['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi'])
			smifile = open('cdx.smi',"r")

			m=0
			for line in smifile:
				name = 'comp' + str(m)+'_'
				mol = Chem.MolFromSmiles(line)
				conformer_generation(mol,name,args)
				m += 1

			# creating the com files from the sdf files created read *_confs.sdf
			# READING THE SDF

		sdf_to_gjf_files = glob.glob('*_confs.sdf')


		for lot in args.level_of_theory:
			for bs in args.basis_set:
				for bs_gcp in args.basis_set_genecp_atoms:
					if args.single_point ==  True:
						folder = 'sp/' + str(lot) + '-' + str(bs)

						try:
							os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder):
								pass
							else:
								raise

							#writing the com files
						for file in sdf_to_gjf_files:
							energies = []
							f = open(file,"r")
							readlines = f.readlines()
							for i in range(len(readlines)):
								if readlines[i].find('>  <Energy>') > -1:
									energies.append(float(readlines[i+1].split()[0]))
							f.close()
							name = os.path.splitext(file)[0]
							write_gaussian_input_file(file,name,lot, bs, bs_gcp, energies, args)
					else:

						folder = 'gaussian/' + str(lot) + '-' + str(bs)


						try:
							os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder):
								pass
							else:
								raise
							#writing the com files
						for file in sdf_to_gjf_files:
							energies = []
							f = open(file,"r")
							readlines = f.readlines()
							for i in range(len(readlines)):
								if readlines[i].find('>  <Energy>') > -1:
									energies.append(float(readlines[i+1].split()[0]))
							f.close()
							print(energies)
							name = os.path.splitext(file)[0]
							write_gaussian_input_file(file,name,lot, bs, bs_gcp, energies, args)

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
					#print(log_files)
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
