'''
### ADD DESCRIPTION
Authors: Shree Sowndarya S. V., Juan V. Alegre Requena,
please report any bugs to svss@colostate.edu or juanvi89@hotmail.com.
'''

from db_gen_functions import *
from db_gen_variables import *

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")
	parser.add_argument("--input",help="Input smi file pr SDF files")
	parser.add_argument("-c","--compute", action="store_true", default=False, help="Create input files for Gaussian")
	parser.add_argument("-a","--analysis", action="store_true", default=False, help="Fix and analyze Gaussian output files")
	parser.add_argument("-r","--resubmit", action="store_true", default=False, help="Resubmit Gaussian input files")
	parser.add_argument("-s","--secondrun", action="store_true", default=False, help="Set true for second run analysis")
	parser.add_argument("-n","--nmr",action="store_true",default=False, help="Create Files for Single Point which includes NMR calculation after DFT Optimization")
	parser.add_argument("-b","--boltz", action="store_true", default=False, help="Boltzmann factor for each conformers from Gaussian output files")
	parser.add_argument("-d","--combine", action="store_true", default=False, help="Combine files of differnt molecules including boltzmann weighted energies")
	#pass the argument for path for the gaussian folder.
	parser.add_argument("--path", help="Path for analysis/boltzmann factor/combining files where the gaussian folder created is present")
	parser.add_argument("-v","--verbose",action="store_true",default=False, help="verbose output")
	args = parser.parse_args()

### check if it's working with *.smi (if it isn't, we need to change this a little bit)
### no it isnt working for *.smi, same or .sdf. We can add that feature in the end.
	if args.compute == True:
		#checking if .smi
		if os.path.splitext(args.input)[1] == '.smi':
			smifile = open(args.input)
			root = os.path.splitext(args.input)[0]

			for line in smifile:
				toks = line.split()
				smi = toks[0]
				name = ''.join(toks[1:])

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
			for mol in suppl:
				name = IDs[i]
				#doing confomer genertion for each mol object
				conformer_generation(mol,name,args)
				i += 1

		#provide a .csv file with one columns SMILES and code_name
		elif os.path.splitext(args.input)[1] == '.csv':
			csv_smiles = pd.read_csv(args.input)
			for i in range(len(csv_smiles)):
				#assigning names and smi i  each loop
				name = csv_smiles.loc[i, 'code_name']
				smi = csv_smiles.loc[i, 'SMILES']
				mol = Chem.MolFromSmiles(smi)
				conformer_generation(mol,name,args)

		#checking if .cdx
		elif os.path.splitext(args.input)[1] == '.cdx':
			#converting to smiles from chemdraw
			subprocess.run(['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi'])
			smifile = open('cdx.smi',"r")

			i=1
			for line in smifile:
				name = 'mol_' + str(i)
				mol = Chem.MolFromSmiles(line)
				conformer_generation(mol,name,args)
				i += 1

			# creating the com files from the sdf files created read *_confs.sdf
			# READING THE SDF

		sdf_to_gjf_files = glob.glob('*_confs.sdf')

		for lot in level_of_theory:
			for bs in basis_set:
				if single_point ==  True:
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
						name = os.path.splitext(file)[0]
						write_gaussian_input_file(file,name,lot, bs)
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
						name = os.path.splitext(file)[0]
						write_gaussian_input_file(file,name,lot, bs)

	if args.analysis == True:
		# Sets the folder and find the log files to analyze
		for lot in level_of_theory:
			for bs in basis_set:
				if args.secondrun != True:
					w_dir = args.path + str(lot) + '-' + str(bs)+'/'
					os.chdir(w_dir)
					log_files = glob.glob('*.log')
					output_analyzer(log_files, w_dir, lot, bs, nprocs, mem, args)

				else:
					w_dir = args.path + str(lot) + '-' + str(bs)+'/New_Gaussian_Input_Files/'
					if os.path.isdir(w_dir):
						os.chdir(w_dir)
						log_files = glob.glob('*.log')
						output_analyzer(log_files, w_dir, lot, bs, nprocs, mem, args)
					else:
						pass

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.resubmit == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in level_of_theory:
			for bs in basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs)+'/New_Gaussian_Input_Files'
				if  os.path.isdir(w_dir):
					os.chdir(w_dir)
					cmd = 'qsub_summit' + ' *.com'
					os.system(cmd)
				else:
					pass

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.nmr == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in level_of_theory:
			for bs in basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs)+'/Finished'
				if os.path.isdir(w_dir):
					os.chdir(w_dir)
					log_files = glob.glob('*.log')
					output_analyzer(log_files, w_dir, lot, bs, nprocs, mem, args)
				else:
					pass

	#once all files are finished are in the Finished folder
	if args.boltz == True:
		# Sets the folder and find the log files to analyze
		for lot in level_of_theory:
			for bs in basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs)+'/'+'Finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				for i in molecules:
					#grab all the corresponding files make sure to renamme
					log_files = glob.glob('RE_'+ str(i)+ '_confs_' + '*.log')
					#print(log_files)
					val = ' '.join(log_files)
					boltz_calculation(val,i)

	if args.combine == True:
		#combines the files and gives the boltzmann weighted energies
		for lot in level_of_theory:
			for bs in basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) + '/Finished'
				os.chdir(w_dir)
				#read the csv log_files
				csv_files = glob.glob('Goodvibes*.csv')
				combine_files(csv_files, lot, bs, args)
