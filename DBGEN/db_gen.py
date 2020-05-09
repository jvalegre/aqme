"""######################################################################################
#########################################################################################
###																					  ###
###  RotaConfort is a tool that allows to carry out automated:						  ###
###  (1) Conformational searches and creation of COM files using RDKit, xTB and ANI1  ###
###  (2) LOG file processing (detects imaginary freqs and error terminations		  ###
###      and creates new COM files)													  ###
###  (3) Use LOG files to create new COM files with new keywords (i.e. single-point   ###
###      corrections after geometry optimization)									  ###
###  																				  ###
#########################################################################################
###  																				  ###
###  Version: v1.0, Release date: 24-April-2020										  ###
###  																				  ###
#########################################################################################
###  																				  ###
###  Authors: Shree Sowndarya S. V., Juan V. Alegre Requena, Robert S. Paton		  ###
###  																				  ###
###  Please, report any bugs or suggestions to:										  ###
###  svss@colostate.edu or juanvi89@hotmail.com  									  ###
###																					  ###
#########################################################################################
######################################################################################"""


from DBGEN.db_gen_functions import *
from DBGEN.confgen import *
#from DBGEN.Template import *

def main():
	parser = argparse.ArgumentParser(description="Generate conformers depending on type of optimization (change parameters in db_gen_PATHS.py file).")

	#Input details
	parser.add_argument("--varfile", dest="varfile", default=None, help="Parameters in YAML format")
	parser.add_argument("-i", "--input", help="File containing molecular structure(s)")
	parser.add_argument("--output_name", dest="output_name", default="output", metavar="output_name", help="Change output filename to DBGEN_\"output\".dat")
	parser.add_argument("--output", dest="output", default=".sdf", metavar="output", help="The extension of the SDF files written")

	#metal complex
	parser.add_argument("--metal_complex", action="store_true", default=False, help="Request metal complex with coord. no. 4, 5 or 6")
	parser.add_argument("--metal",  help="Specify metallic element", default=["Ir"], dest="metal", type=str)
	parser.add_argument("--complex_spin",  help="Multiplicity of metal complex", default="1", dest="complex_spin", type=int)
	parser.add_argument("--complex_coord", help="Coord. no. of metal complex (automatically updates)", default=[], dest="complex_coord", type=int)
	parser.add_argument("--complex_type",  help="Geometry about metal (e.g. octahedral)", default="octahedral", dest="complex_type", type=str)
	parser.add_argument("--m_oxi",  help="Metal oxidation state", default=["3"], dest="m_oxi", type=int)
	parser.add_argument("--metal_idx",  help="Metal index (automatically updates)", default=[], dest="metal_idx", type=int)
	parser.add_argument("--charge",  help="Charge of metal complex (automatically updates)", default=[], dest="charge", type=int)
	parser.add_argument("--charge_default",  help="Charge default to be considered", default="0", dest="charge_default", type=int)
	parser.add_argument("--metal_sym",  help="Symbols of metals to be considered from list (automatically updates)", default=[], dest="metal_sym", type=str)

	#NCI complex
	parser.add_argument("--nci_complex", action="store_true", default=False, help="Request NCI complexes")
	parser.add_argument("-m", "--maxnumber", help="Number of compounds", type=int, metavar="maxnumber")
	parser.add_argument("--prefix", help="Prefix for naming files", default=None, metavar="prefix")

	#work the script has to do
	parser.add_argument("-w", "--compute", action="store_true", default=False, help="Perform conformational analysis")
	parser.add_argument("--write_gauss", action="store_true", default=False, help="Create input files for Gaussian")
	parser.add_argument("-a", "--analysis", action="store_true", default=False, help="Fix and analyze Gaussian outputs")
	parser.add_argument("-r", "--resubmit", action="store_true", default=False, help="Resubmit Gaussian input files")

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
	parser.add_argument("--opt_fmax", action="store",default=0.05, help="fmax value used in xTB and AN1 optimizations", type=float)
	parser.add_argument("--opt_steps", action="store",default=1000, help="max cycles used in xTB and AN1 optimizations", type=int)
	parser.add_argument("--opt_steps_RDKit", action="store",default=1000, help="max cycles used in RDKit optimizations", type=int)
	parser.add_argument("--time","-t",action='store_true',default=False,help="request program runtime")
	parser.add_argument("--heavyonly", help="only consider torsion angles involving heavy (non H) elements (default=True)", default=True, metavar="heavyonly")
	parser.add_argument("--constraints", action="store_true",default=None, help="distance constraint")
	parser.add_argument("--nodihedrals", action="store_true", default=False, help="turn off dihedral scan")
	parser.add_argument("-d","--degree", type=float,help="Amount, in degrees, to enumerate torsions by (default 30.0)",default=30.0)
	parser.add_argument("--etkdg", dest="etkdg",action="store_true",default=False, help="use new ETKDG knowledge-based method instead of distance geometry")
	parser.add_argument("--max_torsions",type=int,help="Skip any molecules with more than this many torsions (default 5)",default=5)
	parser.add_argument("--sample", help="number of conformers to sample to get non-torsional differences (default 100)", default=100, type=int, metavar="sample")
	parser.add_argument("--auto_sample", help="final factor to multiply in the auto mode for the sample option (default 20)", default=20, type=int, metavar="auto_sample")
	parser.add_argument("--ff", help="force field (MMFF or UFF)", default="MMFF", metavar="ff")
	parser.add_argument("--seed", help="random seed (default 062609)", default="062609", type=int, metavar="s")
	parser.add_argument("--rms_threshold", help="cutoff for considering sampled conformers the same (default 0.25)", default=0.25, type=float, metavar="R")
	parser.add_argument("--max_matches_RMSD", help="iteration cutoff for considering  matches in sampled conformers the same (default 1000000 )", default=1000000 , type=float, metavar="max_matches_RMSD")
	parser.add_argument("--energy_threshold", dest="energy_threshold",action="store",default=0.05, help="energy difference between unique conformers")
	parser.add_argument("--initial_energy_threshold", dest="initial_energy_threshold",action="store",default=0.01, help="energy difference between unique conformers for the first filter of only E")
	parser.add_argument("--max_MolWt", help="Max. molecular weight of molecule", default=1000, type=int, metavar="max_MolWt")
	parser.add_argument("--num_rot_bonds", help="Max. number of rotatable bonds in a molecule", default=20, type=int, metavar="num_rot_bonds")
	parser.add_argument("--large_sys", action="store_true",default=False, help="Large systems for xtb optimizations")
	parser.add_argument("--STACKSIZE", help="STACKSIZE for optimization of large systems", default="500m")

	#arguments for gaussian files Creation
	parser.add_argument("-l", "--level_of_theory",help="Level of Theory", default=['wB97xd'], dest="level_of_theory", type=str, nargs='*')
	parser.add_argument("--basis_set",  help="Basis Set", default=['6-31g*'], dest="basis_set", type=str, nargs='*')
	parser.add_argument("--basis_set_genecp_atoms",default=['LANL2DZ'], help="Basis Set genecp/gen: The length has to be the same as basis_set", dest="basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--genecp_atoms",  help="genecp atoms",default=[], dest="genecp_atoms",type=str, nargs='*')
	parser.add_argument("--gen_atoms",  help="gen atoms",default=[], dest="gen_atoms",type=str, nargs='*')
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

	#define the logging object
	log = Logger("DBGEN", args.output_name)

	start_time = time.time()

	### Variables will be updated from YAML file ###
	if args.varfile != None:
		if os.path.exists(args.varfile):
			if os.path.splitext(args.varfile)[1] == '.yaml':
				log.write("\no  IMPORTING VARIABLES FROM " + args.varfile)
				with open(args.varfile, 'r') as file:
					param_list = yaml.load(file, Loader=yaml.FullLoader)
	for param in param_list:
		if hasattr(args, param):
			if getattr(args, param) != param_list[param]:
				log.write("o  RESET " + param + " from " + str(getattr(args, param)) + " to " + str(param_list[param]))
				setattr(args, param, param_list[param])
			else:
				log.write("o  DEFAULT " + param + " : " + str(getattr(args, param)))

### check if it's working with *.smi (if it isn't, we need to change this a little bit)
### no it isnt working for *.smi, same or .sdf. We can add that feature in the end.

	if args.compute == True: # this will perform conformational analysis and create inputs for Gaussian

		# input file format specified
		[file_name, file_format] = os.path.splitext(args.input)

		if file_format not in ['.smi', '.sdf', '.cdx', '.csv','.com','.gjf']:
			log.write("\nx  INPUT FILETYPE NOT CURRENTLY SUPPORTED!"); sys.exit()
		else:
			mol_objects = [] # a list of mol objects that will be populated

		# writing the list of DUPLICATES
		if args.nodihedrals == True:
			if args.xtb != True and args.ANI1ccx != True:
				dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','time (seconds)','Overall charge'])
			elif args.xtb == True:
				dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','xTB-Initial-samples','xTB-initial_energy_threshold','xTB-RMS-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
			elif args.ANI1ccx == True:
				dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-initial_energy_threshold','ANI1ccx-RMS-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])
		else:
			if args.xtb != True and args.ANI1ccx != True:
				dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','time (seconds)','Overall charge'])
			elif args.xtb == True:
				dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','xTB-Initial-samples','xTB-initial_energy_threshold','xTB-RMS-and-energy-duplicates','xTB-Unique-conformers','time (seconds)','Overall charge'])
			elif args.ANI1ccx == True:
					dup_data =  pd.DataFrame(columns = ['Molecule','RDKIT-Initial-samples', 'RDKit-energy-duplicates','RDKit-RMS-and-energy-duplicates','RDKIT-Unique-conformers','RDKIT-Rotated-conformers','RDKIT-Rotated-Unique-conformers','ANI1ccx-Initial-samples','ANI1ccx-initial_energy_threshold','ANI1ccx-RMS-and-energy-duplicates','ANI1ccx-Unique-conformers','time (seconds)','Overall charge'])

		ori_ff = args.ff

		if file_format == '.smi': # SMILES input specified
			smifile = open(args.input)

			for i, line in enumerate(smifile):
				args.ff = ori_ff
				args.metal_idx = []
				args.complex_coord = []
				args.metal_sym = []

				toks = line.split()

				#editing part
				smi = toks[0]

				#taking largest component for salts
				pieces = smi.split('.')
				if len(pieces) > 1:
					smi = max(pieces, key=len) #take largest component by length

				if args.prefix == False: name = ''.join(toks[1:])
				else: name = args.prefix+str(i)+'_'+''.join(toks[1:])

				# Converts each line to a rdkit mol object
				if args.verbose: log.write("   -> Input Molecule {} is {}".format(i, smi))

				if args.metal_complex == True:
					mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(smi,args,log)
				else:
					mol = Chem.MolFromSmiles(smi)

				# get manually for square planar and SQUAREPYRIMIDAL
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
					if len(args.metal_idx) == 1:
						file_template = path.dirname(path.abspath(__file__)) +'/Template/template-4-and-5.sdf'
						temp = Chem.SDMolSupplier(file_template)
						mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed_sp(mol,temp,name,args,log)
						for i in range(len(mol_objects_from_template)):
							mol_objects.append([mol_objects_from_template[i],name[i],coord_Map[i],alg_Map[i],mol_template[i]])
						for [mol, name, coord_Map,alg_Map,mol_template] in mol_objects:
							conformer_generation(mol,name,start_time,args,log,dup_data,i,coord_Map,alg_Map,mol_template)
					else:
						log.write("x  Cannot use templates for complexes involving more than 1 metal.")
						sys.exit()
				else:
					conformer_generation(mol,name,start_time,args,log,dup_data,i)
			dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

		elif os.path.splitext(args.input)[1] == '.csv': # CSV file with one columns SMILES and code_name
			csv_smiles = pd.read_csv(args.input)
			m = 0
			for i in range(len(csv_smiles)):
				args.ff = ori_ff
				args.metal_idx = []
				args.complex_coord = []
				args.metal_sym = []
				#assigning names and smi i  each loop
				if args.prefix == False: name = csv_smiles.loc[i, 'code_name']
				else: name = 'comp_'+str(m)+'_'+csv_smiles.loc[i, 'code_name']

				smi = csv_smiles.loc[i, 'SMILES']
				#checking for salts
				#taking largest component for salts
				pieces = smi.split('.')
				if len(pieces) > 1:
					smi = max(pieces, key=len) #take largest component by length

				# Converts each line to a rdkit mol object
				if args.verbose: log.write("   -> Input Molecule {} is {}".format(i, smi))

				if args.metal_complex == True:
					mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(smi,args,log)
				else:
					mol = Chem.MolFromSmiles(smi)

				# get manually for square planar and SQUAREPYRIMIDAL
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
					if len(args.metal_idx) == 1:
						file_template = path.dirname(path.abspath(__file__)) +'/Template/template-4-and-5.sdf'
						temp = Chem.SDMolSupplier(file_template)
						mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed_sp(mol,temp,name,args,log)
						for i in range(len(mol_objects_from_template)):
							mol_objects.append([mol_objects_from_template[i],name[i],coord_Map[i],alg_Map[i],mol_template[i]])
						for [mol, name, coord_Map,alg_Map,mol_template] in mol_objects:
							conformer_generation(mol,name,start_time,args,log,dup_data,i,coord_Map,alg_Map,mol_template)
					else:
						log.write("x  Cannont use templates for complexes invovling more than 1 metal.")
						sys.exit()
				else:
					conformer_generation(mol,name,start_time,args,log,dup_data,i)
				m += 1
			dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

		elif os.path.splitext(args.input)[1] == '.cdx': # CDX file
				#converting to smiles from chemdraw
			subprocess.run(['obabel', '-icdx', args.input, '-osmi', '-O', 'cdx.smi'])
			smifile = open('cdx.smi',"r")

			for i, smi in enumerate(smifile):
				args.ff = ori_ff
				name = 'comp' + str(i)+'_'

				#taking largest component for salts
				pieces = smi.split('.')
				if len(pieces) > 1:
					smi = max(pieces, key=len) #take largest component by length

				# Converts each line to a rdkit mol object
				if args.verbose: log.write("   -> Input Molecule {} is {}".format(i, smi))

				if args.metal_complex == True:
					mol,args.metal_idx,args.complex_coord,args.metal_sym = substituted_mol(smi,args,log)
				else:
					mol = Chem.MolFromSmiles(smi)

				# get manually for square planar and SQUAREPYRIMIDAL
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
					if len(args.metal_idx) == 1:
						file_template = path.dirname(path.abspath(__file__)) +'/Template/template-4-and-5.sdf'
						temp = Chem.SDMolSupplier(file_template)
						mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed_sp(mol,temp,name,args,log)
						for i in range(len(mol_objects_from_template)):
							mol_objects.append([mol_objects_from_template[i],name[i],coord_Map[i],alg_Map[i],mol_template[i]])
					else:
						log.write("x  Cannont use templates for complexes invovling more than 1 metal.")
						sys.exit()
				else:
					mol_objects.append([mol, name])

		elif os.path.splitext(args.input)[1] == '.com' or os.path.splitext(args.input)[1] == '.gjf': # COM file
			#converting to sdf from comfile to preserve geometry

			#### ADD args.ff = ori_ff
			comfile = open(args.input,"r")
			comlines = comfile.readlines()

			emptylines=[]

			for i in range(0,len(comlines)):
				if len(comlines[i].strip()) == 0: emptylines.append(i)

			#assigning the charges
			args.charge_default = comlines[(emptylines[1]+1)].split(' ')[0]

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

				if args.prefix == False: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				# get manually for square planar and SQUAREPYRIMIDAL
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
					file_template = path.dirname(path.abspath(__file__)) +'/Template/template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template,name, coord_Map, alg_Map, mol_template = template_embed_sp(mol,temp,name,args,log)
					for i in range(len(mol_objects_from_template)):
						mol_objects.append([mol_objects_from_template[i],name[i],coord_Map[i],alg_Map[i],mol_template[i]])
				else:
					mol_objects.append([mol, name])

#------------------ Check for metals ----------------------------------------------

		elif file_format == '.sdf': # SDF input specified


			#### ADD args.ff = ori_ff
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

				if args.prefix == False: name = IDs[i]
				else: name = args.prefix+str(m)+'_'+IDs[i]
				if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
					file_template = path.dirname(path.abspath(__file__)) +'/Template/template-4-and-5.sdf'
					temp = Chem.SDMolSupplier(file_template)
					mol_objects_from_template = template_embed_sp(mol,temp,name,args,log)
					mol_objects.append(mol_objects_from_template)
				else:
					mol_objects.append([mol, name])

#------------------------------------------------------------------------------------------
		# if args.complex_type == 'squareplanar' or args.complex_type == 'squarepyrimidal':
		# 	dup_data_idx = 0
		# 	for [mol, name, coord_Map,alg_Map,mol_template] in mol_objects: # Run confomer generation for each mol object
		# 		conformer_generation(mol,name,start_time,args,log,dup_data,dup_data_idx,coord_Map,alg_Map,mol_template)
		# 		dup_data_idx += 1
		# 	dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)
		# else:
		# 	dup_data_idx = 0
		# 	for [mol, name] in mol_objects: # Run confomer generation for each mol object
		# 		conformer_generation(mol,name,start_time,args,log,dup_data,dup_data_idx)
		# 		dup_data_idx += 1
		# 	dup_data.to_csv(args.input.split('.')[0]+'-Duplicates Data.csv',index=False)

	#applying rule to get the necessary conformers only
	if args.exp_rules == True:
		if args.verbose == True: log.write("   ----- Applying experimental rules to write the new confs file -----")
		### do 2 cases, for RDKit only and RDKIt+xTB
		#grab all the gaussian files
		if args.xtb != True and args.ANI1ccx != True:
			conf_files =  glob.glob('*_rdkit.sdf')
		elif args.xtb == True:
			conf_files =  glob.glob('*_xtb.sdf')
		elif args.ANI1ccx == True:
			conf_files =  glob.glob('*_ani.sdf')

		for file in conf_files:
			allmols = Chem.SDMolSupplier(file, removeHs=False)
			if allmols is None:
				log.write("Could not open "+ file)
				sys.exit(-1)

			sdwriter = Chem.SDWriter(file.split('.')[0]+args.exp_rules_output_ext)

			for mol in allmols:
				check_mol = True
				check_mol = exp_rules_output(mol,args,log)
				if check_mol == True:
					sdwriter.write(mol)
			sdwriter.close()

	if args.write_gauss == True:
		if args.exp_rules == True:
			conf_files =  glob.glob('*_rules.sdf')
		#grad all the gaussian files
		elif args.xtb != True and args.ANI1ccx != True and args.nodihedrals == True:
			conf_files =  glob.glob('*_rdkit.sdf')
		elif args.xtb == True:
			conf_files =  glob.glob('*_xtb.sdf')
		elif args.ANI1ccx == True:
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

					if args.single_point ==  True: # only create this directory if single point calculation is requested
						folder = sp_dir + '/' + str(lot) + '-' + str(bs)
						log.write("\no  PREPARING SINGLE POINT INPUTS in {}".format(folder))
						try: os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder): pass
							else: raise

						#writing the com files
						for file in conf_files: # check conf_file exists, parse energies and then write dft input
							if os.path.exists(file):
								if args.verbose: log.write("   -> Converting from {}".format(file))
								energies = read_energies(file,log)
								name = os.path.splitext(file)[0]
								write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args,log,charge_data)

					else: # else create the directory for optimizations
						folder = g_dir + '/' + str(lot) + '-' + str(bs)
						log.write("\no  PREPARING GAUSSIAN INPUTS in {}".format(folder))
						try: os.makedirs(folder)
						except OSError:
							if  os.path.isdir(folder): pass
							else: raise

						#writing the com files
						print(conf_files)
						for file in conf_files: # check conf_file exists, parse energies and then write dft input

							if os.path.exists(file):
								if args.verbose: log.write("   -> Converting from {}".format(file))
								energies = read_energies(file,log)
								name = os.path.splitext(file)[0]

								write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, args,log,charge_data)

	#moving files arefter compute and write_gauss or only after compute
	#moving all the sdf files to a separate folder after writing gaussian files
	src = os.getcwd()
	if args.xtb == True:
		all_xtb_conf_files = glob.glob('*_xtb.sdf')
		destination_xtb = src +'/xTB_minimised_generated_SDF_files'
		for file in all_xtb_conf_files:
			try:
				os.makedirs(destination_xtb)
				shutil.move(os.path.join(src, file), os.path.join(destination_xtb, file))
			except OSError:
				if  os.path.isdir(destination_xtb):
					shutil.move(os.path.join(src, file), os.path.join(destination_xtb, file))
				else:
					raise
	elif args.ANI1ccx == True:
		all_ani_conf_files = glob.glob('*_ani.sdf')
		destination_ani = src +'/ANI1ccx_minimised_generated_SDF_files'
		for file in all_ani_conf_files:
			try:
				os.makedirs(destination_ani)
				shutil.move(os.path.join(src, file), os.path.join(destination_ani, file))
			except OSError:
				if  os.path.isdir(destination_ani):
					shutil.move(os.path.join(src, file), os.path.join(destination_ani, file))
				else:
					raise
	all_name_conf_files = glob.glob('*.sdf')
	destination_rdkit = 'RDKit_generated_SDF_files'
	for file in all_name_conf_files:
		try:
			os.makedirs(destination_rdkit)
			shutil.move(os.path.join(src, file), os.path.join(destination_rdkit, file))
		except OSError:
			if  os.path.isdir(destination_rdkit):
				shutil.move(os.path.join(src, file), os.path.join(destination_rdkit, file))
			else:
				raise

	if args.analysis == True:
		#adding in for general analysis
		#need to specify the lot, bs as arguments for each analysis
		if args.path == '':
			log_files = glob.glob('*.log')
			w_dir = os.getcwd()
			w_dir_fin = w_dir+'/Finished'
			output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,log)
		#taking the coorct path
		else:
			# Sets the folder and find the log files to analyze
			for lot in args.level_of_theory:
				for bs in args.basis_set:
					for bs_gcp in args.basis_set_genecp_atoms:
						w_dir = args.path + str(lot) + '-' + str(bs) +'/'
						#check if New_Gaussian_Input_Files folder exists
						w_dir = check_for_final_folder(w_dir,log)
						#assign the path to the finished directory.
						w_dir_fin = args.path + str(lot) + '-' + str(bs) +'/finished'
						#log.write(w_dir)
						os.chdir(w_dir)
						log.write(w_dir)
						log_files = glob.glob('*.log')
						output_analyzer(log_files, w_dir, lot, bs, bs_gcp, args, w_dir_fin,log)

	#adding the part to check for resubmission of the newly created gaussian files.
	if args.qsub == True:
		#chceck if ech level of theory has a folder New gaussin FILES
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'
				#check if New_Gaussian_Input_Files folder exists
				w_dir = check_for_final_folder(w_dir,log)
				os.chdir(w_dir)
				cmd = args.submission_command + ' *.com'
				if args.qsub == True:
					os.system(cmd)

	#once all files are finished are in the Finished folder
	if args.dup == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				try:
					log_files = glob.glob('*.log')
					if len(log_files) != 0:
						val = ' '.join(log_files)
						dup_calculation(val,w_dir,args,log)
					else:
						log.write(' Files for are not there!')
					#log.write(log_files)
				except:
					pass

	#once all files are finished are in the Finished folder
	if args.boltz == True:
		# Sets the folder and find the log files to analyze
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) +'/'+'finished'
				os.chdir(w_dir)
				#can change molecules to a range as files will have codes in a continous manner
				for i in range(args.maxnumber):
					#grab all the corresponding files make sure to renamme prefix when working with differnet files
					try:
						log_files = glob.glob('RE' + '_' + str(i)+'_'+'confs_low.log')
						if len(log_files) != 0:
							val = ' '.join(log_files)
							boltz_calculation(val,i,log)
						else:
							log.write(' Files for {} are not there!'.format(i))
						#log.write(log_files)
					except:
						pass

	if args.combine == True:
		#combines the files and gives the boltzmann weighted energies
		for lot in args.level_of_theory:
			for bs in args.basis_set:
				w_dir = args.path + str(lot) + '-' + str(bs) + '/finished'
				os.chdir(w_dir)
				#read the csv log_files
				csv_files = glob.glob('Goodvibes*.csv')
				combine_files(csv_files, lot, bs, args,log)

if __name__ == "__main__":
	main()
