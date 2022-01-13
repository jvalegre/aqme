#########################################################.
#        This file stores all the functions             #
#          used in writing SDF and COM files,           #
#              as well as the logger and                #
#                 yaml file importer                    #
#########################################################.

import os
import sys
import glob
import pandas as pd
from aqme.utils import load_from_yaml, Logger
from aqme.argument_parser import set_options


# template classes
class qprep():

	"""
	Class containing all the functions from the QPREP module related to Gaussian input files.

	Parameters
	----------
	mol : mol object with 3D coordinates
		Optionally, starts from a mol object to prepare the input QM file
	destination : str
		Directory to create the input file
	molecule : str
		Name of the input file (without the extension i.e. NAME for NAME.com)
	charge : int
		Charge of the calculations used in the following input files
	mult : int
		Multiplicity of the calculations used in the following input files
	chk : bool
		True to include the initial %chk line in Gaussian
	mem : str
		Memory used in the input file. Formats: GB, MB or MW (i.e. 4GB, 800MB or 16MW).
	nprocs : int
		Number of processors used in the input file
	suffix : str
		Suffix for the new input files
	atom_types : list of strings
		(If no mol is used) List containing the atoms of the system
	cartesians : list of lists
		(If no mol is used) Cartesian coordinates used for further processing
	bs_gen : str
		Basis set for the Gen(ECP) section
	bs : str
		Basis set for regular atoms when Gen(ECP) is used
	gen_atoms : str
		Atoms considered for Gen(ECP). Format: ATOM1,ATOM2,ATOM3... (i.e. C,H,O)
	program : str
		Target program to generate input files. Options: 'gaussian' and 'orca'
	qm_input : string
		Keyword line of the input file
	qm_end : str
		Final line(s) of the input file
	yaml_file : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""
	
	def __init__(self, mol=None, destination=os.getcwd(), molecule='', charge=0, mult=1, chk=False,
				mem='8GB', nprocs=4, suffix='', atom_types=[], cartesians=[], bs_gen='', bs='', gen_atoms=[],
				program='gaussian',	qm_input='', qm_end='', yaml_file=None, **kwargs):
				
		if mol != None:
			for _,atom in enumerate(mol.GetAtoms()):
				atom_types.append(atom.GetSymbol())

			cartesians = mol.GetConformers()[0].GetPositions()

		self.destination = destination
		self.molecule = molecule
		self.chk = chk
		self.mem = mem
		self.nprocs = nprocs
		self.suffix = suffix
		self.charge = charge
		self.mult = mult
		self.n_atoms = len(atom_types)
		self.atom_types = atom_types
		self.cartesians = cartesians
		self.qm_input = qm_input
		self.bs_gen = bs_gen
		self.bs = bs
		self.gen_atoms = gen_atoms
		self.qm_end = qm_end
		self.program = program

		if 'options' in kwargs:
			self.args = kwargs['options']
		else:
			self.args = set_options(kwargs)

		self.args.varfile = yaml_file

		if yaml_file is not None:
			self.args, self.log = load_from_yaml(self.args, self.log)

		# start a log file to track the QCORR module
		if not os.path.isdir(self.destination+'/dat_files/'):
			os.makedirs(self.destination+'/dat_files/')
		self.log = Logger(self.destination+'/dat_files/pyCONFORT','QPREP')
		self.log.write("\no  Creating input files in {}\n".format(self.destination))

		if self.qm_input == '':
			self.log.write('x  No keywords line was specified! (qm_input=KEYWORDS_LINE).')
			sys.exit('x  No keywords line was specified! (qm_input=KEYWORDS_LINE).')

		if molecule == '':
			self.log.write('x  No name was specified! (molecule=NAME).')
			sys.exit('x  No name was specified! (molecule=NAME).')
		
		self.write()
		
	def get_header(self):
		'''
		Gets the part of the input file above the molecular coordinates.
		'''

		txt = ''

		if self.program.lower() == 'gaussian':
			if self.chk: 
				txt += f'%chk={self.molecule}.chk\n'
			txt += f'%nprocshared={self.nprocs}\n'
			txt += f'%mem={self.mem}\n'
			txt += f'# {self.qm_input}'
			txt += f'\n\n'
			txt += f'{self.molecule}\n\n'
			txt += f'{self.charge} {self.mult}\n'

		elif self.program.lower() == 'orca':
			txt += f'# {self.molecule}\n'
			if self.mem.find('GB'):
				mem_orca = int(self.mem.split('GB')[0])*1000
			elif self.mem.find('MB'):
				mem_orca = self.mem.split('MB')[0]
			elif self.args.mem.find('MW'):
				mem_orca = self.mem.split('MW')[0]
			txt += f'%maxcore {mem_orca}\n'
			txt += f'%pal nprocs {self.nprocs} end\n'
			txt += f'! {self.qm_input}\n'
			txt += f'* xyz {self.charge} {self.mult}\n'

		return txt

	def get_tail(self):
		'''
		Gets the part of the input file below the molecular coordinates.
		'''

		txt = ''

		if self.program.lower() == 'gaussian':
			if len(self.gen_atoms) > 0:
				# writes part for Gen/GenECP
				ecp_used, ecp_not_used, gen_type  = [],[],'gen'
				if self.qm_input.lower().find('genecp') > -1:
					gen_type = 'genecp'

				for _,element_ecp in enumerate(self.atom_types):
					if element_ecp in self.gen_atoms and element_ecp not in ecp_used:
						ecp_used.append(element_ecp)
					elif element_ecp not in self.gen_atoms and element_ecp not in ecp_not_used:
						ecp_not_used.append(element_ecp)

				if len(ecp_not_used) > 0:
					elements_not_used = ' '.join([f'{sym}' for sym in ecp_not_used])
					txt += f'{elements_not_used} 0\n{self.bs}\n****\n'
				if len(ecp_used) > 0:
					elements_used = ' '.join([f'{sym}' for sym in ecp_used])
					txt += f'{elements_used} 0\n{self.bs_gen}\n****\n'
					
				if gen_type == 'genecp' and len(ecp_used) > 0:
					txt += '\n'
					txt += f'{elements_used} 0\n{self.bs_gen}\n****\n'

				txt += '\n'

			# writes final section if selected
			if self.qm_end != '': 
				txt += f'{self.qm_end}\n\n'

		return txt


	def write(self):
		
		if self.program.lower() == 'gaussian':
			extension = 'com'
		elif self.program.lower() == 'orca':
			extension = 'inp'
		if self.suffix != '':
			comfile = f'{self.destination}/{self.molecule}_{self.suffix}.{extension}'
		else:
			comfile = f'{self.destination}/{self.molecule}.{extension}'
		
		if os.path.exists(comfile):
			os.remove(comfile)
		
		header = self.get_header()
		tail = self.get_tail()

		fileout = open(comfile, "w")
		fileout.write(header)
		
		for atom in range(0,self.n_atoms):
			fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(self.atom_types[atom], self.cartesians[atom][0], self.cartesians[atom][1], self.cartesians[atom][2]))
			if atom != self.n_atoms-1:
				fileout.write("\n")

		if self.program.lower() == 'gaussian':
			fileout.write("\n\n")
		elif self.program.lower() == 'orca':
			fileout.write("\n*")

		fileout.write(tail)
		fileout.close()


# Aux Functions for QM input generation
def get_molecule_list(filepath,lowest_only=False,
						lowest_n=False, energy_threshold=0.0):
	out_molecules = []

	molecules = [mol for mol in pybel.readfile('sdf',filepath)]
	energies = [mol.energy for mol in molecules]
	min_energy = energies[0]
	for mol,energy in zip(molecules,energies):
		is_in_threshold = energy - min_energy < energy_threshold
		title,i = mol.title.strip().rsplit(maxsplit=1)
		mol.title = f"{title} {int(i):03d}"
		if lowest_n and is_in_threshold:
			out_molecules.append(mol)
		elif lowest_n:
			break
		elif lowest_only:
			out_molecules.append(mol)
			break
		else:
			out_molecules.append(mol)

	return out_molecules


def load_charge_data(filepath,backup_files):
	#read in dup_data to get the overall charge of MOLECULES
	invalid_files = []
	try:
		charge_data = pd.read_csv(filepath, usecols=['Molecule','Overall charge'])
	except:
		charge_data = pd.DataFrame()
		for i,sdf_file in enumerate(backup_files):
			if not(Path(sdf_file).exists()):
				invalid_files.append(sdf_file)
				maxsplit = 1
				if 'filter' in sdf_file:
					maxsplit += 1
				name = sdf_file.rsplit('_',maxsplit[0])
				charge = 'Invalid'
			else:
				mol = next(pybel.readfile(sdf_file))
				name = mol.title.split(maxsplit=1)[0]
				charge = mol.data['Real charge']
			charge_data.at[i,'Molecule'] = name
			charge_data.at[i,'Overall charge'] = charge
	return charge_data,invalid_files


# MAIN QPREP FUNCTION
def qprep_main(w_dir_initial,args,log):

	if len(args.geom_rules) >= 1:
		conf_files =  glob.glob('*_rules.sdf')
	# define the SDF files to convert to COM Gaussian files
	elif args.CMIN == 'xtb':
		conf_files =  glob.glob('*_xtb.sdf')
	elif args.CMIN=='ani':
		conf_files =  glob.glob('*_ani.sdf')
	elif args.CSEARCH=='rdkit':
		conf_files =  glob.glob('*_rdkit.sdf')
	elif args.CSEARCH=='summ':
		conf_files =  glob.glob('*_summ.sdf')
	elif args.CSEARCH=='fullmonte':
		conf_files =  glob.glob('*_fullmonte.sdf')
	else:
		conf_files =  glob.glob('*.sdf')

	if args.com_from_xyz:
		xyz_files =  glob.glob('*.xyz')
		for file in xyz_files:
			mol = next(pybel.readfile('xyz',file))
			stem = Path(file).stem
			mol.write('sdf',f'{stem}.sdf')
		conf_files =  glob.glob('*.sdf')

	if not conf_files:
		log.write('\nx  No SDF files detected to convert to gaussian COM files')
		return

	# names for directories created
	if args.QPREP == 'gaussian':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/G16')
	elif args.QPREP == 'orca':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/ORCA')
	elif args.QPREP == 'turbomole':
		qm_folder = Path(f'{w_dir_initial}/QMCALC/TURBOMOLE')

	csv_name = args.input.split('.')[0]
	csv_file = f"{w_dir_initial}/CSEARCH/csv_files/{csv_name}-CSEARCH-Data.csv"
	charge_data, invalid_files = load_charge_data(csv_file,conf_files)

	# remove the invalid files and non-existing files
	accept_file = lambda x: x not in invalid_files and Path(x).exists()
	conf_files = [file for file in conf_files if accept_file(file) ]

	# Prepare the list of molecules that are to be written
	molecules = []
	for file in conf_files:
		filepath = f'{w_dir_initial}/{file}'
		new_mols = get_molecule_list(filepath,
							lowest_only=args.lowest_only,lowest_n=args.lowest_n,
							energy_threshold=args.energy_threshold_for_gaussian)
		molecules.extend(new_mols)

		# this variable keeps track of folder creation
		qm_folder.mkdir(parents=True,exist_ok=True)

		# writing the com files
		# check conf_file exists, parse energies and then write DFT input
		write_qm_input_files(args, qm_folder, molecules, 
							charge_data, None, args.mult, 'mol', None, None, None, args.qm_keywords)
