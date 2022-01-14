#########################################################.
#        This file stores all the functions             #
#          used in writing SDF and COM files,           #
#              as well as the logger and                #
#                 yaml file importer                    #
#########################################################.

import subprocess
import os
import sys
import glob
import shlex
import pandas as pd
from aqme.turbomole import TurbomoleInput
from aqme.utils import periodic_table,load_from_yaml,Logger
from aqme.argument_parser import set_options

# Hotfixes
def fix_obabel_isotopes(gjf_lines):
	"""
	Ensures that gaussian input lines with isotopes ('H  (Iso=2)  x  y  z')
	have the appropiate format for gaussian inputs ('H(Iso=2)  x  y  z').
	Modifies in-place the input list.

	Parameters
	----------
	gjf_lines : list
		list of strings corresponding to a gaussian input file
	"""
	for i,line in enumerate(gjf_lines):
		if '(Iso' in line:
			assert len(line.split()) == 5
			gjf_lines[i] = ''.join(line.split(maxsplit=1))

# Basis set class
class BasisSet(object):
	"""
	Representation of the neccesary information related with basis sets. May
	contain information about auxiliary basis sets and ECPs.

	Parameters
	----------
	basis : dict, optional
		dictionary with element symbols as keys and their corresponding basis
		set as value. 'all' and 'default' are accepted as keys, by default None
	auxbasis : dict, optional
		dictionary with element symbols as keys and their corresponding auxiliary
		basis set as value. 'all' and 'default' are accepted as keys, by default
		None
	ecp : dict, optional
		dictionary with element symbols as keys and their corresponding ECP
		as value. 'all' and 'default' are not accepted as keys, by default None.
	"""
	_counter = 0

	def __init__(self,basis=None,auxbasis=None,ecp=None,name=None):
		self.basis = basis
		if basis is None:
			self._basis = None
		else:
			self._basis = dict()
			self._basis.update(self.basis)
		self.auxbasis = auxbasis
		if auxbasis is None:
			self._auxbasis = None
		else:
			self._auxbasis = dict()
			self._auxbasis.update(self.auxbasis)
		self.elements = set()
		self._elements = set() # Hardcopy to allow reset behaviour
		reserved_keys = ['all','default']
		if not(ecp is None) and any(key in reserved_keys for key in ecp):
			raise ValueError("'all' or 'default' key found in the ecp parameter")
		self.ecp = ecp
		items = [basis,auxbasis,ecp]
		for item in items:
			if item is None:
				continue
			keys = [key for key in item if key not in reserved_keys]
			self.elements.update(keys)
			self._elements.update(keys)
		if name is None:
			default_name = f'basisset_{self._counter:02d}'
			name = self.basis('all','') or self.basis('default','')
			if not name:
				name = default_name
				self._counter += 1
		self.name = self._homogenize_name(name)

	def _homogenize_name(self,name):
		"""
		Ensures that the name of the basis set can be used as a folder name.

		Parameters
		----------
		name : str
			initialized name

		Returns
		-------
		str
			modified name
		"""
		if name.find('**') > -1:
			name = name.replace('**','(d,p)')
		elif name.find('*') > -1:
			name = name.replace('*','(d)')

		if str(name).find('/') > -1:
			new_name =  str(name).split('/')[0]
		else:
			new_name = name
		return new_name

	def update_atomtypes(self,molecules):
		"""
		Updates the elements of the basisset with an iterable of molecules.

		Parameters
		----------
		molecules : list
			list of RDKit Molecule objects
		"""
		elements = set()
		for molecule in molecules:
			if isinstance(molecule,rdkit.Chem.rdchem.Mol):
				atoms = [atom for atom in molecule.GetAtoms()]
				syms = [atom.GetSymbol() for atom in atoms]
			elif isinstance(molecule,pybel.Molecule):
				atoms = [atom for atom in molecule.atoms]
				atnums = [atom.OBAtom.GetAtomicNum() for atom in atoms]
				syms = [periodic_table()[atnum] for atnum in atnums]
			else:
				raise ValueError(f'{molecule.__class__} is not a valid molecule object')
			elements.update(syms)
		self.elements.update(elements)
		# Now update the basis,auxbasis and ecp dictionaries
		for element in elements:
			for item in [self.basis,self.auxbasis]:
				if item is None:
					continue
				val = item.get(element,'') or item.get('default','') or item.get('all','')
				item[element] = val

	def clear_atomtypes(self):
		"""
		Resets the state of the Basisset to its initialization state.

		Parameters
		----------
		molecules : list
			list of RDKit Molecule objects
		"""
		self.elements = set()
		self.elements.update(self._elements)
		for attr in ['basis','auxbasis']:
			val = getattr(self,attr)
			if val is None:
				continue
			backup = getattr(self,f'_{attr}')
			new_val = dict()
			new_val.update(backup)
			setattr(self,attr,new_val)

	def write_basisfile(self,filepath):
		"""
		Writes a basis set file that can be read back into a basis set object.

		Parameters
		----------
		filepath : str or Path
			path to the desired basis set file.
		"""
		lines = []
		for attr in ['basis','auxbasis','ecp']:
			mappable = getattr(self,attr)
			if mappable is None:
				continue
			lines.append(f'${attr}')
			k_all = mappable.get('all',None)
			k_default = mappable.get('default',None)
			if k_all is not None:
				value = k_all
				lines.append(f'all {value}')
			elif k_default is not None:
				value = k_default
				lines.append(f'default {value}')
			for element in self.elements:
				item = mappable.get(element,value)
				lines.append(f'{element} {item}')
		lines.append('$end')
		with open(filepath,'w') as F:
			F.write('\n'.join(lines))

	@classmethod
	def from_file(cls,filepath):
		"""
		Reads a file containing information on basis sets, auxiliar basis sets and
		ECPs. The file follows the following format:
		$basis
		sym1 basis1
		sym2 basis2
		...
		$auxbasis
		sym1 auxbasis1
		sym2 auxbasis2
		...
		$ecp
		sym1 ecp1
		sym2 ecp2
		...
		$end

		Parameters
		----------
		filepath : str or Path
			a valid path to a file with the previously specified format.

		Returns
		-------
		dict or None
			basis, auxbasis, ecp
		"""
		regex = r'\$([^\$]*)'
		regex = re.compile(regex)
		name = Path(filepath).stem
		with open(filepath,'r') as F:
			txt = F.read()
		keywords = dict()
		kwblocks = regex.findall(txt)
		for block in kwblocks:
			lines = [line for line in block.split('\n') if line]
			kw_line = lines[0]
			keyword = kw_line.split(' ')[0]
			if keyword == 'end':
				continue
			keywords[keyword] = dict()
			for line in lines[1:]:
				sym,basis = line.split(' ',1)
				keywords[keyword][sym] = basis.strip()
		basis = keywords.get('basis',None)
		auxbasis = keywords.get('auxbasis',None)
		ecp = keywords.get('ecp',None)
		return cls(basis,auxbasis,ecp,name)


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
	atom_types : list of strings
		(If no mol is used) List containing the atoms of the system
	cartesians : list of lists
		(If no mol is used) Cartesian coordinates used for further processing
	bs_gen : bool
		Basis set for the Gen(ECP) section
	qm_input : string
		Keyword line of the input file
	yaml_file : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""
	
	def __init__(self, mol=None, destination=os.getcwd(), molecule='', charge=0,
				mult=1, atom_types=[], cartesians=[], bs_gen='',
				qm_input='', yaml_file=None, **kwargs):
				
		if mol != None:
			for _,atom in enumerate(mol.GetAtoms()):
				atom_types.append(atom.GetSymbol())

			cartesians = mol.GetConformers()[0].GetPositions()

		self.destination = destination
		self.molecule = molecule
		self.charge = charge
		self.mult = mult
		self.n_atoms = len(atom_types)
		self.atom_types = atom_types
		self.cartesians = cartesians
		self.qm_input = qm_input
		self.bs_gen = bs_gen

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
		self.log = Logger(self.destination+'/dat_files/aqme','QPREP')
		self.log.write("\no  Creating input files in {}\n".format(self.destination))

		if self.qm_input == '':
			self.log.write('x  No keywords line was specified! (qm_input=KEYWORDS_LINE).')
			sys.exit('x  No keywords line was specified! (qm_input=KEYWORDS_LINE).')

		if molecule == '':
			self.log.write('x  No name was specified! (molecule=NAME).')
			sys.exit('x  No name was specified! (molecule=NAME).')

		self.write()
		
	def get_header(self):
		txt = ''
		#Link0
		if self.args.chk: 
			txt += f'%chk={self.molecule}.chk\n'
		txt += f'%nprocshared={self.args.nprocs}\n'
		txt += f'%mem={self.args.mem}\n'
		txt += f'# {self.qm_input}'
		return txt

	def get_tail(self):
		txt = ''
		# writes part for Gen/GenECP
		ecp_used, ecp_not_used, gen_atoms = [],[],[]
		if self.args.gen_atoms != None:
			gen_atoms = self.args.gen_atoms.split(',')
		elif self.args.gen_atoms != None:
			gen_atoms = self.args.gen_atoms.split(',')

		for _,element_ecp in enumerate(gen_atoms):
			if element_ecp in self.atom_types and element_ecp not in ecp_used:
				ecp_used.append(element_ecp)
			elif element_ecp not in self.atom_types and element_ecp not in ecp_used:
				ecp_not_used.append(element_ecp)

		if len(ecp_used) > 0:
			elements_used = ' '.join([f'-{sym}' for sym in ecp_used])
			txt += f'{elements_used} 0\n{self.bs_gen}\n****\n'
		if len(ecp_not_used) > 0:
			elements_not_used = ' '.join([f'-{sym}' for sym in ecp_not_used])
			txt += f'{elements_not_used} 0\n{self.bs_gen}\n****\n'
		if self.args.gen_atoms != None and len(ecp_used) > 0:
			txt += '\n'
			txt += f'{elements_used} 0\n{self.bs_gen}\n****\n'
		if self.args.gen_atoms != None or self.args.gen_atoms != None:
			txt += '\n'

		# writes final section if selected
		if self.args.qm_input_end != '': 
			txt += f'{self.args.qm_input_end}\n\n'

		return txt


	def write(self):

			comfile = f'{self.destination}/{self.molecule}.com'
			
			if os.path.exists(comfile):
				os.remove(comfile)

			header = self.get_header()
			tail = self.get_tail()

			fileout = open(comfile, "w")
			fileout.write(header)
			fileout.write('\n\n')
			fileout.write(self.molecule+"\n\n")
			fileout.write(str(self.charge)+' '+str(self.mult)+'\n')
			
			for atom in range(0,self.n_atoms):
				fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(self.atom_types[atom], self.cartesians[atom][0], self.cartesians[atom][1], self.cartesians[atom][2]))
				fileout.write("\n")
			fileout.write("\n")

			fileout.write(tail)
			fileout.close()

class OrcaTemplate(object):
	"""
	[summary]

	Parameters
	----------
	basisset : [type], optional
		[description], by default None
	functional : [type], optional
		[description], by default None
	memory : str, optional
		[description], by default '24GB'
	nprocs : int, optional
		[description], by default 12
	extra_commandline : [type], optional
		[description], by default None
	solvation : str, optional
		[description], by default 'gas_phase'
	solvent : str, optional
		[description], by default ''
	extra_cpcm_input : [type], optional
		[description], by default None
	scf_iters : int, optional
		[description], by default 500
	mdci : list, optional
		additional input for the MDCI (Matrix driven correlation) section of
		ORCA input, by default None
	print_mini : bool, optional
		[description], by default True
	"""

	def __init__(self,basisset=None,functional=None,memory='24GB',nprocs=12,
				extra_commandline=None,solvation='gas_phase',solvent='',
				extra_cpcm_input=None,scf_iters=500,mdci=None,
				print_mini=True):
		self.basisset = basisset
		self.functional = functional
		# calculate memory for ORCA input
		mem_orca = int(memory[:-2])
		is_GB = 'gb' in memory.lower()
		if is_GB:
			mem_orca *= 1000
		else:
			# assume MB
			pass
		self.memory = mem_orca # in MB
		self.nprocs = nprocs
		if extra_commandline is None:
			extra_commandline = ''
		self.extra_commandline = extra_commandline
		self.solvation = solvation
		self.solvent = solvent
		if extra_cpcm_input is None:
			extra_cpcm_input = []
		self.extra_cpcm_input = extra_cpcm_input
		self.scf_iters = scf_iters
		if mdci is None:
			mdci = ['Density None']
		self.mdci = mdci
		self.print_mini = print_mini

	def to_dict(self):
		kwargs = dict()
		kwargs['functional'] = self.functional
		if self.basisset is not None:
			kwargs['basis'] = self.basisset.basis
			kwargs['auxbasis'] = self.basisset.auxbasis
			kwargs['ecp'] = self.basisset.ecp

		kwargs['ecp_list'] = self.aux_atoms_orca
		kwargs['bs_gcp'] = self.aux_genecp_bs
		kwargs['bs_gcp_fit'] = self.aux_fit_gen_atoms
		kwargs['orca_aux_section'] = self.enable_aux_section
		kwargs['memory'] = self.memory
		kwargs['nprocs'] = self.nprocs
		kwargs['extra_input'] = self.extra_input
		kwargs['solvation'] = self.solvation
		kwargs['solvent'] = self.solvent
		kwargs['cpcm_input_orca'] = self.cpcm_input_orca
		kwargs['scf_iters_orca'] = self.scf_iters_orca
		kwargs['orca_mdci'] = self.mdci_orca
		kwargs['print_mini'] = self.print_mini

		return kwargs

	@property
	def enable_aux_section(self):
		return bool(self.basisset.elements) or (self.basisset.ecp is not None)
	@property
	def basis(self):
		if self.enable_aux_section:
			return self.basisset.basis.get('default','')
		basis_all = self.basisset.basis.get('all','')
		if basis_all:
			return basis_all
		raise RuntimeError("Could not determine the basis")

	@classmethod
	def from_args(cls,args):
		new = cls()
		# Standard stuff
		new.memory = args.mem
		new.nprocs = args.nprocs
		new.extra_commandline = args.qm_input
		new.solvation = args.solvent_model
		new.solvent = args.solvent_name
		new.extra_cpcm_input = args.cpcm_input
		new.scf_iters = args.orca_scf_iters
		new.mdci = args.mdci_orca
		new.print_mini = args.print_mini_orca

		return new

	def get_aux_section(self):
		txt = '%basis\n'
		# basis set block
		def_basis = self.basis
		items = []
		for sym in self.basisset.elements:
			basis = self.basisset.basis.get(sym,'')
			if basis and basis == def_basis:
				continue
			items.append((sym,basis))
		if items:
			items = sorted(items,key=itemgetter(1))
			for basis,pairs in groupby(items,key=itemgetter(1)):
				elements = ' '.join([f'{sym}' for sym,_ in pairs])
				txt += f'NewGTO {elements} {basis} end\n'
		# Auxiliary Basis set block
		# TODO
		if self.basisset.ecp is None:
			txt += 'end'
			return txt
		# ECP block
		items = []
		for sym in self.basisset.ecp:
			ecp = self.basisset.ecp.get(sym,'')
			if basis and basis == def_basis:
				continue
			items.append((sym,ecp))
		if items:
			items = sorted(items,key=itemgetter(1))
			for ecp,pairs in groupby(items,key=itemgetter(1)):
				elements = ' '.join([f'{sym}' for sym,_ in pairs])
				txt += f'NewECP {elements} {ecp} end\n'
		txt += 'end'
		return txt
	def get_solvation_section(self):
		txt = ''
		if self.solvation.lower() == 'smd':
			txt += '%cpcm\n'
			txt += 'smd true\n'
			txt += f'SMDsolvent "{self.solvent}"\n'
		elif self.solvation.lower() == 'cpcm':
			txt += f'! CPCM({self.solvent})\n'
			if self.cpcm_input_orca != 'None':
				txt += '%cpcm\n'
				for cpcm_line in self.cpcm_input_orca:
					txt += f'{cpcm_line}\n'
		txt += 'end'
		return txt
	def get_command_line(self):
		commandline = ' '.join(f'{self.basis} {self.functional}'.split())
		if self.extra_commandline != 'None':
			commandline += f' {self.extra_commandline}'
		return commandline

	def get_header(self,title):
		txt =  f'# {title}\n'
		txt +=  '# Memory per core\n'
		txt += f'%maxcore {self.args.memory}\n'
		txt +=  '# Number of processors\n'
		txt += f'%pal nprocs {self.args.nprocs} end\n'
		txt += f'! {self.get_command_line()}\n'

		sections = []
		if self.enable_aux_section:
			sections.append(self.get_aux_section())

		if self.solvation != 'gas_phase':
			sections.append(self.get_solvation_section())

		sections.append(f'%scf maxiter {self.scf_iters}\nend')

		if self.mdci:
			mdci_lines = ['% mdci',] + self.mdci + ['end',]
			sections.append('\n'.join(mdci_lines))

		if self.print_mini:
			mini_lines = [  '%output',
							'printlevel mini',
							'print[ P_SCFInfo ] 1',
							'print[ P_SCFIterInfo ] 1',
							'print[ P_OrbEn ] 0',
							'print[ P_Cartesian ] 0',
							'end',
							'%elprop',
							'Dipole False',
							'end']
			sections.append('\n'.join(mini_lines))
		txt += '\n'.join(sections)
		txt += '\n'
		return txt
	def write(self,destination,molecule):
		stem = molecule.title.lstrip().replace(' ','_')
		infile = f'{destination}/{stem}.inp'
		header = self.get_header(molecule.title)
		xyz_lines = molecule.write('orcainp').split('\n')[3:]

		with open(infile,'w') as F:
			F.write(header)
			F.write('\n'.join(xyz_lines))

		return infile

class TurbomoleTemplate(object):
	def __init__(self,basisset=None,functional=None,dispersion='off',
				grid='m4',epsilon='gas',maxcore=200,ricore=200,cavity='none'):
		self.basisset = basisset
		self.functional = functional
		self.dispersion =  dispersion
		self.grid =  grid
		self.epsilon =  epsilon
		self.maxcore = maxcore
		self.ricore = ricore
		self.cavity =  cavity

	def to_dict(self):
		kwargs = dict()
		kwargs['functional'] = self.functional
		if self.basisset is not None:
			kwargs['basis'] = self.basisset.basis
			kwargs['auxbasis'] = self.basisset.auxbasis
			kwargs['ecp'] = self.basisset.ecp

		kwargs['dispersion'] = self.dispersion
		kwargs['grid'] = self.grid
		kwargs['epsilon'] = self.epsilon
		kwargs['maxcore'] = self.maxcore
		kwargs['ricore'] = self.ricore
		kwargs['cavity'] = self.cavity
		return kwargs

	@classmethod
	def from_args(cls,args):
		new = cls()
		new.dispersion = args.tmdispersion
		new.grid = args.tmgrid
		new.epsilon = args.tmepsilon
		new.maxcore = args.tmmaxcore
		new.ricore = args.tmricore
		new.cavity = args.tmcavity

		return new
	def write(self,destination,molecule):
		name,i = molecule.title.strip().rsplit(maxsplit=1)
		stem = f'{name}_{i}'
		folder = Path(f'{destination}/{stem}')
		assert folder.parent.exists()
		folder.mkdir(exist_ok=True,parents=False)
		xyz_file = str(folder/f'{stem}.xyz')
		molecule.write('xyz',xyz_file)
		coord_file = str(folder/f'coord')
		molecule.write('tmol',coord_file)
		kwargs = self.to_dict()
		t_input = TurbomoleInput(folder,charge=molecule.charge, multiplicity=molecule.spin,
		title=molecule.title,**kwargs)
		t_input.generate()
		return folder

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
		qm_folder = Path(f'{w_dir_initial}/QCALC/G16')
	elif args.QPREP == 'orca':
		qm_folder = Path(f'{w_dir_initial}/QCALC/ORCA')
	elif args.QPREP == 'turbomole':
		qm_folder = Path(f'{w_dir_initial}/QCALC/TURBOMOLE')

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

# FUNCTIONS FROM THE OLD CODE JUST TO MAKE RAUL'S VERSION WORK!

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

# MAIN FUNCTION TO CREATE GAUSSIAN JOBS
def write_gaussian_input_file(file, name, lot, bs, bs_gcp, energies, keywords_line, args, log, charge_data, w_dir_initial):

	# get the names of the SDF files to read from depending on the optimizer and their suffixes. Also, get molecular charge
	if module_type == 'qprep':
		charge_com = get_name_and_charge(name,charge_data)
	elif module_type == 'qcorr':
		charge_com = charge_qcorr

	if charge_com != 'Invalid':

		com = '{0}_.com'.format(name)
		com_low = '{0}_low.com'.format(name)

		header = header_com(name, args.mem, args.nprocs, args, keywords_line)

		# defining path to place the new COM files
		if module_type == 'qprep':
			program_input == args.qprep
			variable_path = '/QCALC'
		elif module_type == 'qcorr':
			program_input == args.qcorr
			variable_path = '/QCORR'

			if program_input == 'gaussian':
				path_write_input_files = variable_path+'/G16/'
			elif program_input == 'orca':
				path_write_input_files = variable_path+'/ORCA/'

			os.chdir(w_dir_initial+path_write_input_files)
			convert_sdf_to_com(w_dir_initial,file,com,com_low,energies,header,args,log)


		com_files = glob.glob('{0}_*.com'.format(name))

		for file in com_files:
			#patch for Isotopes
			cmd = shlex.split(r"sed -i 's/\([a-zA-Z]\{1,3\}\)[^(]*\([(]Iso\)/\1\2/g' "+file)
			subprocess.call(cmd)
			ecp_list,ecp_gen_atoms,ecp_gen_atoms = [],False,False
			read_lines = open(file,"r").readlines()
			rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)
			read_lines = open(file,"r").readlines()

			# Detect if there are atoms to use genecp or not (to use gen)
			atom_types = []
			for i in range(4,len(read_lines)):
				if read_lines[i].split(' ')[0] not in atom_types and read_lines[i].split(' ')[0] in possible_atoms:
					atom_types.append(read_lines[i].split(' ')[0])

			# write genecp/gen part
			type_gen = 'qprep'
			if args.QPREP == 'gaussian':
				# define genecp/gen atoms
				ecp_list,ecp_gen_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(atom_types,args,'analysis','gaussian')

				#error if both genecp and gen are
				if ecp_gen_atoms and ecp_gen_atoms:
					sys.exit("x  ERROR: Can't use Gen and GenECP at the same time")
				fileout = open(file, "a")

				if genecp != 'None':
					write_genecp(atom_types,type_gen,fileout,genecp,ecp_list,ecp_gen_atoms,ecp_gen_atoms,bs,lot,bs_gcp,args,w_dir_initial,path_write_input_files)

				if args.qm_input_end != 'None':
					fileout.write(args.qm_input_end+'\n\n')

				fileout.close()

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
				ecp_list,ecp_gen_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(atom_types,args,'analysis','orca')

				rename_file_name = rename_file_and_charge_chk_change(read_lines,file,args,charge_com)

				#create input file
				orca_file_gen(read_lines,rename_file_name,bs,lot,genecp,args.aux_atoms_orca,args.aux_genecp_bs,args.aux_fit_gen_atoms,charge_com,args.mult,orca_aux_section,args,args.qm_input,args.solvent_model,args.solvent_name,args.cpcm_input,args.orca_scf_iters,args.mdci_orca,args.print_mini_orca)

			# submitting the input file on a HPC
			if args.qsub:
				cmd_qsub = [args.submission_command, rename_file_name]
				subprocess.call(cmd_qsub)

	os.chdir(w_dir_initial)

	return keywords_line


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
