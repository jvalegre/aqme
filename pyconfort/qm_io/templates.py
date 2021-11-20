from operator import itemgetter
from itertools import groupby
from pathlib import Path

from .gaussian import fix_obabel_isotopes
from .turbomole import TurbomoleInput

class GaussianTemplate(object): 
	def __init__(self,basisset=None,functional=None,memory='24GB',nprocs=12,
				commandline=None,solvation='gas_phase',solvent='',
				max_opt_cycles=500,dispersion=None,enable_chk=False, 
				optimization=True,frequencies=True,calcfc=False,
				last_line_for_input=''): 
		self.basisset = basisset
		self.functional = functional
		#Link0
		self.enable_chk = enable_chk
		self.memory = memory
		self.nprocs = nprocs
		# opt options
		self.optimization = optimization
		self.calcfc = calcfc
		self.max_opt_cycles = max_opt_cycles
		# freq options
		self.frequencies = frequencies
		# scrf options
		self.solvation = solvation
		self.solvent = solvent
		# dispersion
		if dispersion is None: 
			dispersion = ''
		self.dispersion = dispersion
		# other commands
		if commandline is None:
			commandline = ''
		self._commandline = commandline
		# other tail
		self.last_line_for_input = last_line_for_input

	def to_dict(self):
		kwargs = dict()
		kwargs['functional'] = self.functional
		if self.basisset is not None: 
			kwargs['basis'] = self.basisset.basis
			kwargs['auxbasis'] = self.basisset.auxbasis
			kwargs['ecp'] = self.basisset.ecp
			
		kwargs['memory'] = self.memory
		kwargs['nprocs'] = self.nprocs
		kwargs['extra_input'] = self.extra_input
		kwargs['enable_chk'] = self.enable_chk
		kwargs['frequencies'] = self.frequencies
		kwargs['dispersion'] = self.dispersion
		kwargs['calcfc'] = self.calcfc
		kwargs['solvation'] = self.solvation
		kwargs['solvent'] = self.solvent
		kwargs['last_line_for_input'] = self.last_line_for_input
		kwargs['max_opt_cycles'] = self.max_opt_cycles
		return kwargs

	@property
	def genecp(self):
		ecp_genecp_atoms = self.basisset.ecp is not None
		ecp_gen_atoms = bool(self.basisset.elements)
		if ecp_genecp_atoms:
			return 'genecp'
		elif ecp_gen_atoms:
			return 'gen'
		else:
			return ''
	@property
	def basis(self):
		if self.genecp:
			return self.genecp
		basis_all = self.basisset.basis.get('all','')
		if basis_all:
			return basis_all
		basis_default = self.basisset.basis.get('default','')
		if basis_default:
			return basis_default
		raise RuntimeError("""Could not determine basis keyword: 
		No 'all' nor 'default' basis set specified and no element specified""")
	
	@property
	def commandline(self): 
		if self._commandline:
			return self._commandline
		return self.get_commandline()

	def get_commandline(self):
		opt_k = ''
		if self.optimization: 
			opt_sub = f'(maxcycles={self.max_opt_cycles}'
			if self.calcfc:
				opt_sub += ',calcfc'
			opt_sub += ')'
			opt_k = f'opt={opt_sub}'
		freq_k = ''
		if self.frequencies: 
			freq_k = 'freq=noraman'
		disp_k = ''
		if self.dispersion:
			disp_k = f'empiricaldispersion={self.dispersion}'
		scrf_k = ''
		if self.solvation != 'gas_phase': 
			scrf_k = f'scrf=({self.solvation},solvent={self.solvent})'
		commandline = f'{self.functional} {self.basis} {opt_k} {freq_k} {disp_k} {scrf_k}'
		return ' '.join(commandline.split())

	@classmethod
	def from_args(cls,args):
		new = cls()

		# Basis set creation
		new.genecp_atoms_sp = args.genecp_atoms_sp
		new.gen_atoms_sp = args.gen_atoms_sp
		new.genecp_atoms = args.genecp_atoms
		new.gen_atoms = args.gen_atoms
		
		#Standard stuff
		new.memory = args.mem
		new.nprocs = args.nprocs
		new.enable_chk = args.chk
		new.extra_input = args.set_input_line
		new.frequencies = args.frequencies
		if args.empirical_dispersion != None: 
			new.dispersion = args.empirical_dispersion
		new.calcfc = args.calcfc
		new.solvation = args.solvent_model
		new.solvent = args.solvent_name
		new.last_line_of_input = args.last_line_for_input
		new.max_opt_cycles = args.max_cycle_opt
		return new

	def get_header(self,filename='',title=''):
		txt = ''
		#Link0
		if self.enable_chk: 
			txt += f'%chk={filename}.chk\n'
		txt += f'%nprocshared={self.nprocs}\n'
		txt += f'%mem={self.memory}\n'
		txt += f'# {self.commandline}'
		return txt
	def get_tail(self):
		txt = ''
		# write basis set section
		basisset = self.basisset.basis
		if self.genecp:
			def_basis = basisset.get('all','')
			def_basis = basisset.get('default',def_basis)
			items = [(sym,basisset.get(sym,def_basis)) for sym in self.basisset.elements]
			items = sorted(items,key=itemgetter(1))
			for basis,pairs in groupby(items,key=itemgetter(1)):
				elements = ' '.join([f'-{sym}' for sym,_ in pairs])
				txt += f'{elements} 0\n{basis}\n****\n'
			txt += '\n'
		# write ecp section
		if self.genecp == 'genecp': 
			items = sorted(list(self.basisset.ecp.items()),key=itemgetter(1))
			for basis,pairs in groupby(items,key=itemgetter(1)):
				elements = ' '.join([f'-{sym}' for sym,_ in pairs])
				txt += f'{elements} 0\n{basis}\n****\n'
			txt += '\n'
		if self.last_line_for_input: 
			txt += f'{self.last_line_for_input}\n'
		txt += '\n'
		return txt
	def write(self,destination,molecule): 
		stem = molecule.title.lstrip().replace(' ','_')
		comfile = f'{destination}/{stem}.com'
		
		header = self.get_header(filename=stem)

		lines = molecule.write('gjf',opt=dict(k=header)).split('\n')
		fix_obabel_isotopes(lines)
		tail = self.get_tail()
		lines.append(tail)
		with open(comfile,'w') as F:
			F.write('\n'.join(lines))

		return comfile

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
	optimization : bool, optional
		[description], by default True
	frequencies : bool, optional
		[description], by default True
	"""

	def __init__(self,basisset=None,functional=None,memory='24GB',nprocs=12,
				extra_commandline=None,solvation='gas_phase',solvent='',
				extra_cpcm_input=None,scf_iters=500,mdci=None,
				print_mini=True,optimization=True,frequencies=True):
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
		self.optimization = optimization
		self.frequencies = frequencies
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
		kwargs['bs_gcp'] = self.aux_basis_set_genecp_atoms
		kwargs['bs_gcp_fit'] = self.aux_fit_genecp_atoms
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
		new.extra_commandline = args.set_input_line 
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
		if self.optimization: 
			commandline += ' OPT'
		if self.frequencies: 
			commandline += ' FREQ'
		if self.extra_commandline != 'None':
			commandline += f' {self.extra_commandline}'
		return commandline

	def get_header(self,title):
		txt =  f'# {title}\n'
		txt +=  '# Memory per core\n'
		txt += f'%maxcore {self.memory}\n'
		txt +=  '# Number of processors\n'
		txt += f'%pal nprocs {self.nprocs} end\n'
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
