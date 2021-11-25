from pyconfort.exp_rules import passes_custom_rules
from pyconfort.qcorr_gaussian import *
from pyconfort.utils import periodic_table
import numpy

try:
	import pybel
except ImportError:
	from openbabel import pybel

INPUT_SUFFIXES = ['.com','.gjf']
OUTPUT_SUFFIXES = ['.log','.LOG','.out','.OUT','.json']

class GaussianOutputFile(object): 
	def __init__(self,file): 
		with open(file,'r') as F: 
			lines = F.readlines()
		self.lines = lines
		
		# Start Parsing 
		name,charge,multiplicity = self._get_charge_and_multiplicity()
		self.name = name
		self.charge = charge
		self.multiplicity = multiplicity
		termination,errortype = self._get_termination_type()
		self.termination = termination
		self.errortype = errortype
		rms,rms_line = self._get_rms()
		self.rms = rms
		atoms,cartesians = self._get_last_geometry(before_line=rms_line)
		self.atomtypes = atoms
		self.natoms = len(atoms)
		self.cartesians = cartesians
		self._cartesians = cartesians  # backup in case of modifying the coordinates
		frequencies,reducedmass,forceconst,normalmodes = self._get_frequencies()
		self.freqs = frequencies
		self.nfreqs = len(frequencies)
		self.im_freqs = [i for i,freq in enumerate(frequencies) if freq < 0]
		self.readmass = reducedmass
		self.forceconst = forceconst
		self.normalmode = normalmodes

	def _get_charge_and_multiplicity(self): # TODO This is easier finding the matching input file
		stop_name = 0
		# only for name an and charge
		for i,line in enumerate(self.lines):
			if stop_name == 2:
				break
			# Get the name of the compound (specified in the title)
			if line.find('Symbolic Z-matrix:') > -1:
				name = self.lines[i-2]
				stop_name += 1
			# Determine charge and multiplicity
			if line.find("Charge = ") > -1:
				charge = int(line.split()[2])
				multiplicity = int(line.split()[5].rstrip("\n"))
				stop_name += 1
		return name, charge, multiplicity
	def _get_termination_type(self): 
		# use reversed loops to find type of termination (faster than forward loops)
		lines = reversed(self.lines[-16:])
		termination = 'unfinished'
		errortype = 'unknown'
		for line in lines:
			# Determine the kind of job termination
			if 'Normal termination' in line:
				termination = "normal"
			elif 'Error termination' in line:
				termination = "error"
			elif 'Atomic number out of range' in line:
				errortype = "atomicbasiserror"
			elif 'basis sets are only available' in line:
				errortype = "atomicbasiserror"
			elif ' '.join(['SCF Error',]*8) in line: 
				errortype = "SCFerror"
		return termination,errortype
	def _get_last_geometry(self,before_line=-1):
		if self.termination != 'normal': 
			return [],[]
		
		# Find the starting and ending lines of the last orientation
		stop = -1
		for i,line in reversed(enumerate(self.lines[:before_line+1])):
			if 'Standard orientation' in line or 'Input orientation' in line:
				start = i + 5
				break
			is_ori_end = 'Rotational constants' in line or 'Distance matrix' in line
			if stop < 0 and is_ori_end:
				stop = i - 1
		else:
			return [], []
		# If the orientation is found parse it
		cartesians = []
		atoms = []
		for line in self.lines[start:stop]:
			_, atnum, _, x,y,z = line.strip().split()
			atoms.append(periodic_table.index(int(atnum)))
			cartesians.append(list(map(float,[x,y,z])))
		return atoms,cartesians
	def _get_frequencies(self): 
		for i,line in reversed(enumerate(self.lines)): 
			if 'Harmonic frequencies' in line: 
				start = i + 4
				break
			if '- Thermochemistry -' in line: 
				stop = i - 3
		else:
			return [], [], [], []

		frequencies = []
		reducedmass = []
		forceconst =  []
		normalmodes = []
		_iter = self.lines[start:stop].__iter__()
		keep_reading = True
		while keep_reading:
			freqs,redmass,consts,modes = self._get_frequency_block(_iter)
			frequencies.extend(freqs)
			reducedmass.extend(redmass)
			forceconst.extend(consts)
			normalmodes.extend(modes)
			keep_reading = any(freqs,redmass,consts,modes)
		return frequencies,reducedmass,forceconst,normalmodes
	def _get_frequency_block(self,iterable):
		try:
			# Find the frequencies line
			line = next(iterable)
		except StopIteration: 
			return [], [], [], []
		
		while 'Frequencies --' not in line:
			line = next(line)
		frequencies = [float(i) for i in line.strip().split()[2:]]
		# Get the reduced masses
		line = next(iterable)
		reducedmass = [float(i) for i in line.strip().split()[3:]]
		# Get the force constants
		line = next(iterable)
		forceconst = [float(i) for i in line.strip().split()[3:]]
		# Get the normal modes coordinates
		_ = next(iterable)
		_ = next(iterable)
		normal_coords = []
		for _ in range(len(self.atomtypes)):
			line = next(iterable)
			items = line.strip().split()[2:]
			n = items%3
			assert n in [1,2,3]
			normal_coords.append([items[i:i+3] for i in range(n)])
		# unpack the normal coordinates per frequency
		modes = [i for i in zip(*normal_coords)]

		return frequencies,reducedmass,forceconst,modes
	def _get_rms(self):
		if self.termination == 'normal': 
			return self._get_rms_error_termination()
		else:
			return self._get_rms_normal_termination()
	def _get_rms_normal_termination(self):
		for i,line in reversed(enumerate(self.lines)):
			if 'Cartesian Forces:  Max' in line:
				last_line = i
				try:
					rms = float(line.strip().split()[5])
				except ValueError: 
					rms = 10000
				break
		else:
			return 10000, -1
		return rms,last_line

	def _get_rms_error_termination(self):
		for i,line in reversed(enumerate(self.lines)):
			if 'Cartesian Forces:  Max' in line:
				last_line = i
				try:
					rms = float(line.strip().split()[5])
				except ValueError: 
					rms = 10000
				break
		else:
			return 10000, -1
		return rms,last_line

	def fix_imaginary_frequencies(self,amplitude=0.2): 
		xyz = numpy.array(self.cartesians)
		assert xyz.shape[1] == 3
		for im_freq in self.im_freqs: 
			mode = numpy.array(self.normalmode[im_freq])
			assert mode.shape[1] == 3
			xyz = xyz + amplitude*mode
		cartesians = xyz.tolist()
		return cartesians

	# ADAPT THIS FUNCTION TO THE OTHER XYZ GENERATOR!
	def to_xyz(self,cartesians=None):
		if cartesians is None: 
			cartesians = self.cartesians
		n_atoms = len(cartesians)
		xyz = []
		for sym,(x,y,z) in zip(self.atomtypes,cartesians): 
			xyz.append(f'{sym}    {x: 0.6f}    {y: 0.6f}    {z: 0.6f}')
		xyz_print = str(n_atoms)+'\n'+str(self.name)+'\n'+'\n'.join(xyz)+'\n'
		return xyz_print


### Functions copied from qprep

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
			if atomtype not in ecp_list and atomtype in periodic_table:
				ecp_list.append(atomtype)
			if atomtype in genecp_atoms_include:
				ecp_genecp_atoms = True
			if atomtype in gen_atoms_include:
				ecp_gen_atoms = True
		
		if program_gen == 'turbomole': 
			# Do stuff? 
			pass

	if ecp_gen_atoms:
		genecp = 'gen'
	if ecp_genecp_atoms:
		genecp = 'genecp'

	return ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section

# write genecp/gen part
def write_genecp(type_gen,fileout,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,bs_com,lot_com,bs_gcp_com,args,w_dir_initial,new_gaussian_input_files):

	for _,element_ecp in enumerate(ecp_list):
		if type_gen == 'sp':
			if element_ecp not in (args.genecp_atoms_sp or args.gen_atoms_sp):
				fileout.write(element_ecp+' ')
		else:
			if element_ecp not in (args.genecp_atoms or args.gen_atoms):
				fileout.write(element_ecp+' ')
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
		for _,element_ecp in enumerate(ecp_list):
			if type_gen == 'sp':
				if element_ecp in (args.genecp_atoms_sp or args.gen_atoms_sp):
					fileout.write(element_ecp+' ')
			else:
				if element_ecp in (args.genecp_atoms or args.gen_atoms):
					fileout.write(element_ecp+' ')

		fileout.write('0\n')
		fileout.write(bs_gcp_com+'\n')
		fileout.write('****\n\n')
		if ecp_genecp_atoms:
			for _,element_ecp in enumerate(ecp_list):
				if type_gen == 'sp':
					if element_ecp in args.genecp_atoms_sp:
						fileout.write(element_ecp+' ')
				else:
					if element_ecp in args.genecp_atoms:
						fileout.write(element_ecp+' ')

def orca_file_gen(mol,rename_file_name,bs,lot,
					ecp_list,bs_gcp,bs_gcp_fit,
					orca_aux_section,memory,nprocs,
					extra_input,solvation,solvent_orca,
					cpcm_input_orca,scf_iters_orca,orca_mdci,print_mini):
		
	title = mol.title
	lines = [f'# {title}',
			f'# Memory per core']

	# calculate memory for ORCA input
	mem_orca = int(memory[:-2])
	is_GB = 'gb' in memory.lower() 
	if is_GB:
		mem_orca *= 1000
	else:
		# assume MB
		pass
	lines.append(f'%maxcore {mem_orca}')
	
	lines.append('# Number of processors')
	lines.append(f'%pal nprocs {nprocs} end')

	commandline = f'! {bs} {lot}' 
	if extra_input != 'None':
		commandline += f' {extra_input}'
	lines.append(commandline)

	if orca_aux_section:
		lines.append("%basis")

		ecp_atoms = ' '.join([str(periodic_table.index(atom)) for atom in ecp_list])
		gen_orca_line_1 = f'NewGTO {ecp_atoms}'
		gen_orca_line_2 = f'NewAuxCGTO {ecp_atoms}'

		if len(bs_gcp) != 0:
			gen_orca_line_1 += f' "{bs_gcp[0]}" end'
		if len(bs_gcp_fit) != 0:
			gen_orca_line_2 += f' "{bs_gcp_fit[0]}" end'

		lines.append(gen_orca_line_1)
		lines.append(gen_orca_line_2)
		lines.append('end')

	if solvation != 'gas_phase':
		if solvation.lower() == 'smd':
			lines.append('%cpcm') 
			lines.append('smd true')
			lines.append(f'SMDsolvent "{solvent_orca}"')
		elif solvation.lower() == 'cpcm':
			lines.append(f'! CPCM({solvent_orca})')
			if cpcm_input_orca != 'None':
				lines.append('%cpcm')
				for cpcm_line in cpcm_input_orca:
					lines.append(cpcm_line)
		lines.append('end')
	
	lines.append(f'%scf maxiter {scf_iters_orca}')
	lines.append(f'end')

	if orca_mdci != 'None':
		mdci_lines = ['% mdci',] + orca_mdci + ['end',]
		lines.extend(mdci_lines)
	
	if print_mini:
		mini_lines = [  '% output',
						'printlevel mini',
						'print[ P_SCFInfo ] 1',
						'print[ P_SCFIterInfo ] 1',
						'print[ P_OrbEn ] 0',
						'print[ P_Cartesian ] 0',
						'end',
						'% elprop',
						'Dipole False',
						'end']
		lines.extend(mini_lines)

	# this part gets the coordinates from above
	xyz_lines = mol.write('orcainp').split('\n')[3:]
	lines.extend(xyz_lines)

	# Write to file 
	with open(rename_file_name,'w') as F: 
		F.write('\n'.join(lines))

### End of qprep functions 
### Functions copied from mainf
#creation of csv for QCORR
def creation_of_ana_csv(args):

	columns_list = ['Total files', 'Normal termination', 'Imaginary frequencies', 
					'SCF error', 'Basis set error', 'Other errors', 'Unfinished']
	if args.dup:
		columns_list.append('Duplicates')
	if len(args.exp_rules) >= 1:
		columns_list.append('Exp_rules filter')
	if args.check_geom:
		columns_list.append('Geometry changed')
	ana_data = pd.DataFrame(columns = columns_list)

	return ana_data
#finding the file type to move for analysis
def get_com_or_log_out_files(type,name):
	files = []
	if type =='output':
		formats = ['*.log','*.LOG','*.out','*.OUT','*json']
	elif type =='input':
		formats =['*.com','*.gjf']
	for _,format in enumerate(formats):
		if name is None:
			all_files = enumerate(glob.glob(format))
		else:
			all_files = enumerate(glob.glob(name+format))
		for _,file in all_files:
			if file not in files:
				files.append(file)
	return files
### End of mainf functions

## pybel utilities 
def get_btab(OBMol):
	"""
	Gets the bond table of the specified molecule. Each entry in the bond table
	contains the indices of the atoms that participate in the bond and the bond 
	order. 

	Parameters
	----------
	OBMol : pybel.ob.OBMol
		An openbabel molecule object. 

	Returns
	-------
	list
		list of tuples with atom indices and bond orders. 
	"""

	n_bonds = OBMol.NumBonds()
	btab = [] 

	for i in range(n_bonds):
		bond = OBMol.GetBond(i) 
		at1 = bond.GetBeginAtomIdx()
		at2 = bond.GetEndAtomIdx()
		order = bond.GetBondOrder()
		btab.append((at1,at2,order))

	return btab

# Filters
def passes_geometry_check(probe,target,factor=1.0):

	btab = get_btab(target)

	for i,j,_ in btab: # id1,id2,order
		probe_at1 = probe.atoms[i].OBAtom
		probe_at2 = probe.atoms[j].OBAtom
		probe_bond = probe.GetBond(probe_at1,probe_at2)
		probe_length = probe_bond.GetLength()
		
		target_at1 = target.atoms[i].OBAtom
		target_at2 = target.atoms[j].OBAtom
		target_bond = target.GetBond(target_at1,target_at2)
		target_length = target_bond.GetLength()
		
		if target_length > target_bond: 
			large_change = target_length > factor*probe_length
		else:
			large_change = probe_length > factor*target_length
		if large_change:
			return False
	return True

# File utilities
def read_com_as_xyz(file):
	"""
	Reads a gaussian input file and returns a string in the xyz format that 
	could be written to a valid xyz file. 

	Parameters
	----------
	file : str or pathlib.Path
		path to an existing gaussian input file

	Returns
	-------
	str
		string with the geometry of the molecule in xyz format
	"""
	with open(file) as F: 
		_iter = F.__iter__()
		# Find the first empty line after the command line
		line = next(_iter).strip()
		while line:
			line = next(_iter).strip()
		# Assume 1 line of title
		title = [next(_iter).strip(),]
		# Find next empty line
		line = next(_iter).strip()
		while line: 
			title.append(line)
			line = next(_iter).strip()
		# ignore charge and spin
		_ = next(_iter).strip()
		# Store the coordinate lines until the next empty line
		xyz = []
		line = next(_iter).strip()
		while line:
			xyz.append(line)
			line = next(_iter).strip()
	
	xyz_print = str(len(xyz))+'\n'+str(title)+'\n'+'\n'.join(xyz)+'\n'
	return xyz_print
	

# Aux functions
def check_for_final_folder(w_dir):
	base_folder = Path(w_dir)
	inputs_folder = base_folder/'input_files'
	if not inputs_folder.exists():
		return w_dir, 1
	folders = [item for item in inputs_folder.iterdir() if item.isdir() and 'run_' in item.stem]
	get_number = lambda x: int(x.stem.rsplit('_',1)[1]) 
	last_folder = sorted(folders,key=get_number)[-1]
	last_number = get_number(last_folder)
	return str(last_folder), last_number

def output_analyzer(duplicates,log_files,com_files, w_dir, w_dir_main,lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, log, ana_data, round_num):

	input_route = input_route_line(args)
	finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr = 0,0,0,0,0,0,0,0

	if round_num == 1:
		#moves the comfiles to respective folder
		for file in com_files:
			source = w_dir+'/'+file
			destination = w_dir_main +'/input_files/run_'+str(round_num)
			moving_files(source, destination)

	for file in log_files: # Ideally the in_file and out_file should be iterated in parallel

		filepath = Path(f'{w_dir}/{file}')
		GOF = GaussianOutputFile(filepath)

		# Options: 
		# If it has imaginary frequencies -- do stuff
		# If it finished with an error -- do stuff
		# If terminated normally + no img freqs 
		# Move files to their appropiate places
		# If SP or NCIS requested, create them

		if GOF.im_freqs: 
			cartesians = GOF.fix_imaginary_frequencies()
		else:
			cartesians = GOF.cartesians
		
		xyz_str = GOF.to_xyz(cartesians)
		out_mol = pybel.readstring('xyz',xyz_str)

		passes_rules = True
		if args.exp_rules:
			passes_rules = passes_custom_rules(out_mol,args)
		
		if args.check_geom and passes_rules:
			xyz_str = read_com_as_xyz()
			in_mol = pybel.readstring('xyz',xyz_str)
			passes_geom = passes_geometry_check(out_mol,in_mol,factor=args.length_criteria)
		

		# this part filters off conformers based on user-defined exp_rules
		passing_rules = True
		valid_mol_gen = True
		if len(args.exp_rules) >= 1:
			if TERMINATION == "normal" and IM_FREQS == 0:
				log.write("  ----- Exp_rules filter(s) will be applied to the output file -----\n")
				try:
					mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
					print_error_exp_rules=False
					if ob_compat and rdkit_compat:
						passing_rules = exp_rules_output(mol,args,log,file,print_error_exp_rules,ob_compat,rdkit_compat)
					os.remove(file.split('.')[0]+'.mol')
				except AttributeError:
					valid_mol_gen = False
					os.remove(file.split('.')[0]+'.mol')
					log.write("The file could not be converted into a mol object, exp_rules filter(s) will be disabled\n")

		passing_geom = True
		if args.check_geom and passing_rules and valid_mol_gen:
			if TERMINATION == "normal" and IM_FREQS == 0:
				log.write("  ----- Geometrical check will be applied to the output file -----\n")
				# this creates a mol object from the optimized log file
				mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
				# this creates a mol object from the input file
				try:
					os.chdir(w_dir_main +'/input_files/run_'+str(round_num))
					com_2_xyz_2_sdf(args,os.path.splitext(file)[0]+'.com')
					mol2,ob_compat,rdkit_compat = output_to_mol(file,'xyz',file)
					passing_geom = check_geom_filter(mol,mol2,args)
					# remove created files
					os.remove(file.split('.')[0]+'.xyz')
					os.remove(file.split('.')[0]+'.sdf')
					os.remove(file.split('.')[0]+'.mol')

				except FileNotFoundError:
					log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")

				os.chdir(w_dir)
				os.remove(file.split('.')[0]+'.mol')

		elif args.check_geom and not valid_mol_gen:
			log.write("The file could not be converted into a mol object, check_geom test will be disabled\n")

		# This part places the calculations in different folders depending on the type of termination
		finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr = create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,passing_rules,passing_geom,check_geom_qcorr)

		# check if gen or genecp are active
		# right now, QCORR is only working with Gaussian output files

		ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')

		# create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
		if IM_FREQS > 0 or TERMINATION != "normal" and not os.path.exists(w_dir_main+'/failed/run_'+str(round_num)+'/error/basis_set_error/'+file):
			create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT, orca_aux_section)

		# adding in the NMR componenet only to the finished files after reading from normally finished log files
		if TERMINATION == "normal" and IM_FREQS == 0 and passing_rules and passing_geom:
			if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole' or args.nics:

				#get coordinates
				ATOMTYPES, CARTESIANS,stand_or = get_coords_normal(outlines, stand_or, NATOMS, periodic_table, ATOMTYPES, CARTESIANS)

				if args.sp == 'gaussian':
					# creating new folder with new input gaussian files
					single_point_input_files = w_dir_fin+'/../G16-SP_input_files'

				elif args.sp == 'orca':
					# creating new folder with new input gaussian files
					single_point_input_files = w_dir_fin+'/../ORCA-SP_input_files'
				
				elif args.sp == 'turbomole':
					# creating new folder with new input gaussian files
					single_point_input_files = w_dir_fin+'/../TURBOMOLE-SP_input_files'

				if args.nics:
					nics_input_files = w_dir_fin+'/../G16-NICS_input_files'

				# Options for genecp
				if args.sp == 'gaussian':
					ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','gaussian')

				elif args.sp == 'orca':
					ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','orca')

				elif args.sp == 'turbomole':
					# Check for gen or genecp
					pass

				# this avoids problems related to genecp
				if genecp == None:
					basis_set_for_genecp = args.basis_set_sp
				elif genecp == 'genecp' or genecp == 'gen':
					basis_set_for_genecp = args.basis_set_genecp_atoms_sp

				# Sets the folder and find the log files to analyze
				for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,basis_set_for_genecp):
					if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole':
						if str(bs_sp).find('/') > -1:
							log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp.split('/')[0]))
						else:
							log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp))
					if args.nics:
						log.write('-> Creating NICS input files for {0} in {1}/{2}-{3}'.format(file,nics_input_files,lot_sp,bs_sp))

					# eliminates * from the name (since folders cannot contain * in their names)
					if bs_sp.find('**') > -1:
						bs_sp = bs_sp.replace('**','(d,p)')
					elif bs_sp.find('*') > -1:
						bs_sp = bs_sp.replace('*','(d)')

					if str(bs_sp).find('/') > -1:
						dir_name = str(lot_sp) + '-' + str(bs_sp.split('/')[0])
					else:
						dir_name = str(lot_sp) + '-' + str(bs_sp)

					keywords_opt = ''
					if args.sp == 'gaussian':
						if genecp == 'genecp' or  genecp == 'gen':
							keywords_opt = lot_sp + '/' + genecp
						else:
							keywords_opt = lot_sp + '/' + bs_sp
						if args.set_input_line_sp != 'None':
							keywords_opt += ' {0}'.format(args.set_input_line_sp)
						if args.empirical_dispersion_sp != 'None':
							keywords_opt += ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
						if args.solvent_model_sp != 'gas_phase':
							keywords_opt += ' scrf=({0},solvent={1})'.format(args.solvent_model_sp,args.solvent_name_sp)

					if args.charge_sp != 'None':
						CHARGE = args.charge_sp

					if args.mult_sp != 'None':
						MULT = args.mult_sp

					if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole':
						if not os.path.isdir(single_point_input_files+'/'+dir_name):
							os.makedirs(single_point_input_files+'/'+dir_name)
						os.chdir(single_point_input_files+'/'+dir_name)
						new_com_file('sp',w_dir_initial,log,single_point_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

					if args.nics:
						if not os.path.isdir(nics_input_files+'/'+dir_name):
							os.makedirs(nics_input_files+'/'+dir_name)
						os.chdir(nics_input_files+'/'+dir_name)
						new_com_file('nics',w_dir_initial,log,nics_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

	#write to csv ana_data
	if duplicates=='None':
		ana_data.at[0,'Total files'] = len(log_files)
	else:
		ana_data.at[0,'Total files'] = len(log_files)+int(duplicates) # since duplicates are moved before anything else
	ana_data.at[0,'Normal termination'] = finished
	ana_data.at[0,'Imaginary frequencies'] = imag_freq
	ana_data.at[0,'SCF error'] = scf_error
	ana_data.at[0,'Basis set error'] =  atom_error
	ana_data.at[0,'Other errors'] = other_error
	ana_data.at[0,'Unfinished'] = unfinished
	if args.dup:
		ana_data.at[0,'Duplicates'] = duplicates
	if len(args.exp_rules) >= 1:
		ana_data.at[0,'Exp_rules filter'] = exp_rules_qcorr
	if args.check_geom:
		ana_data.at[0,'Geometry changed'] = check_geom_qcorr

	if not os.path.isdir(w_dir_main+'/csv_files/'):
		os.makedirs(w_dir_main+'/csv_files/')
	ana_data.to_csv(w_dir_main+'/csv_files/Analysis-Data-QCORR-run_'+str(round_num)+'.csv',index=False)

def analyze_single_output(file,log,w_dir,args):
	################# Output Parsing ##############
	# read the file
	log.write(file)
	filepath = Path(f'{w_dir}/{file}')
	GOF = GaussianOutputFile(filepath)
	
	if GOF.im_freqs: 
		cartesians = GOF.fix_imaginary_frequencies()
	else:
		cartesians = GOF.cartesians

	############################################################################
	######### Handle the logic of each output file and what to do with them
	
	# this part filters off conformers based on user-defined exp_rules
	passing_rules = True
	valid_mol_gen = True
	xyz_str = GOF.to_xyz(cartesians)
	mol = pybel.readstring('xyz',xyz_str)
	if len(args.exp_rules) >= 1:
		if TERMINATION == "normal" and IM_FREQS == 0:
			log.write("  ----- Exp_rules filter(s) will be applied to the output file -----\n")
			try:
				mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
				print_error_exp_rules=False
				if ob_compat and rdkit_compat:
					passing_rules = exp_rules_output(mol,args,log,file,print_error_exp_rules,ob_compat,rdkit_compat)
				os.remove(file.split('.')[0]+'.mol')
			except AttributeError:
				valid_mol_gen = False
				os.remove(file.split('.')[0]+'.mol')
				log.write("The file could not be converted into a mol object, exp_rules filter(s) will be disabled\n")

	passing_geom = True
	if args.check_geom and passing_rules and valid_mol_gen:
		if TERMINATION == "normal" and IM_FREQS == 0:
			log.write("  ----- Geometrical check will be applied to the output file -----\n")
			# this creates a mol object from the optimized log file
			mol,ob_compat,rdkit_compat = output_to_mol(file,'log',file)
			# this creates a mol object from the input file
			try:
				os.chdir(w_dir_main +'/input_files/run_'+str(round_num))
				com_2_xyz_2_sdf(args,os.path.splitext(file)[0]+'.com')
				mol2,ob_compat,rdkit_compat = output_to_mol(file,'xyz',file)
				passing_geom = check_geom_filter(mol,mol2,args)
				# remove created files
				os.remove(file.split('.')[0]+'.xyz')
				os.remove(file.split('.')[0]+'.sdf')
				os.remove(file.split('.')[0]+'.mol')

			except FileNotFoundError:
				log.write("x  No com file were found for "+file+", the check_geom test will be disabled for this calculation")

			os.chdir(w_dir)
			os.remove(file.split('.')[0]+'.mol')

	elif args.check_geom and not valid_mol_gen:
		log.write("The file could not be converted into a mol object, check_geom test will be disabled\n")

	# This part places the calculations in different folders depending on the type of termination
	finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,check_geom_qcorr = create_folder_move_log_files(w_dir,w_dir_main,round_num,file,IM_FREQS,TERMINATION,ERRORTYPE,w_dir_fin,finished,unfinished,atom_error,scf_error,imag_freq,other_error,exp_rules_qcorr,passing_rules,passing_geom,check_geom_qcorr)

	# check if gen or genecp are active
	# right now, QCORR is only working with Gaussian output files

	ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'analysis','gaussian')

	# create folders and set level of theory in COM files to fix imaginary freqs or not normal terminations
	if IM_FREQS > 0 or TERMINATION != "normal" and not os.path.exists(w_dir_main+'/failed/run_'+str(round_num)+'/error/basis_set_error/'+file):
		create_folder_and_com(w_dir,w_dir_main,round_num,log,NATOMS,ATOMTYPES,CARTESIANS,args,TERMINATION,IM_FREQS,w_dir_fin,file,lot,bs,bs_gcp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,ERRORTYPE,input_route,w_dir_initial,name,CHARGE, MULT, orca_aux_section)

	# adding in the NMR componenet only to the finished files after reading from normally finished log files
	if TERMINATION == "normal" and IM_FREQS == 0 and passing_rules and passing_geom:
		if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole' or args.nics:

			#get coordinates
			ATOMTYPES, CARTESIANS,stand_or = get_coords_normal(outlines, stand_or, NATOMS, periodic_table, ATOMTYPES, CARTESIANS)

			if args.sp == 'gaussian':
				# creating new folder with new input gaussian files
				single_point_input_files = w_dir_fin+'/../G16-SP_input_files'

			elif args.sp == 'orca':
				# creating new folder with new input gaussian files
				single_point_input_files = w_dir_fin+'/../ORCA-SP_input_files'
				
			elif args.sp == 'turbomole':
				# creating new folder with new input gaussian files
				single_point_input_files = w_dir_fin+'/../TURBOMOLE-SP_input_files'

			if args.nics:
				nics_input_files = w_dir_fin+'/../G16-NICS_input_files'

			# Options for genecp
			if args.sp == 'gaussian':
				ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','gaussian')

			elif args.sp == 'orca':
				ecp_list,ecp_genecp_atoms,ecp_gen_atoms,genecp,orca_aux_section = check_for_gen_or_genecp(ATOMTYPES,args,'sp','orca')

			elif args.sp == 'turbomole':
				# Check for gen or genecp
				pass

			# this avoids problems related to genecp
			if genecp == None:
				basis_set_for_genecp = args.basis_set_sp
			elif genecp == 'genecp' or genecp == 'gen':
				basis_set_for_genecp = args.basis_set_genecp_atoms_sp

			# Sets the folder and find the log files to analyze
			for lot_sp,bs_sp,bs_gcp_sp in zip(args.level_of_theory_sp,args.basis_set_sp,basis_set_for_genecp):
				if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole':
					if str(bs_sp).find('/') > -1:
						log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp.split('/')[0]))
					else:
						log.write('-> Creating new single point files for {0} in {1}/{2}-{3}'.format(file,single_point_input_files,lot_sp,bs_sp))
				if args.nics:
					log.write('-> Creating NICS input files for {0} in {1}/{2}-{3}'.format(file,nics_input_files,lot_sp,bs_sp))

				# eliminates * from the name (since folders cannot contain * in their names)
				if bs_sp.find('**') > -1:
					bs_sp = bs_sp.replace('**','(d,p)')
				elif bs_sp.find('*') > -1:
					bs_sp = bs_sp.replace('*','(d)')

				if str(bs_sp).find('/') > -1:
					dir_name = str(lot_sp) + '-' + str(bs_sp.split('/')[0])
				else:
					dir_name = str(lot_sp) + '-' + str(bs_sp)

				keywords_opt = ''
				if args.sp == 'gaussian':
					if genecp == 'genecp' or  genecp == 'gen':
						keywords_opt = lot_sp + '/' + genecp
					else:
						keywords_opt = lot_sp + '/' + bs_sp
					if args.set_input_line_sp != 'None':
						keywords_opt += ' {0}'.format(args.set_input_line_sp)
					if args.empirical_dispersion_sp != 'None':
						keywords_opt += ' empiricaldispersion={0}'.format(args.empirical_dispersion_sp)
					if args.solvent_model_sp != 'gas_phase':
						keywords_opt += ' scrf=({0},solvent={1})'.format(args.solvent_model_sp,args.solvent_name_sp)

				if args.charge_sp != 'None':
					CHARGE = args.charge_sp

				if args.mult_sp != 'None':
					MULT = args.mult_sp

				if args.sp == 'gaussian' or args.sp == 'orca' or args.sp == 'turbomole':
					if not os.path.isdir(single_point_input_files+'/'+dir_name):
						os.makedirs(single_point_input_files+'/'+dir_name)
					os.chdir(single_point_input_files+'/'+dir_name)
					new_com_file('sp',w_dir_initial,log,single_point_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

				if args.nics:
					if not os.path.isdir(nics_input_files+'/'+dir_name):
						os.makedirs(nics_input_files+'/'+dir_name)
					os.chdir(nics_input_files+'/'+dir_name)
					new_com_file('nics',w_dir_initial,log,nics_input_files+'/'+dir_name,file,args,keywords_opt,name,CHARGE,MULT,NATOMS,ATOMTYPES,CARTESIANS,genecp,ecp_list,ecp_genecp_atoms,ecp_gen_atoms,TERMINATION,IM_FREQS,bs_sp,lot_sp,bs_gcp_sp,orca_aux_section)

def classify_files(inputs,outputs,output_objects):
	finished = []
	im_freq = []
	unfinished = []
	atombasis_error = []
	scf_error = []
	unknown_error = []
	custom_rules = []
	geom_check = []
	for ifile, ofile, ofile_o in zip(inputs,outputs,output_objects):
		has_im_freq = bool(ofile_o.im_freq)
		normal_term = ofile_o.termination == 'normal'
		rule_pass = passes_custom_rules() 
		geom_pass = passes_geometry_check()
		good_end = normal_term and (not has_im_freq) and rule_pass and geom_pass
		if good_end: 
			finished.append((ifile, ofile, ofile_o))
		elif has_im_freq: 
			im_freq.append((ifile, ofile, ofile_o))
		elif not normal_term:
			if ofile_o.termination == 'unfinished': 
				unfinished.append((ifile, ofile, ofile_o))
			elif ofile_o.errortype == 'atomicbasiserror': 
				atombasis_error.append((ifile, ofile, ofile_o))
			elif ofile_o.errortype == 'SCFerror':
				scf_error.append((ifile, ofile, ofile_o))
			else:
				unknown_error.append((ifile, ofile, ofile_o))
		elif not rule_pass:
			custom_rules.append((ifile, ofile, ofile_o))
		elif not geom_pass:
			geom_check.append((ifile, ofile, ofile_o))
	classification = dict()
	classification['finished'] = finished
	classification['imag_freq'] = im_freq
	classification['atombasis_error'] = atombasis_error
	classification['scf_error'] = scf_error
	classification['unknown_error'] = unknown_error
	classification['unfinished'] = unfinished
	classification['custom_rules'] = custom_rules
	classification['geom_check'] = geom_check
	return classification
def handle_unfinished_calculations(ofile,ofile_object,destination):
	pass



# Main function
def qcorr_gaussian_main(duplicates,w_dir_initial,args,log):

	# when you run analysis in a folder full of output files
	qmcalc_folder = Path(f'{w_dir_initial}/QMCALC')
	if not qmcalc_folder.exists():
		w_dir_main = os.getcwd()
		w_dir_fin = w_dir_main+'/success/output_files'
		for lot,bs,bs_gcp in zip(args.level_of_theory, args.basis_set,args.basis_set_genecp_atoms):
			if not os.path.isdir(w_dir_main+'/dat_files/'):
				os.makedirs(w_dir_main+'/dat_files/')
			w_dir,n_run = check_for_final_folder(w_dir_main)
			os.chdir(w_dir)
			log = Logger(f'{w_dir_main}/dat_files/pyCONFORT-QCORR-run_{n_run}', args.output_name)
			ana_data = creation_of_ana_csv(args)
			log.write("\no  Analyzing output files in {}\n".format(w_dir))
			log_files = get_com_or_log_out_files('output',None)
			if len(log_files) == 0:
				log.write('x  There are no output files in this folder.')
			com_files = get_com_or_log_out_files('input',None)
			output_analyzer(duplicates,log_files, com_files, w_dir, w_dir_main, 
							lot, bs, bs_gcp, args, w_dir_fin, w_dir_initial, 
							log, ana_data, n_run)
			os.chdir(w_dir_main)
	# when you specify multiple levels of theory
	else:
		folders = [item for item in (qmcalc_folder/'G16').iterdir() if item.isdir()]
		dat_folders = []
		for folder in folders: 
			dat_folder = folder/'dat_files'
			dat_folder.mkdir(exist_ok=True)
			last_run, n_run = check_for_final_folder(str(folder))
			dat_file = dat_folder/f'pyCONFORT-QCORR-run_{n_run:02d}'
			log = Logger(str(dat_file),args.output_name)
			outputs_folder = folder/'success/output_files'
			ana_data = creation_of_ana_csv(args)
			log.write(f"\no  Analyzing output files in {last_run}\n")
			is_output = lambda x: x.is_file() and x.suffix in OUTPUT_SUFFIXES
			is_input = lambda x: x.is_file() and x.suffix in INPUT_SUFFIXES
			outputs = [file for file in folder.iterdir() if is_output(file)]
			inputs = [file for file in folder.iterdir() if is_input(file)]
			output_objects = [GaussianOutputFile(file) for file in outputs]
			classification = classify_files(inputs,outputs,output_objects)
			
			# Move files to their proper places
			# create folders
			inputs_folder = folder/f'input_files/run_{n_run}'
			imag_freq_folder = folder/f'failed/run_{n_run}/imag_freq' 
			errors_folder = folder/f'failed/run_{n_run}/error'
			scf_error_folder = errors_folder/f'scf_error'
			basisset_error_folder = errors_folder/f'basis_set_error'
			unknown_error_folder = errors_folder/f'unknown_error'
			unfinished_folder = folder/f'failed/run_{n_run}/unfinished'
			rules_folder = folder/f'failed/run_{n_run}/exp_rules_filter'
			geomcheck_folder = folder/f'failed/run_{n_run}/geometry_changed'

			# Prepare next run for imag_freqs files 
			for ifile, ofile, ofile_o in classification['imag_freq']: 
				# move the inputs
				# move the outputs
				ofile_o.cartesians = ofile_o.fix_imaginary_frequencies()
				# write the new inputs

			### Repeat for unknown_error  
			#classification['unknown_error']
			### Repeat for scf_error  
			#classification['scf_error']
			### Repeat for unfinished  
			#classification['unfinished']

			# Good endings 
			for ifile,ofile,ofile_o in classification['finished']:
				# move the inputs
				# move the outputs
				pass # Check if SP or NICS... and use appropiate function/template etc etc
