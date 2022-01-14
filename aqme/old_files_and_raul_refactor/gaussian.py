import numpy 
from aqme.argument_parser import periodic_table

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
		termination,errortype = self._get_qcorr_params()
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
	def _get_qcorr_params(self): 
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
	def to_xyz(self,cartesians=None):
		if cartesians is None: 
			cartesians = self.cartesians
		n_atoms = len(cartesians)
		xyz = []
		for sym,(x,y,z) in zip(self.atomtypes,cartesians): 
			xyz.append(f'{sym}    {x: 0.6f}    {y: 0.6f}    {z: 0.6f}')
		return f"{n_atoms}\n{self.name}\n{'\n'.join(xyz)}\n"

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