import shutil
import os
import re
from pathlib import Path 

import rdkit
from rdkit import Chem

try:
	import pybel
except ImportError:
	from openbabel import pybel

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
				syms = [possible_atoms[atnum] for atnum in atnums]
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

#  Functions brought from qprep

# WRITE SDF FILES FOR xTB AND ANI1
def write_confs(conformers, energies,selectedcids, name, args, program,log):
	if len(conformers) > 0:
		# name = name.split('_'+args.CSEARCH)[0]# a bit hacky
		sdwriter = Chem.SDWriter(name+'_'+program+args.output)

		write_confs = 0
		for cid in selectedcids:
			sdwriter.write(conformers[cid])
			write_confs += 1

		if args.verbose:
			log.write("o  Writing "+str(write_confs)+ " conformers to file " + name+'_'+program+args.output)
		sdwriter.close()
	else:
		log.write("x  No conformers found!")

# MOVES SDF FILES TO THEIR CORRESPONDING FOLDERS
def moving_files(destination,src,file):
	try:
		os.makedirs(destination)
		shutil.move(os.path.join(src, file), os.path.join(destination, file))
	except OSError:
		if  os.path.isdir(destination):
			shutil.move(os.path.join(src, file), os.path.join(destination, file))
		else:
			raise


items= """X
 H                                                                                                  He
Li Be                                                                            B   C   N   O   F  Ne
Na Mg                                                                           Al  Si   P   S  Cl  Ar
 K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu  Zn  Ga  Ge  As  Se  Br  Kr
Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag  Cd  In  Sn  Sb  Te   I  Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au  Hg  Tl  Pb  Bi  Po  At  Rn
Fr Ra Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo
""" # The "X" is necessary to ensure that index == AtNum
periodic_table = items.replace('\n',' ').strip().split()