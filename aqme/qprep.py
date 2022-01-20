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
from rdkit import Chem
from aqme.utils import load_from_yaml, Logger, move_file
from aqme.argument_parser import set_options
from openbabel import pybel
from pathlib import Path


# template classes
class qprep:

	"""
	Class containing all the functions from the QPREP module related to Gaussian input files.

	Parameters
	----------
	mol : mol object with 3D coordinates
		Optionally, starts from a mol object to prepare the input QM file
	w_dir_main : str
		Working directory
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
	varfile : str
		Option to parse the variables using a yaml file (specify the filename)
	kwargs : argument class
		Specify any arguments from the QCORR module
	"""

	def __init__(
		self,
		mol=None,
		w_dir_main=os.getcwd(),
		destination=None,
		molecule="",
		charge=0,
		mult=1,
		chk=False,
		mem="8GB",
		nprocs=4,
		suffix='',
		atom_types=[],
		cartesians=[],
		bs_gen="",
		bs="",
		gen_atoms=None,
		program="gaussian",
		qm_input="",
		qm_end="",
		varfile=None,
		**kwargs,
	):

		self.mol = mol
		self.w_dir_main = Path(w_dir_main)
		self.molecule = molecule
		self.chk = chk
		self.mem = mem
		self.nprocs = nprocs
		self.suffix = suffix
		self.charge = charge
		self.mult = mult
		self.qm_input = qm_input
		self.bs_gen = bs_gen
		self.bs = bs
		self.gen_atoms = gen_atoms
		self.qm_end = qm_end
		self.program = program

		if destination is None:
			self.destination = Path(os.getcwd()+"/QCALC/")
      
		else:
			self.destination = Path(destination)
     
    dat_folder = self.destination.joinpath("dat_files/")

		if "options" in kwargs:
			self.args = kwargs["options"]
		else:
			self.args = set_options(kwargs)

		self.args.varfile = varfile

		if varfile is not None:
			self.args, self.log = load_from_yaml(self.args, self.log)
			self.w_dir_main = Path(self.args.w_dir_main)
			self.destination = self.args.destination
			if destination is None:
				self.destination = self.w_dir_main.joinpath("QCALC")
			else:
				self.destination = Path(self.args.destination)
			self.charge = self.args.charge
			self.mult = self.args.mult
			if self.charge == None:
				self.charge = 0
			if self.mult == None:
				self.mult = 1
			self.chk = self.args.chk
			self.mem = self.args.mem
			self.nprocs = self.args.nprocs
			self.qm_input = self.args.qm_input
			self.gen_atoms = self.args.gen_atoms
			self.bs_gen = self.args.bs_gen
			self.bs = self.args.bs
			self.qm_end = self.args.qm_end
			self.program = self.args.program
			self.suffix = self.args.suffix

		self.atom_types = atom_types
		self.cartesians = cartesians

		if self.mol != None:
			self.atom_types = [ atom.GetSymbol() for _, atom in enumerate(self.mol.GetAtoms())]
			self.cartesians = self.mol.GetConformers()[0].GetPositions()
		self.n_atoms = len(self.atom_types)

		if self.qm_input == "":
			sys.exit("x  No keywords line was specified! (qm_input=KEYWORDS_LINE).")

		if molecule == "":
			sys.exit("x  No name was specified! (molecule=NAME).")
		
		comfile = self.write()
		move_file(self.destination, self.w_dir_main, comfile)


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
		"""
		Gets the part of the input file below the molecular coordinates.
		"""

		txt = ""

		if self.program.lower() == "gaussian":
			if self.gen_atoms is not None and len(self.gen_atoms) > 0:
				# writes part for Gen/GenECP
				ecp_used, ecp_not_used, gen_type = [], [], "gen"
				if self.qm_input.lower().find("genecp") > -1:
					gen_type = "genecp"

				for _, element_ecp in enumerate(self.atom_types):
					if element_ecp in self.gen_atoms and element_ecp not in ecp_used:
						ecp_used.append(element_ecp)
					elif (
						element_ecp not in self.gen_atoms
						and element_ecp not in ecp_not_used
					):
						ecp_not_used.append(element_ecp)

				if len(ecp_not_used) > 0:
					elements_not_used = " ".join([f"{sym}" for sym in ecp_not_used])
					txt += f"{elements_not_used} 0\n{self.bs}\n****\n"
				if len(ecp_used) > 0:
					elements_used = " ".join([f"{sym}" for sym in ecp_used])
					txt += f"{elements_used} 0\n{self.bs_gen}\n****\n"

				if gen_type == "genecp" and len(ecp_used) > 0:
					txt += "\n"
					txt += f"{elements_used} 0\n{self.bs_gen}\n****\n"

				txt += "\n"

			# writes final section if selected
			if self.qm_end != "":
				txt += f"{self.qm_end}\n\n"

		return txt

	def write(self):

		if self.program.lower() == 'gaussian':
			extension = 'com'
		elif self.program.lower() == 'orca':
			extension = 'inp'
		if self.suffix != '':
			comfile = f'{self.molecule}_{self.suffix}.{extension}'
		else:
			comfile = f'{self.molecule}.{extension}'

		if os.path.exists(comfile):
			os.remove(comfile)

		header = self.get_header()
		tail = self.get_tail()

		fileout = open(comfile, "w")
		fileout.write(header)

		for atom in range(0, self.n_atoms):
			fileout.write(
				"{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}".format(
					self.atom_types[atom],
					self.cartesians[atom][0],
					self.cartesians[atom][1],
					self.cartesians[atom][2],
				)
			)
			if atom != self.n_atoms - 1:
				fileout.write("\n")

		if self.program.lower() == "gaussian":
			fileout.write("\n\n")
		elif self.program.lower() == "orca":
			fileout.write("\n*")

		fileout.write(tail)
		fileout.close()

		return comfile


# Aux Functions for QM input generation
def get_molecule_list(
	filepath, lowest_only=False, lowest_n=False, energy_threshold=0.0
):
	out_molecules = []

	molecules = [mol for mol in Chem.SDMolSupplier(filepath, removeHs=False)]
	energies = [float(mol.GetProp("Energy")) for mol in molecules]
	min_energy = energies[0]
	for mol, energy in zip(molecules, energies):
		is_in_threshold = energy - min_energy < energy_threshold
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


def load_charge_data(filepath, backup_files):
	# read in dup_data to get the overall charge of MOLECULES
	invalid_files = []
	try:
		charge_data = pd.read_csv(filepath, usecols=["Molecule", "Overall charge"])
	except:
		charge_data = pd.DataFrame()
		for i, sdf_file in enumerate(backup_files):
			if not (Path(sdf_file).exists()):
				invalid_files.append(sdf_file)
				maxsplit = 1
				if "filter" in sdf_file:
					maxsplit += 1
				name = sdf_file.rsplit("_", maxsplit[0])
				charge = "Invalid"
			else:
				mol = next(pybel.readfile(sdf_file))
				name = mol.title.split(maxsplit=1)[0]
				charge = mol.data["Real charge"]
			charge_data.at[i, "Molecule"] = name
			charge_data.at[i, "Overall charge"] = charge
	return charge_data, invalid_files
