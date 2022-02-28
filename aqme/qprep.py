######################################################.
#        This file stores the QPREP class            #
######################################################.

import os
import sys
import glob
import pandas as pd
from rdkit import Chem
from aqme.utils import (
	move_file,
	load_variables)
from aqme.qprep_utils import qprep_coords
from openbabel import pybel
from pathlib import Path


class qprep:
	"""
	Class containing all the functions from the QPREP module related to Gaussian input files
	"""

	def __init__(self,**kwargs):

		# load default and user-specified variables
		self.args = load_variables(kwargs,'qprep')

		if self.args.destination is None:
			self.args.destination = self.args.w_dir_main.joinpath("QCALC")
		else:
			self.args.destination = Path(self.args.destination)

		if self.args.qm_input == "":
			print("x  No keywords line was specified! (i.e. qm_input=KEYWORDS_LINE).")
			sys.exit()

		# w_dir_initial is included to avoid problems in jupyter notebooks
		w_dir_initial = Path(os.getcwd())

		# write input files	
		os.chdir(self.args.w_dir_main)
		for file in self.args.files:
			found_coords = True
			if self.args.atom_types == [] or self.args.cartesians == []:
				found_coords = False
				self.atom_types,self.cartesians,found_coords,self = qprep_coords(self.args.w_dir_main,file,found_coords,self)
			else:
				self.atom_types = self.args.atom_types
				self.cartesians = self.args.cartesians
			if not found_coords:
				continue
			self.args.n_atoms = len(self.atom_types)
			self.args.molecule = file.split('.')[0]
			comfile = self.write()
			move_file(self.args.destination, self.args.w_dir_main, comfile)
		os.chdir(w_dir_initial)


	def get_header(self):
		'''
		Gets the part of the input file above the molecular coordinates.
		'''

		txt = ''

		if self.args.program.lower() == 'gaussian':
			if self.args.chk:
				txt += f'%chk={self.args.molecule}.chk\n'
			txt += f'%nprocshared={self.args.nprocs}\n'
			txt += f'%mem={self.args.mem}\n'
			txt += f'# {self.args.qm_input}'
			txt += '\n\n'
			txt += f'{self.args.molecule}\n\n'
			txt += f'{self.charge} {self.mult}\n'

		elif self.args.program.lower() == 'orca':
			txt += f'# {self.args.molecule}\n'
			if self.args.mem.find('GB'):
				mem_orca = int(self.args.mem.split('GB')[0])*1000
			elif self.args.mem.find('MB'):
				mem_orca = self.args.mem.split('MB')[0]
			elif self.args.args.mem.find('MW'):
				mem_orca = self.args.mem.split('MW')[0]
			txt += f'%maxcore {mem_orca}\n'
			txt += f'%pal nprocs {self.args.nprocs} end\n'
			txt += f'! {self.args.qm_input}\n'
			txt += f'* xyz {self.charge} {self.mult}\n'

		return txt

	def get_tail(self):
		"""
		Gets the part of the input file below the molecular coordinates.
		"""

		txt = ""

		if self.args.program.lower() == "gaussian":
			if self.args.gen_atoms is not None and len(self.args.gen_atoms) > 0:
				# writes part for Gen/GenECP
				ecp_used, ecp_not_used, gen_type = [], [], "gen"
				if self.args.qm_input.lower().find("genecp") > -1:
					gen_type = "genecp"

				for _, element_ecp in enumerate(self.atom_types):
					if element_ecp in self.args.gen_atoms and element_ecp not in ecp_used:
						ecp_used.append(element_ecp)
					elif (
						element_ecp not in self.args.gen_atoms
						and element_ecp not in ecp_not_used
					):
						ecp_not_used.append(element_ecp)

				if len(ecp_not_used) > 0:
					elements_not_used = " ".join([f"{sym}" for sym in ecp_not_used])
					txt += f"{elements_not_used} 0\n{self.args.bs}\n****\n"
				if len(ecp_used) > 0:
					elements_used = " ".join([f"{sym}" for sym in ecp_used])
					txt += f"{elements_used} 0\n{self.args.bs_gen}\n****\n"

				if gen_type == "genecp" and len(ecp_used) > 0:
					txt += "\n"
					txt += f"{elements_used} 0\n{self.args.bs_gen}\n****\n"

				txt += "\n"

			# writes final section if selected
			if self.args.qm_end != "":
				txt += f"{self.args.qm_end}\n\n"

		return txt

	def write(self):

		if self.args.program.lower() == 'gaussian':
			extension = 'com'
		elif self.args.program.lower() == 'orca':
			extension = 'inp'
		if self.args.suffix != '':
			comfile = f'{self.args.molecule}_{self.args.suffix}.{extension}'
		else:
			comfile = f'{self.args.molecule}.{extension}'

		if os.path.exists(comfile):
			os.remove(comfile)

		header = self.get_header()
		tail = self.get_tail()

		fileout = open(comfile, "w")
		fileout.write(header)

		for atom in range(0, self.args.n_atoms):
			fileout.write(
				"{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}".format(
					self.atom_types[atom],
					self.cartesians[atom][0],
					self.cartesians[atom][1],
					self.cartesians[atom][2],
				)
			)
			if atom != self.args.n_atoms - 1:
				fileout.write("\n")

		if self.args.program.lower() == "gaussian":
			fileout.write("\n\n")
		elif self.args.program.lower() == "orca":
			fileout.write("\n*")

		fileout.write(tail)
		fileout.close()

		return comfile


# Aux Functions for QM input generation
def get_molecule_list(filepath, lowest_only=False, lowest_n=False, energy_threshold=0.0):
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
	# read in dup_data to get the overall charge of molecules
	invalid_files = []
	try:
		charge_data = pd.read_csv(filepath, usecols=["Molecule", "Overall charge","Mult"])
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
				mult = mol.data["Mult"]
			charge_data.at[i, "Molecule"] = name
			charge_data.at[i, "Overall charge"] = charge
			charge_data.at[i, "Mult"] = mult
	return charge_data, invalid_files