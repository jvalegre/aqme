"""
Parameters
----------
   files : mol object, str or list of str, default=None
      This module prepares input QM file(s). Formats accepted: mol object(s), 
      Gaussian or ORCA LOG/OUT output files, JSON, XYZ, SDF, PDB. Also, 
      lists can be used (i.e. [FILE1.log, FILE2.log] or \*.FORMAT such as \*.json).
   atom_types : list of str, default=[]
      (If files is None) List containing the atoms of the system
   cartesians : list of str, default=[]
      (If files is None) Cartesian coordinates used for further processing
   w_dir_main : str, default=os.getcwd()
      Working directory
   destination : str, default=None,
      Directory to create the input file(s)
   varfile : str, default=None
      Option to parse the variables using a yaml file (specify the filename)
   program : str, default=None
      Program required to create the new input files. Current options: 'gaussian', 'orca'
   qm_input : str, default=''
      Keywords line for new input files (i.e. 'B3LYP/6-31G opt freq')
   qm_end : str, default=''
      Final line(s) in the new input files
   charge : int, default=None
      Charge of the calculations used in the following input files. If charge isn't defined, it defaults to 0
   mult : int, default=None
      Multiplicity of the calculations used in the following input files. If mult isn't defined, it defaults to 1
   suffix : str, default=''
      Suffix for the new input files (i.e. FILENAME_SUFFIX.com for FILENAME.log)
   prefix : str, default=''  
      Prefix added to all the names  
   chk : bool, default=False
      Include the chk input line in new input files for Gaussian calculations
   oldchk : bool, default=False
      Include the oldchk input line in new input files for Gaussian calculations
   chk_path : str, default=''
      PATH to store CHK files. For example, if chk_path='root/user/FILENAME.chk, the chk line of the input file would be
      %chk=root/user/FILENAME.chk
   oldchk_path : str, default=''
      PATH to read CHK files with %oldchk. For example, if oldchk_path='root/user/FILENAME.chk, the oldchk line of the input file would be
      %oldchk=root/user/FILENAME.chk
   mem : str, default='4GB'
      Memory for the QM calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor
   nprocs : int, default=None
      Number of processors used in the QM calculations
   gen_atoms : list of str, default=[]
      Atoms included in the gen(ECP) basis set (i.e. ['I','Pd'])
   bs_gen : str, default=''
      Basis set used for gen(ECP) atoms (i.e. 'def2svp')
   bs_nogen : str, default=''
      Basis set used for non gen(ECP) atoms in gen(ECP) calculations (i.e. '6-31G*')
   lowest_only : bool, default=False
      Only create input for the conformer with lowest energy of the SDF file
   lowest_n : int, default=None
      Only create inputs for the n conformers with lowest energy of the SDF file
   e_threshold_qprep : float, default=None
      Only create inputs for conformers below the energy threshold (to the lowest conformer)
      of the SDF file
"""
######################################################.
#        This file stores the QPREP class            #
######################################################.

import os
import subprocess
import sys
import glob
import time
import json
import pandas as pd
from pkg_resources import resource_filename

from aqme.utils import (
    cclib_atoms_coords,
    QM_coords,
    read_file,
    move_file,
    load_variables,
    read_xyz_charge_mult,
    mol_from_sdf_or_mol_or_mol2,
    add_prefix_suffix,
    check_files,
    check_dependencies,
    set_destination
)

from aqme.csearch.crest import xyzall_2_xyz
from pathlib import Path
from rdkit import Chem

TEMPLATES_PATH = Path(resource_filename("aqme", "templates"))

class qprep:
    """
    Class containing all the functions from the QPREP module related to Gaussian input files
    """

    def __init__(self, create_dat=True, **kwargs):

        start_time_overall = time.time()

        # load default and user-specified variables
        self.args = load_variables(kwargs, "qprep", create_dat=create_dat)

        # check whether dependencies are installed
        _ = check_dependencies(self)

        # retrieves the different files to run in QPREP
        _ = check_files(self,'qprep')

        file_format = os.path.basename(Path(self.args.files[0])).split('.')[-1]
        if file_format.lower() not in ['sdf', 'xyz', 'pdb', 'log', 'out', 'json']:
            self.args.log.write(f"\nx  The format used ({file_format}) is not compatible with QPREP! Formats accepted: sdf, xyz, pdb, log, out, json")
            self.args.log.finalize()
            sys.exit()

        qprep_program = True
        if self.args.program is None:
            qprep_program = False
        if qprep_program:
            if self.args.program.lower() not in ["gaussian", "orca"]:
                qprep_program = False
        if not qprep_program:
            self.args.log.write('\nx  Program not supported for QPREP input file creation! Specify: program="gaussian" (or "orca")')
            self.args.log.finalize()
            sys.exit()

        # set number of processors
        if self.args.nprocs is None:
            self.args.nprocs = 8

        destination = set_destination(self,'QCALC')

        # check if qm_input is not empty
        if self.args.qm_input == "" and create_dat:
            self.args.log.write("x  No keywords line was specified! (i.e. qm_input=KEYWORDS_LINE).")
            self.args.log.finalize()
            sys.exit()

        # check if functionals and basis sets used are correct
        # so far, it only works for Gaussian
        if self.args.program.lower() == 'gaussian' and create_dat:
            _ = self.check_level_of_theory()

        # checks for gen/genecp
        if self.args.gen_atoms != [] and self.args.bs_nogen == "" and create_dat:
            self.args.log.write("x  Atoms for Gen(ECP) were specified (gen_atoms=[ATOM_LIST]) but no basis set was included for non-Gen(ECP) atoms (i.e. bs_nogen=BASIS_SET).")
            self.args.log.finalize()
            sys.exit()

        elif self.args.gen_atoms != [] and self.args.bs_gen == "" and create_dat:
            self.args.log.write("x  Atoms for Gen(ECP) were specified (gen_atoms=[ATOM_LIST]) but no basis set was included for Gen(ECP) atoms (i.e. bs_gen=BASIS_SET).")
            self.args.log.finalize()
            sys.exit()

        # write input files
        for file in self.args.files:
            name = os.path.basename(Path(file)).split(".")[0]
            if file_format.lower() in ["sdf", "xyz", "pdb"]:
                sdf_files = []
                if file_format.lower() == "xyz":
                    # separate the parent XYZ file into individual XYZ files
                    xyzall_2_xyz(file, f"{self.args.w_dir_main}/{name}")
                    for conf_file in glob.glob(
                        f"{self.args.w_dir_main}/{name}_conf_*.xyz"
                    ):
                        if self.args.charge is None:
                            charge_xyz, _ = read_xyz_charge_mult(conf_file)
                        else:
                            charge_xyz = self.args.charge
                        if self.args.mult is None:
                            _, mult_xyz = read_xyz_charge_mult(conf_file)
                        else:
                            mult_xyz = self.args.mult
                        # generate SDF files from XYZ with Openbabel
                        command_xyz = [
                            "obabel",
                            "-ixyz",
                            conf_file,
                            "-osdf",
                            f"-O{conf_file.split('.xyz')[0]}.sdf",
                            "--property",
                            f"Real charge={str(charge_xyz)}",
                            ";",
                            f"Mult={str(mult_xyz)}",
                        ]
                        subprocess.run(
                            command_xyz,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                        )
                        sdf_files.append(f"{conf_file.split('.xyz')[0]}.sdf")
                        # delete individual XYZ files
                        os.remove(conf_file)

                elif file_format.lower() == "pdb":
                    command_pdb = [
                        "obabel",
                        "-ipdb",
                        f"{file}",
                        "-osdf",
                        f"-O{file.split('.pdb')[0]}.sdf",
                    ]
                    subprocess.run(
                        command_pdb,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    sdf_files.append(f"{file.split('.pdb')[0]}.sdf")

                else:

                    sdf_files.append(file)

                for sdf_file in sdf_files:
                    try:
                        self.sdf_2_com(sdf_file, destination, file_format)

                    except OSError:
                        self.args.log.write(f"x  {name} couldn't be processed!")
                        continue

                    if create_dat:
                        self.args.log.write(f"o  {name} successfully processed at {destination}")

                    if file_format.lower() in ["xyz", "pdb"]:
                        # delete SDF files when the input was an XYZ/PDB file
                        os.remove(sdf_file)

            # for Gaussian output files (LOG/OUT), JSON files and MOL objects
            else:
                atom_types, cartesians, charge, mult, found_coords = self.qprep_coords(
                    file, None, file_format
                )

                if not found_coords:
                    continue

                qprep_data = {
                    "atom_types": atom_types,
                    "cartesians": cartesians,
                    "charge": charge,
                    "mult": mult,
                    "name": name,
                }
                comfile = self.write(qprep_data)

                move_file(destination, self.args.w_dir_main, comfile)
                if create_dat:
                    self.args.log.write(f"o  {name} successfully processed at {destination}")

        if create_dat:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"\nTime QPREP: {elapsed_time} seconds\n")
            self.args.log.finalize()

    def sdf_2_com(self, sdf_file, destination, file_format):
        sdf_name = os.path.basename(Path(sdf_file)).split(".")[0]
        # get atom types, atomic coordinates, charge and multiplicity of all the mols in the SDF file

        low_check = None
        if self.args.lowest_only == True:
            low_check='lowest_only'
        if self.args.lowest_n is not None:
            low_check=int(self.args.lowest_n)
        if self.args.e_threshold_qprep is not None:
            low_check=float(self.args.e_threshold_qprep)
        mols = mol_from_sdf_or_mol_or_mol2(sdf_file, "qprep", self.args, low_check=low_check)

        for i, mol in enumerate(mols):
            (
                atom_types,
                cartesians,
                charge,
                mult,
                _,
            ) = self.qprep_coords(sdf_file, mol, file_format)

            if "_conf_" not in sdf_name:
                name_conf = f"{sdf_name}_conf_{i+1}"
            else:
                name_conf = sdf_name

            qprep_data = {
                "atom_types": atom_types,
                "cartesians": cartesians,
                "charge": charge,
                "mult": mult,
                "name": name_conf,
            }

            comfile = self.write(qprep_data)
            move_file(destination, self.args.w_dir_main, comfile)

    def get_header(self, qprep_data):
        """
        Gets the part of the input file above the molecular coordinates.
        """

        txt = ""
        name_file = add_prefix_suffix(qprep_data["name"], self.args)

        if self.args.program.lower() == "gaussian":
            if self.args.chk_path != '':
                txt += f'%chk={self.args.chk_path}\n'
            elif self.args.chk:
                txt += f'%chk={name_file}.chk\n'
            if self.args.oldchk_path != '':
                txt += f'%oldchk={self.args.oldchk_path}\n'
            elif self.args.oldchk:
                txt += f'%oldchk={name_file}.chk\n'
            txt += f"%nprocshared={self.args.nprocs}\n"
            txt += f"%mem={self.args.mem}\n"
            if self.args.qm_input[:2] not in ['p ','P ']:
                txt += f"# {self.args.qm_input}\n\n"
            else: # for #p in Gaussian inputs
                txt += f"#{self.args.qm_input}\n\n"
            txt += f'{name_file}\n\n'
            txt += f'{qprep_data["charge"]} {qprep_data["mult"]}\n'

        elif self.args.program.lower() == "orca":
            txt += f'# {name_file}\n'
            if "GB" in self.args.mem:
                mem_orca = int(self.args.mem.split("GB")[0]) * 1000
            elif "MB" in self.args.mem:
                mem_orca = self.args.mem.split("MB")[0]
            elif "MW" in self.args.mem:
                mem_orca = self.args.mem.split("MW")[0]
            else:
                mem_orca = self.args.mem
            if '%maxcore' not in self.args.qm_input:
                txt += f"%maxcore {mem_orca}\n"
            pal_included = False
            pal_list = ['%pal','pal1','pal3','pal3','pal4','pal5','pal6','pal7','pal8']
            for keyword in self.args.qm_input.split():
                if keyword.rstrip("\n").lower() in pal_list:
                    pal_included = True
            if not pal_included:
                txt += f"%pal nprocs {self.args.nprocs} end\n"
            txt += f"! {self.args.qm_input}\n"
            txt += f'* xyz {qprep_data["charge"]} {qprep_data["mult"]}\n'

        return txt


    def get_tail(self, qprep_data):
        """
        Gets the part of the input file below the molecular coordinates.
        """

        txt = ""
        # if the radius is modified for SMD, it has to be after the genecp info
        modifysph_line = "" 

        if self.args.program.lower() == "gaussian":
            # writes final section if selected
            if self.args.qm_end != "":

                qm_end_local = self.args.qm_end

                # check if the 'modifysph' line is in qm_end
                if "modifysph" in qm_end_local.lower():
                    end_lines = qm_end_local.split("\n")
                    for idx, line in enumerate(end_lines):
                        if "modifysph" in line.lower():
                            modifysph_idx = idx
                            break

                    # Remove empty lines after "modifysph"
                    while end_lines[modifysph_idx+1].strip() == "":
                        del end_lines[modifysph_idx+1]

                    modifysph_line = end_lines[modifysph_idx] + "\n\n" + end_lines[modifysph_idx+1] + "\n\n"
                    del end_lines[modifysph_idx:modifysph_idx+2]
                    qm_end_local = "\n".join(end_lines)

                txt += f"{qm_end_local}\n\n"

            if self.args.gen_atoms != [] and len(self.args.gen_atoms) > 0:
                # writes part for Gen/GenECP
                ecp_used, ecp_not_used, gen_type = [], [], "gen"
                if self.args.qm_input.lower().find("genecp") > -1:
                    gen_type = "genecp"

                for _, element_ecp in enumerate(qprep_data["atom_types"]):
                    if (
                        element_ecp in self.args.gen_atoms
                        and element_ecp not in ecp_used
                    ):
                        ecp_used.append(element_ecp)
                    elif (
                        element_ecp not in self.args.gen_atoms
                        and element_ecp not in ecp_not_used
                    ):
                        ecp_not_used.append(element_ecp)

                if len(ecp_not_used) > 0:
                    elements_not_used = " ".join([f"{sym}" for sym in ecp_not_used])
                    txt += f"{elements_not_used} 0\n{self.args.bs_nogen}\n****\n"
                if len(ecp_used) > 0:
                    elements_used = " ".join([f"{sym}" for sym in ecp_used])
                    txt += f"{elements_used} 0\n{self.args.bs_gen}\n****\n"

                if gen_type == "genecp" and len(ecp_used) > 0:
                    txt += "\n"
                    txt += f"{elements_used} 0\n{self.args.bs_gen}\n"

                txt += "\n"
                
        txt = txt.lstrip('\n')
        txt += modifysph_line
        return txt
        
    def write(self, qprep_data):

        if self.args.program.lower() == "gaussian":
            extension = "com"
        elif self.args.program.lower() == "orca":
            extension = "inp"

        name_file = add_prefix_suffix(qprep_data["name"], self.args)
        comfile = f'{name_file}.{extension}'

        if os.path.exists(self.args.w_dir_main / comfile):
            os.remove(self.args.w_dir_main / comfile)

        header = self.get_header(qprep_data)
        tail = self.get_tail(qprep_data)

        fileout = open(self.args.w_dir_main / comfile, "w")
        fileout.write(header)

        for atom_idx in range(0, len(qprep_data["atom_types"])):
            fileout.write(
                "{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}".format(
                    qprep_data["atom_types"][atom_idx],
                    qprep_data["cartesians"][atom_idx][0],
                    qprep_data["cartesians"][atom_idx][1],
                    qprep_data["cartesians"][atom_idx][2],
                )
            )
            if atom_idx != len(qprep_data["atom_types"]) - 1:
                fileout.write("\n")

        if self.args.program.lower() == "gaussian":
            fileout.write("\n\n")
        elif self.args.program.lower() == "orca":
            fileout.write("\n*")

        fileout.write(tail)
        fileout.close()

        return comfile

    def qprep_coords(self, file, mol, file_format):
        """
        Retrieve atom types and coordinates from multiple formats (LOG, OUT, JSON, MOL)
        """

        found_coords = False
        charge, mult = None, None
        if self.args.atom_types == [] or self.args.cartesians == []:
            if mol is not None:
                atom_types = []
                atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
                for atom, isotope in atom_data:
                    if isotope:
                        atom_types.append(atom.GetSymbol() + "(iso={})".format(isotope))
                    else:
                        atom_types.append(atom.GetSymbol())
                cartesians = mol.GetConformers()[0].GetPositions()
                try:
                    charge = int(mol.GetProp("Real charge"))
                except KeyError:
                    charge = Chem.GetFormalCharge(mol)
                try:
                    mult = int(mol.GetProp("Mult"))
                except KeyError:
                    NumRadicalElectrons = 0
                    for Atom in mol.GetAtoms():
                        NumRadicalElectrons += Atom.GetNumRadicalElectrons()
                    TotalElectronicSpin = NumRadicalElectrons / 2
                    mult = int((2 * TotalElectronicSpin) + 1)

            elif file_format in ["log", "out"]:
                # detect QM program and number of atoms
                if not self.args.command_line:
                    outlines = read_file(os.getcwd(), self.args.w_dir_main, file)
                else:
                    # if command lines are used, the program is already in that folder
                    outlines = read_file(os.getcwd(), os.getcwd(), file)
                n_atoms = 0
                resume_line = 0
                found_n_atoms = False

                for i, line in enumerate(outlines):
                    if line.find("Gaussian, Inc."):
                        program = "gaussian"
                        resume_line = i
                        break
                    elif line[i].find("O   R   C   A"):
                        program = "orca"
                        resume_line = i
                        break

                for i in range(resume_line, len(outlines)):
                    if program == "gaussian":
                        # get charge and mult
                        if outlines[i].find("Charge = ") > -1:
                            charge = int(outlines[i].split()[2])
                            mult = int(outlines[i].split()[5].rstrip("\n"))
                        # get number of atoms
                        elif outlines[i].find("Symbolic Z-matrix:") > -1:
                            for j in range(i + 2, len(outlines)):
                                if len(outlines[j].split()) > 0:
                                    n_atoms += 1
                                else:
                                    found_n_atoms = True
                                    break
                        elif found_n_atoms:
                            break

                atom_types, cartesians = QM_coords(outlines, -1, n_atoms, program, "")

            elif file_format == "json":
                with open(file) as json_file:
                    cclib_data = json.load(json_file)
                try:
                    atom_types, cartesians = cclib_atoms_coords(cclib_data)
                    charge = cclib_data["properties"]["charge"]
                    mult = cclib_data["properties"]["multiplicity"]
                except (AttributeError, KeyError):
                    atom_types, cartesians = [], []

        else:
            atom_types = self.args.atom_types
            cartesians = self.args.cartesians

        try:
            if atom_types == [] or cartesians == []:
                self.args.log.write(f"x  {file} does not contain coordinates and/or atom type information")
            if atom_types != [] and cartesians != []:
                found_coords = True
        except ValueError: # avoids an issue when comparing arrays with == []
            found_coords = True

        # overwrite with user-defined charge and multiplicity (if any)
        # or sets to default charge 0 and mult 1 if the parameters weren't found
        if self.args.charge is not None:
            charge = self.args.charge
        elif charge is None:
            charge = 0
        if self.args.mult is not None:
            mult = self.args.mult
        elif mult is None:
            mult = 1

        return atom_types, cartesians, charge, mult, found_coords


    def check_level_of_theory(self):
        """
        Cross check a chosen functional and basis set against a precompiled list of available options.
        Not necessarily a definitive list!
        """

        # read the predifined list of functionals and basis sets
        functional_csv = TEMPLATES_PATH / Path('functionals.csv')
        basis_set_csv = TEMPLATES_PATH / Path('basis_sets.csv')
        
        functional_data = pd.read_csv(functional_csv)
        functional_data.drop_duplicates(inplace=True)

        basis_set_data = pd.read_csv(basis_set_csv)
        basis_set_data.drop_duplicates(inplace=True)

        functional_list = functional_data[self.args.program].to_numpy().flatten()
        functional_list = [x for x in functional_list if str(x) != 'nan'] # remove NaN
        basis_set_list = basis_set_data[self.args.program].to_numpy().flatten()
        basis_set_list = [x for x in basis_set_list if str(x) != 'nan'] # remove NaN

        found_func, found_basis = False, False

        # first, look for the basis set from gen/genecp, both sets of basis sets used (i.e. for
        # gen atoms and for other atoms)
        if self.args.bs_gen != '':
            if self.args.bs_gen.upper() in (bs.upper() for bs in basis_set_list):
                if self.args.bs_nogen.upper() in (bs.upper() for bs in basis_set_list):
                    found_basis = True

        # for all the keywords in the qm_input option, check if there are compatible functionals and
        # basis sets
        for keyword in self.args.qm_input.split():
            for subkey in keyword.split('/'):
                if subkey.count('opt') == 0 and subkey.count('freq') == 0 and subkey.count('scrf') == 0 and subkey.count('pop') == 0 and subkey.count('gen') == 0:
                    if subkey.upper() in (func.upper() for func in functional_list):
                        found_func = True
                    if subkey.upper() in (bs.upper() for bs in basis_set_list):
                        found_basis = True

        if not found_func:
            self.args.log.write("x  WARNING! Verify that your functional is correct. If it is, please let us know to add it to the list of known functionals.")

        if not found_basis:
            self.args.log.write("x  WARNING! Verify that your basis set(s) is correct. If it is, please let us know to add it to the list of known basis sets.")