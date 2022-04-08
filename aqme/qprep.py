######################################################.
#        This file stores the QPREP class            #
######################################################.

import os
import sys
import time
from rdkit import Chem
import json
from aqme.utils import cclib_atoms_coords, QM_coords, read_file
from aqme.utils import move_file, load_variables
from pathlib import Path
from aqme.crest import xyzall_2_xyz
import glob
import subprocess

try:
    import pybel
except ImportError:
    from openbabel import pybel  # for openbabel>=3.0.0


class qprep:
    """
    Class containing all the functions from the QPREP module related to Gaussian input files
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "qprep")

        if self.args.destination is None:
            self.args.destination = self.args.w_dir_main.joinpath("QCALC")
        else:
            self.args.destination = Path(self.args.destination)

        if self.args.qm_input == "":
            print("x  No keywords line was specified! (i.e. qm_input=KEYWORDS_LINE).")
            self.args.log.write(
                "x  No keywords line was specified! (i.e. qm_input=KEYWORDS_LINE)."
            )
            sys.exit()

        if self.args.program not in ["gaussian", "orca"]:
            print("x  Check the program specified! (i.e. program=gaussian).")
            self.args.log.write(
                "x  Check the program specified! (i.e. program=gaussian)."
            )
            sys.exit()

        # write input files
        os.chdir(self.args.w_dir_main)
        for file in self.args.files:
            found_coords = False
            name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
            if file.split(".")[1].lower() == "sdf":
                try:
                    # get atom types, atomic coordinates, charge and multiplicity of all the mols in the SDF file
                    mols = Chem.SDMolSupplier(str(file), removeHs=False)
                    for i, mol in enumerate(mols):
                        atom_types, cartesians, charge, mult, _ = self.qprep_coords(
                            file, found_coords, mol
                        )
                        if charge == None:
                            charge = 0
                        if mult == None:
                            mult = 1
                        name_conf = f"{name}_conf_{i+1}"
                        qprep_data = {
                            "atom_types": atom_types,
                            "cartesians": cartesians,
                            "charge": charge,
                            "mult": mult,
                            "name": name_conf,
                        }
                        comfile = self.write(qprep_data)
                        move_file(self.args.destination, self.args.w_dir_main, comfile)

                except OSError:
                    print(f"x  {name} couldn't be processed!")
                    self.args.log.write(f"x  {name} couldn't be processed!")
                    continue

                print(f"o  {name} successfully processed")
                self.args.log.write(f"o  {name} successfully processed")

            elif file.split(".")[1].lower() == "xyz":
                xyzall_2_xyz(file, name)
                xyz_files = glob.glob(name + "_conf_*.xyz")
                for i, file1 in enumerate(xyz_files):
                    name1 = file1.split(".xyz")[0]
                    command_xyz = [
                        "obabel",
                        "-ixyz",
                        file1,
                        "-osdf",
                        "-O" + name1 + ".sdf",
                    ]
                    subprocess.call(command_xyz)
                    os.remove(file1)

                sdf_files = glob.glob(name + "*.sdf")
                for i, file2 in enumerate(sdf_files):
                    mol = Chem.SDMolSupplier(file2, removeHs=False)[0]
                    atom_types, cartesians, charge, mult, _ = self.qprep_coords(
                        file, found_coords, mol
                    )
                    if charge == None:
                        charge = 0
                    if mult == None:
                        mult = 1
                    name_conf = f"{name}_conf_{i+1}"
                    qprep_data = {
                        "atom_types": atom_types,
                        "cartesians": cartesians,
                        "charge": charge,
                        "mult": mult,
                        "name": name_conf,
                    }
                    comfile = self.write(qprep_data)
                    move_file(self.args.destination, self.args.w_dir_main, comfile)
                    os.remove(file2)

                print(f"o  {name} successfully processed")
                self.args.log.write(f"o  {name} successfully processed")

            else:

                found_coords = False
                atom_types, cartesians, charge, mult, found_coords = self.qprep_coords(
                    file, found_coords, None
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
                move_file(self.args.destination, self.args.w_dir_main, comfile)

                print(f"o  {name} successfully processed")
                self.args.log.write(f"o  {name} successfully processed")

        # this is added to avoid path problems in jupyter notebooks
        os.chdir(self.args.initial_dir)

        elapsed_time = round(time.time() - start_time_overall, 2)
        print(f"\nTime QPREP: {elapsed_time} seconds\n")
        self.args.log.write(f"\nTime QPREP: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def get_header(self, qprep_data):
        """
        Gets the part of the input file above the molecular coordinates.
        """

        txt = ""

        if self.args.program.lower() == "gaussian":
            if self.args.chk:
                txt += f'%chk={qprep_data["name"]}.chk\n'
            txt += f"%nprocshared={self.args.nprocs}\n"
            txt += f"%mem={self.args.mem}\n"
            txt += f"# {self.args.qm_input}"
            txt += "\n\n"
            txt += f'{qprep_data["name"]}\n\n'
            txt += f'{qprep_data["charge"]} {qprep_data["mult"]}\n'

        elif self.args.program.lower() == "orca":
            txt += f'# {qprep_data["name"]}\n'
            if self.args.mem.find("GB"):
                mem_orca = int(self.args.mem.split("GB")[0]) * 1000
            elif self.args.mem.find("MB"):
                mem_orca = self.args.mem.split("MB")[0]
            elif self.args.args.mem.find("MW"):
                mem_orca = self.args.mem.split("MW")[0]
            txt += f"%maxcore {mem_orca}\n"
            txt += f"%pal nprocs {self.args.nprocs} end\n"
            txt += f"! {self.args.qm_input}\n"
            txt += f'* xyz {qprep_data["charge"]} {qprep_data["mult"]}\n'

        return txt

    def get_tail(self, qprep_data):
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

    def write(self, qprep_data):

        if self.args.program.lower() == "gaussian":
            extension = "com"
        elif self.args.program.lower() == "orca":
            extension = "inp"
        if self.args.suffix != "":
            comfile = f'{qprep_data["name"]}_{self.args.suffix}.{extension}'
        else:
            comfile = f'{qprep_data["name"]}.{extension}'

        if os.path.exists(comfile):
            os.remove(comfile)

        header = self.get_header(qprep_data)
        tail = self.get_tail(qprep_data)

        fileout = open(comfile, "w")
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

    def qprep_coords(self, file, found_coords, mol):
        """
        Retrieve atom types and coordinates from multiple formats (LOG, OUT, JSON, MOL)
        """

        charge, mult = None, None
        if self.args.atom_types == [] or self.args.cartesians == []:
            if mol is not None:
                atom_types = [atom.GetSymbol() for _, atom in enumerate(mol.GetAtoms())]
                cartesians = mol.GetConformers()[0].GetPositions()
                try:
                    charge = int(mol.GetProp("Real charge"))
                except KeyError:
                    pass
                try:
                    mult = int(mol.GetProp("Mult"))
                except KeyError:
                    pass

            elif file.split(".")[1] in ["log", "out"]:
                # detect QM program and number of atoms
                outlines = read_file(self.args.w_dir_main, file)
                n_atoms = 0
                resume_line = 0

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
                            symbol_z_line = i
                        elif outlines[i].find("GradGrad") > -1:
                            gradgrad_line = i
                            n_atoms = gradgrad_line - symbol_z_line - 4
                            break

                atom_types, cartesians = QM_coords(outlines, 0, n_atoms, program)

            elif file.split(".")[1] == "json":
                with open(file) as json_file:
                    cclib_data = json.load(json_file)
                try:
                    atom_types, cartesians = cclib_atoms_coords(cclib_data)
                    charge = cclib_data["properties"]["charge"]
                    mult = cclib_data["properties"]["multiplicity"]
                except (AttributeError, KeyError):
                    print(
                        f"x  {file} does not contain coordinates and/or atom type information"
                    )
                    atom_types, cartesians = [], []

        else:
            atom_types = self.args.atom_types
            cartesians = self.args.cartesians

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

        if atom_types != [] and cartesians != []:
            found_coords = True

        return atom_types, cartesians, charge, mult, found_coords
