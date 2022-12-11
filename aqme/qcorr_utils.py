######################################################.
#        This file stores functions related          #
#               to the QCORR module                  #
######################################################.

import os
import glob
import pandas as pd
import json
import cclib
from pathlib import Path
from aqme.utils import move_file, read_file, Logger
import numpy as np

# Bondi VDW radii in Angstrom
bondi = {
    "H": 1.09,
    "He": 1.40,
    "Li": 1.81,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Ne": 1.54,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Ar": 1.88,
    "K": 2.75,
    "Ca": 2.31,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
    "Ga": 1.87,
    "Ge": 2.11,
    "As": 1.85,
    "Se": 1.90,
    "Br": 1.83,
    "Kr": 2.02,
    "Rb": 3.03,
    "Sr": 2.49,
    "Pd": 1.63,
    "Ag": 1.72,
    "Cd": 1.58,
    "In": 1.93,
    "Sn": 2.17,
    "Sb": 2.06,
    "Te": 2.06,
    "I": 1.98,
    "Xe": 2.16,
    "Cs": 3.43,
    "Ba": 2.68,
    "Pt": 1.72,
    "Au": 1.66,
    "Hg": 1.55,
    "Tl": 1.96,
    "Pb": 2.02,
    "Bi": 2.07,
    "Po": 1.97,
    "At": 2.02,
    "Rn": 2.20,
    "Fr": 3.48,
    "Ra": 2.83,
    "U": 1.86,
}

# covalent radii in Angstrom (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
rcov = {
    "H": 0.32,
    "He": 0.46,
    "Li": 1.33,
    "Be": 1.02,
    "B": 0.85,
    "C": 0.75,
    "N": 0.71,
    "O": 0.63,
    "F": 0.64,
    "Ne": 0.67,
    "Na": 1.55,
    "Mg": 1.39,
    "Al": 1.26,
    "Si": 1.16,
    "P": 1.11,
    "S": 1.03,
    "Cl": 0.99,
    "Ar": 0.96,
    "K": 1.96,
    "Ca": 1.71,
    "Sc": 1.48,
    "Ti": 1.36,
    "V": 1.34,
    "Cr": 1.22,
    "Mn": 1.19,
    "Fe": 1.16,
    "Co": 1.11,
    "Ni": 1.10,
    "Zn": 1.18,
    "Ga": 1.24,
    "Ge": 1.21,
    "As": 1.21,
    "Se": 1.16,
    "Br": 1.14,
    "Kr": 1.17,
    "Rb": 2.10,
    "Sr": 1.85,
    "Y": 1.63,
    "Zr": 1.54,
    "Nb": 1.47,
    "Mo": 1.38,
    "Tc": 1.28,
    "Ru": 1.25,
    "Rh": 1.25,
    "Pd": 1.20,
    "Ag": 1.28,
    "Cd": 1.36,
    "In": 1.42,
    "Sn": 1.40,
    "Sb": 1.40,
    "Te": 1.36,
    "I": 1.33,
    "Xe": 1.31,
}


def detect_linear(errortype, atom_types, cclib_data):
    """
    Check whether a linear molecule contain the right number of frequencies
    """

    linear_options_3 = [
        ["I", "I", "I"],
        ["N", "N", "N"],
        ["N", "C", "H"],
        ["O", "C", "O"],
        ["O", "C", "S"],
        ["S", "C", "S"],
        ["F", "Be", "F"],
        ["F", "Xe", "F"],
        ["O", "C", "N"],
        ["S", "C", "N"],
    ]
    linear_options_4 = [["C", "C", "H", "H"]]

    if len(atom_types) == 3:
        for linear_3 in linear_options_3:
            if (
                sorted(atom_types) == sorted(linear_3)
                and len(cclib_data["vibrations"]["frequencies"]) != 4
            ):
                errortype = "linear_mol_wrong"
                break
    elif len(atom_types) == 4:
        for linear_4 in linear_options_4:
            if (
                sorted(atom_types) == sorted(linear_4)
                and len(cclib_data["vibrations"]["frequencies"]) != 7
            ):
                errortype = "linear_mol_wrong"
                break

    return errortype


def full_check(w_dir_main=os.getcwd(), destination_fullcheck="", files="*.json", log=None):
    """
    Checks that multiple calculations were done following the same protocols, including
    program and version, grid size, level of theory, dispersion and solvation model.
    Parameters
    ----------
    w_dir_main : str
        Working directory
    destination_fullcheck : str
        Destination to create the file with the full check
    files : list of str
        json files to compare (glob.glob('*.json') and '*.json are both valid inputs to
        include all the json files from a folder)
    log : aqme.utils.Logger
        Logging instance where the status of the calculation will be written.
        If none provided it will default to aqme.utils.Logger('QCORR','fullcheck')
        and it will create the file QCORR_fullcheck.dat in the working directory.
    """

    if log is None: 
        log = Logger('QCORR','fullcheck')

    initial_dir = os.getcwd()
    w_dir_main = Path(w_dir_main)
    os.chdir(w_dir_main)

    if not isinstance(files, list):
        files = glob.glob(files)

    df_fullcheck = pd.DataFrame(
        columns=[
            "file",
            "program",
            "grid_type",
            "level_of_theory",
            "dispersion",
            "solvation",
        ]
    )

    for file in files:
        file_name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
        with open(file) as json_file:
            cclib_data = json.load(json_file)

        program = cclib_data["metadata"]["QM program"]
        solvation = cclib_data["metadata"]["solvation"]
        dispersion = cclib_data["metadata"]["dispersion model"]
        grid_type = cclib_data["metadata"]["grid type"]
        functional = cclib_data["metadata"]["functional"]
        bs = cclib_data["metadata"]["basis set"]
        if functional != "" or bs != "":
            level_of_theory = "/".join([functional, bs])
        else:
            level_of_theory = ""
        # designed to detect G4 calcs
        if level_of_theory == "HF/GFHFB2":
            level_of_theory = "G4"
        df_fullcheck.loc[len(df_fullcheck.index)] = [
            file_name,
            program,
            grid_type,
            level_of_theory,
            dispersion,
            solvation,
        ]

    fullcheck_file = "--QCORR_Fullcheck_Analysis--.dat"
    fullcheck_txt = "\n-- Full check analysis --"
    for prop in df_fullcheck.columns:
        if prop != "file":
            unique_props = df_fullcheck[prop].unique()
            if len(unique_props) > 1:
                fullcheck_txt += f"\nx  Different {prop} used in the calculations:"
                for unique_prop in unique_props:
                    file_names = df_fullcheck["file"].loc[
                        df_fullcheck[prop] == unique_prop
                    ]
                    fullcheck_txt += f"\n     * {unique_prop} in:"
                    for file_name in file_names:
                        adapted_name = file_name.replace("/", "\\").split("\\")[-1]
                        fullcheck_txt += f"\n       - {adapted_name}"
            else:
                fullcheck_txt += (
                    f"\no  Same {prop} ({unique_props[0]}) used in all the calculations"
                )

    fullcheck_analysis = open(fullcheck_file, "w")
    fullcheck_analysis.write(fullcheck_txt)
    fullcheck_analysis.close()
    log.write(fullcheck_txt)

    if destination_fullcheck == "":
        destination_fullcheck = w_dir_main.joinpath("success/json_files/")
    else:
        destination_fullcheck = Path(destination_fullcheck)
    move_file(destination_fullcheck, w_dir_main, fullcheck_file)

    os.chdir(initial_dir)


def check_isomerization(isom_data, file):
    """
    Inputs two molecules with the atoms in the same order and checks if any bond
    is too different between them.

    Bonds are considered when the distance between two atoms is smaller than
    either the sum of their adjusted VDW radii or covalent radii (dist = n*R1 + n*R2, where
    n is a user defined parameter).

    Bonds forming part of TSs are removed.

    Parameters
    ----------
    isom_data : dict
            Contains data related to coordinates, atoms and connectivity for input and output files
    file : string
            Filename (with extension)

    Returns
    -------
    isomerized : bool
            True if there is a clearly distorted bond within the geometries
    """

    # load connectivity matrix from the starting points and convert string into matrix
    if not isom_data["Initial csv"].empty:
        filename = file.replace("_" + file.split("_")[-1], "")
        init_connectivity_string = isom_data["Initial csv"][
            isom_data["Initial csv"]["code_name"] == filename
        ]["initial_connectiv"][0]
        init_connectivity = json.loads(
            init_connectivity_string.replace(".", ",")
            .replace(",]", "],")
            .replace("],]", "]]")
        )
        isom_data["Atoms input"] = init_connectivity[0]

    else:
        init_connectivity = gen_connectivity(
            isom_data, isom_data["Atoms input"], isom_data["Coords input"]
        )

    # in case the systems are not the same
    if len(isom_data["Atoms output"]) != len(isom_data["Atoms input"]):
        isomerized = True
    else:
        final_connectivity = gen_connectivity(
            isom_data, isom_data["Atoms output"], isom_data["Coords output"]
        )

        # check connectivity differences from initial structure
        diff = final_connectivity - init_connectivity

        # remove bonds involved in TSs from connectivity matrixes
        if not isom_data["Initial csv"].empty:
            if "TS_atom_idx" in isom_data["Initial csv"].columns:
                ts_atoms = isom_data["Initial csv"][
                    isom_data["Initial csv"]["code_name"] == filename
                ]["TS_atom_idx"][0].split(",")
                for i, ts_idx in enumerate(ts_atoms):
                    for j, ts_idx_2 in enumerate(ts_atoms):
                        if j > i:
                            diff[int(ts_idx)][int(ts_idx_2)] = 0
                            diff[int(ts_idx_2)][int(ts_idx)] = 0

        isomerized = np.any(diff)

    return isomerized


def gen_connectivity(isom_data, atom_types_conn, COORDINATES_conn):
    """
    Use VDW radii to infer a connectivity matrix
    """

    conn_mat = np.zeros((len(atom_types_conn), len(atom_types_conn)))
    for i, elem_i in enumerate(atom_types_conn):
        for j, elem_j in enumerate(atom_types_conn):
            if j > i:
                vdw_ij = bondi[elem_i] + bondi[elem_j]
                rcov_ij = rcov[elem_i] + rcov[elem_j]
                dist_ij = np.linalg.norm(
                    np.array(COORDINATES_conn[i]) - np.array(COORDINATES_conn[j])
                )
                if dist_ij / vdw_ij < float(
                    isom_data["VdW radii fraction"]
                ) or dist_ij / rcov_ij < float(isom_data["Covalent radii fraction"]):
                    conn_mat[i][j] = 1
                else:
                    pass

    return conn_mat


def get_json_data(self, file, cclib_data):
    """
    Get metadata and GoodVibes data for the json file (for older versions of cclib)
    """

    outlines = read_file(os.getcwd(), self.args.w_dir_main, file)

    # initial loop just to detect the QM program
    for i, line in enumerate(outlines):
        # get program
        if line.strip() == "Cite this work as:":
            cclib_data["metadata"] = {}
            qm_program = outlines[i + 1]

            cclib_data["metadata"]["QM program"] = qm_program[1:-2]
            for j in range(i, i + 60):
                if "**********" in outlines[j]:
                    run_date = outlines[j + 2].strip()
                    cclib_data["metadata"]["run date"] = run_date
                    break
            break

        elif "* O   R   C   A *" in line:
            for j in range(i, i + 100):
                if "Program Version" in line.strip():
                    cclib_data["metadata"] = {}
                    version_program = "ORCA version " + line.split()[2]
                    cclib_data["metadata"]["QM program"] = version_program
                    break

    if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:

        cclib_data["properties"]["rotational"] = {}
        for i, line in enumerate(outlines):
            # Extract memory
            if "%mem" in line:
                mem = line.strip().split("=")[-1]
                cclib_data["metadata"]["memory"] = mem

            # Extract number of processors
            elif "%nprocs" in line:
                nprocs = int(line.strip().split("=")[-1])
                cclib_data["metadata"]["processors"] = nprocs

            # Extract keywords line, solvation, dispersion and calculation type
            elif "#" in line and not hasattr(cclib_data, "keywords_line"):
                keywords_line = ""
                for j in range(i, i + 10):
                    if "----------" in outlines[j]:
                        break
                    else:
                        keywords_line += outlines[j].rstrip("\n")[1:]
                cclib_data["metadata"]["keywords line"] = keywords_line[2:]
                qm_solv, qm_disp = "gas_phase", "none"
                calc_type = "ground_state"
                calcfc_found, ts_found = False, False
                for keyword in keywords_line.split():
                    if keyword.lower().find("opt") > -1:
                        if keyword.lower().find("calcfc") > -1:
                            calcfc_found = True
                        if keyword.lower().find("ts") > -1:
                            ts_found = True
                    elif keyword.lower().startswith("scrf"):
                        qm_solv = keyword
                    elif keyword.lower().startswith("emp"):
                        qm_disp = keyword
                if calcfc_found and ts_found:
                    calc_type = "transition_state"
                cclib_data["metadata"]["solvation"] = qm_solv
                cclib_data["metadata"]["dispersion model"] = qm_disp
                cclib_data["metadata"]["ground or transition state"] = calc_type

            # Basis set name
            elif line[1:15] == "Standard basis":
                cclib_data["metadata"]["basis set"] = line.split()[2]

            # functional
            if not hasattr(cclib_data, "BOMD") and line[1:9] == "SCF Done":
                t1 = line.split()[2]
                if t1 == "E(RHF)":
                    cclib_data["metadata"]["functional"] = "HF"
                else:
                    cclib_data["metadata"]["functional"] = t1[
                        t1.index("(") + 2 : t1.rindex(")")
                    ]
                break

            # Extract grid type
            elif line[1:8] == "ExpMin=":
                grid_lookup = {
                    1: "sg1",
                    2: "coarse",
                    4: "fine",
                    5: "ultrafine",
                    7: "superfine",
                }
                IRadAn = int(line.strip().split()[-3])
                grid = grid_lookup[IRadAn]
                cclib_data["metadata"]["grid type"] = grid

        # Keeps track of convergence during Freq calcs
        if "optimization" in cclib_data:
            cclib_data["optimization"]["times converged"] = 1

        # Extract <S**2> before and after spin annihilation, energy, and convergence in freq calc
        for i in reversed(range(0, len(outlines) - 30)):
            # For time dependent (TD) calculations
            if "E(TD-HF/TD-DFT)" in outlines[i]:
                td_e = float(line.strip().split()[-1])
                cclib_data["properties"]["energy"][
                    "TD energy"
                ] = cclib.parser.utils.convertor(td_e, "hartree", "eV")

            # For G4 calculations look for G4 energies (Gaussian16a bug prints G4(0 K) as DE(HF)) --Brian modified to work for G16c-where bug is fixed.
            elif line.strip().startswith("E(ZPE)="):  # Overwrite DFT ZPE with G4 ZPE
                zero_point_corr = float(line.strip().split()[1])
            elif line.strip().startswith("G4(0 K)"):
                G4_energy = float(line.strip().split()[2])
                G4_energy -= zero_point_corr  # Remove G4 ZPE
                cclib_data["properties"]["energy"][
                    "G4 energy"
                ] = cclib.parser.utils.convertor(G4_energy, "hartree", "eV")

            # For ONIOM calculations use the extrapolated value rather than SCF value
            elif "ONIOM: extrapolated energy" in line.strip():
                oniom_e = float(line.strip().split()[4])
                cclib_data["properties"]["energy"][
                    "ONIOM energy"
                ] = cclib.parser.utils.convertor(oniom_e, "hartree", "eV")

            elif "S**2 before annihilation" in outlines[i]:
                cclib_data["properties"]["S2 after annihilation"] = float(
                    outlines[i].strip().split()[-1]
                )
                cclib_data["properties"]["S2 before annihilation"] = float(
                    outlines[i].strip().split()[-3][:-1]
                )

            # Extract symmetry point group
            elif "Full point group" in outlines[i]:
                point_group = outlines[i].strip().split()[3]
                cclib_data["properties"]["rotational"][
                    "symmetry point group"
                ] = point_group
                break

            elif "Stationary point found" in outlines[i]:
                cclib_data["optimization"]["times converged"] = 2

            # Extract symmetry number, rotational constants and rotational temperatures
            elif "Rotational symmetry number" in outlines[i]:
                symmno = int(outlines[i].strip().split()[3].split(".")[0])
                cclib_data["properties"]["rotational"]["symmetry number"] = symmno

            elif outlines[i].find("Rotational constants (GHZ):") > -1:
                try:
                    roconst = [
                        float(outlines[i].strip().replace(":", " ").split()[3]),
                        float(outlines[i].strip().replace(":", " ").split()[4]),
                        float(outlines[i].strip().replace(":", " ").split()[5]),
                    ]
                except ValueError:
                    if outlines[i].find("********") > -1:
                        roconst = [
                            float(outlines[i].strip().replace(":", " ").split()[4]),
                            float(outlines[i].strip().replace(":", " ").split()[5]),
                        ]
                cclib_data["properties"]["rotational"]["rotational constants"] = roconst

            elif outlines[i].find("Rotational temperature ") > -1:
                rotemp = [float(outlines[i].strip().split()[3])]
                cclib_data["properties"]["rotational"][
                    "rotational temperatures"
                ] = rotemp

            elif outlines[i].find("Rotational temperatures") > -1:
                try:
                    rotemp = [
                        float(outlines[i].strip().split()[3]),
                        float(outlines[i].strip().split()[4]),
                        float(outlines[i].strip().split()[5]),
                    ]
                except ValueError:
                    if outlines[i].find("********") > -1:
                        rotemp = [
                            float(outlines[i].strip().split()[4]),
                            float(outlines[i].strip().split()[5]),
                        ]
                cclib_data["properties"]["rotational"][
                    "rotational temperatures"
                ] = rotemp

            elif outlines[i].find("SCF GIAO Magnetic shielding tensor (ppm)") > -1:
                nmr_iso = []
                nmr_anis = []
                nmr_eigen = []
                cclib_data["properties"]["NMR"] = {}
                for j in range(i, len(outlines)):
                    if outlines[j].find("Isotropic") > -1:
                        nmr_iso.append(float(outlines[j].split()[4]))
                        nmr_anis.append(float(outlines[j].split()[7]))
                    elif outlines[j].find("Eigenvalues") > -1:
                        nmr_eigen.append(
                            [
                                float(outlines[j].split()[1]),
                                float(outlines[j].split()[2]),
                                float(outlines[j].split()[3]),
                            ]
                        )
                    elif outlines[j].find("*************************") > -1:
                        break
                cclib_data["properties"]["NMR"]["NMR anisotopic tensors"] = nmr_anis
                cclib_data["properties"]["NMR"]["NMR eigenvalues"] = nmr_eigen
                cclib_data["properties"]["NMR"]["NMR isotopic tensors"] = nmr_iso


    elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
        for i in reversed(range(0, outlines)):
            if outlines[i][:25] == "FINAL SINGLE POINT ENERGY":
                # in eV to match the format from cclib
                orca_e = float(outlines[i].split()[-1])
                cclib_data["properties"]["energy"][
                    "final single point energy"
                ] = cclib.parser.utils.convertor(orca_e, "hartree", "eV")
                break

    if cclib_data != {}:
        with open(f'{file.split(".")[0]}.json', "w") as outfile:
            json.dump(cclib_data, outfile, indent=1)

    return cclib_data
