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
    "H": 1.09, "He": 1.40, "Li": 1.81, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
    "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80,
    "Cl": 1.75, "Ar": 1.88, "K": 2.75, "Ca": 2.31, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87,
    "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.83, "Kr": 2.02, "Rb": 3.03, "Sr": 2.49, "Pd": 1.63,
    "Ag": 1.72, "Cd": 1.58, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "Pt": 1.72, "Au": 1.66, "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
    "Po": 1.97, "At": 2.02, "Rn": 2.20, "Fr": 3.48, "Ra": 2.83, "U": 1.86,
}

# Covalent radii in Angstrom (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
rcov = {
    "H": 0.32, "He": 0.46, "Li": 1.33, "Be": 1.02, "B": 0.85, "C": 0.75, "N": 0.71, "O": 0.63,
    "F": 0.64, "Ne": 0.67, "Na": 1.55, "Mg": 1.39, "Al": 1.26, "Si": 1.16, "P": 1.11, "S": 1.03,
    "Cl": 0.99, "Ar": 0.96, "K": 1.96, "Ca": 1.71, "Sc": 1.48, "Ti": 1.36, "V": 1.34, "Cr": 1.22,
    "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Zn": 1.18, "Ga": 1.24, "Ge": 1.21, "As": 1.21,
    "Se": 1.16, "Br": 1.14, "Kr": 1.17, "Rb": 2.10, "Sr": 1.85, "Y": 1.63, "Zr": 1.54, "Nb": 1.47,
    "Mo": 1.38, "Tc": 1.28, "Ru": 1.25, "Rh": 1.25, "Pd": 1.20, "Ag": 1.28, "Cd": 1.36, "In": 1.42,
    "Sn": 1.40, "Sb": 1.40, "Te": 1.36, "I": 1.33, "Xe": 1.31,
}


def _get_linear_molecule_templates():
    """Get templates for known linear molecules.
    
    Returns:
        tuple: (linear_options_3, linear_options_4) - atom combinations for 3 and 4 atom linear molecules
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
    
    return linear_options_3, linear_options_4


def _check_linear_3atom(atom_types, cclib_data, linear_options_3):
    """Check if 3-atom molecule is linear with correct frequency count.
    
    Linear triatomic molecules should have 4 vibrational modes (3N-5).
    
    Args:
        atom_types (list): Element symbols for atoms
        cclib_data (dict): Parsed cclib data with vibfreqs
        linear_options_3 (list): Known linear 3-atom combinations
    
    Returns:
        bool: True if linear molecule has wrong frequency count
    """
    for linear_3 in linear_options_3:
        if (sorted(atom_types) == sorted(linear_3) and 
            len(cclib_data["vibfreqs"]) != 4):
            return True
    return False


def _check_linear_4atom(atom_types, cclib_data, linear_options_4):
    """Check if 4-atom molecule is linear with correct frequency count.
    
    Linear 4-atom molecules should have 7 vibrational modes (3N-5).
    
    Args:
        atom_types (list): Element symbols for atoms
        cclib_data (dict): Parsed cclib data with vibfreqs
        linear_options_4 (list): Known linear 4-atom combinations
    
    Returns:
        bool: True if linear molecule has wrong frequency count
    """
    for linear_4 in linear_options_4:
        if (sorted(atom_types) == sorted(linear_4) and 
            len(cclib_data["vibfreqs"]) != 7):
            return True
    return False


def detect_linear(errortype, atom_types, cclib_data):
    """Check whether a linear molecule has the correct number of frequencies.
    
    Linear molecules have 3N-5 vibrational modes (vs 3N-6 for non-linear).
    Validates that known linear molecules have expected frequency count.
    
    Args:
        errortype (str): Current error type classification
        atom_types (list): Element symbols for atoms in molecule
        cclib_data (dict): Parsed cclib data containing vibfreqs
    
    Returns:
        str: Updated errortype ('linear_mol_wrong' if incorrect, otherwise unchanged)
    """
    linear_options_3, linear_options_4 = _get_linear_molecule_templates()
    
    if len(atom_types) == 3:
        if _check_linear_3atom(atom_types, cclib_data, linear_options_3):
            errortype = "linear_mol_wrong"
    
    elif len(atom_types) == 4:
        if _check_linear_4atom(atom_types, cclib_data, linear_options_4):
            errortype = "linear_mol_wrong"
    
    return errortype


def _load_json_files(files):
    """Load and parse JSON files from directory.
    
    Args:
        files (str or list): File pattern or list of files
    
    Returns:
        list: List of file paths
    """
    if not isinstance(files, list):
        files = glob.glob(files)
    return files


def _extract_file_metadata(file):
    """Extract metadata from a single JSON file.
    
    Args:
        file (str): Path to JSON file
    
    Returns:
        dict: Metadata containing file, program, grid_type, level_of_theory, dispersion, solvation
    """
    file_name = os.path.basename(Path(file)).split(".")[0]
    
    with open(file) as json_file:
        cclib_data = json.load(json_file)
    
    program = cclib_data["metadata"]["QM program"]
    solvation = cclib_data["metadata"]["solvation"]
    dispersion = cclib_data["metadata"]["dispersion model"]
    grid_type = cclib_data["metadata"]["grid type"]
    functional = cclib_data["metadata"]["functional"]
    bs = cclib_data["metadata"]["basis set"]
    
    # Construct level of theory
    if functional != "" or bs != "":
        level_of_theory = "/".join([functional, bs])
    else:
        level_of_theory = ""
    
    # Detect G4 calculations
    if level_of_theory == "HF/GFHFB2":
        level_of_theory = "G4"
    
    return {
        "file": file_name,
        "program": program,
        "grid_type": grid_type,
        "level_of_theory": level_of_theory,
        "dispersion": dispersion,
        "solvation": solvation,
    }


def _build_fullcheck_dataframe(files):
    """Build DataFrame with metadata from all JSON files.
    
    Args:
        files (list): List of JSON file paths
    
    Returns:
        pd.DataFrame: DataFrame with file metadata
    """
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
        metadata = _extract_file_metadata(file)
        df_fullcheck.loc[len(df_fullcheck.index)] = [
            metadata["file"],
            metadata["program"],
            metadata["grid_type"],
            metadata["level_of_theory"],
            metadata["dispersion"],
            metadata["solvation"],
        ]
    
    return df_fullcheck


def _generate_fullcheck_report(df_fullcheck):
    """Generate full check analysis report text.
    
    Args:
        df_fullcheck (pd.DataFrame): DataFrame with file metadata
    
    Returns:
        str: Formatted report text
    """
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
    
    return fullcheck_txt


def _write_and_move_report(fullcheck_txt, w_dir_main, destination_fullcheck, log):
    """Write fullcheck report and move to destination.
    
    Args:
        fullcheck_txt (str): Report text content
        w_dir_main (Path): Working directory
        destination_fullcheck (str): Destination folder
        log (Logger): Logger instance
    """
    fullcheck_file = "--QCORR_Fullcheck_Analysis--.dat"
    
    with open(fullcheck_file, "w") as fullcheck_analysis:
        fullcheck_analysis.write(fullcheck_txt)
    
    log.write(fullcheck_txt)
    
    if destination_fullcheck == "":
        destination_fullcheck = w_dir_main.joinpath("success/json_files/")
    else:
        destination_fullcheck = Path(destination_fullcheck)
    
    move_file(destination_fullcheck, w_dir_main, fullcheck_file)


def full_check(w_dir_main=os.getcwd(), destination_fullcheck="", files="*.json", log=None):
    """Check that calculations follow consistent protocols.
    
    Validates that multiple calculations use the same program, version,
    grid size, level of theory, dispersion, and solvation model.
    
    Parameters
    ----------
    w_dir_main : str or Path
        Working directory containing JSON files
    destination_fullcheck : str
        Destination folder for the fullcheck report file
    files : str or list
        JSON file pattern (e.g., '*.json') or list of file paths
    log : aqme.utils.Logger, optional
        Logger for output. If None, creates default logger.
    """
    if log is None:
        log = Logger('QCORR', 'fullcheck')
    
    initial_dir = os.getcwd()
    w_dir_main = Path(w_dir_main)
    os.chdir(w_dir_main)
    
    # Load files and build dataframe
    files = _load_json_files(files)
    df_fullcheck = _build_fullcheck_dataframe(files)
    
    # Generate and write report
    fullcheck_txt = _generate_fullcheck_report(df_fullcheck)
    _write_and_move_report(fullcheck_txt, w_dir_main, destination_fullcheck, log)
    
    os.chdir(initial_dir)


def _load_initial_connectivity(isom_data, file):
    """Load initial connectivity matrix from CSV or generate it.
    
    Args:
        isom_data (dict): Isomerization data including Initial csv
        file (str): Filename with extension
    
    Returns:
        np.ndarray: Initial connectivity matrix
    """
    if not isom_data["Initial csv"].empty:
        filename = file.replace("_" + file.split("_")[-1], "")
        init_connectivity_string = isom_data["Initial csv"][
            isom_data["Initial csv"]["code_name"] == filename
        ]["initial_connectiv"][0]
        
        # Parse JSON string to matrix
        init_connectivity = json.loads(
            init_connectivity_string.replace(".", ",")
            .replace(",]", "],")
            .replace("],]", "]]")
        )
        isom_data["Atoms input"] = init_connectivity[0]
        return np.array(init_connectivity)
    else:
        return gen_connectivity(
            isom_data, isom_data["Atoms input"], isom_data["Coords input"]
        )


def _remove_ts_bonds_from_diff(diff, isom_data, file):
    """Remove bonds involved in transition states from connectivity difference matrix.
    
    Args:
        diff (np.ndarray): Connectivity difference matrix
        isom_data (dict): Isomerization data including Initial csv
        file (str): Filename with extension
    
    Returns:
        np.ndarray: Modified difference matrix with TS bonds zeroed
    """
    if not isom_data["Initial csv"].empty:
        if "TS_atom_idx" in isom_data["Initial csv"].columns:
            filename = file.replace("_" + file.split("_")[-1], "")
            ts_atoms = isom_data["Initial csv"][
                isom_data["Initial csv"]["code_name"] == filename
            ]["TS_atom_idx"][0].split(",")
            
            for i, ts_idx in enumerate(ts_atoms):
                for j, ts_idx_2 in enumerate(ts_atoms):
                    if j > i:
                        diff[int(ts_idx)][int(ts_idx_2)] = 0
                        diff[int(ts_idx_2)][int(ts_idx)] = 0
    
    return diff


def check_isomerization(isom_data, file):
    """Check if molecule has isomerized by comparing connectivity matrices.
    
    Compares input and output geometries using connectivity matrices based
    on VDW and covalent radii. Bonds in transition states are excluded.
    
    Parameters
    ----------
    isom_data : dict
        Contains coordinates, atoms, and connectivity data for input/output
    file : str
        Filename with extension
    
    Returns
    -------
    bool
        True if geometry has isomerized (different connectivity)
    """
    # Load initial connectivity matrix
    init_connectivity = _load_initial_connectivity(isom_data, file)
    
    # Check for atom count mismatch
    if len(isom_data["Atoms output"]) != len(isom_data["Atoms input"]):
        return True
    
    # Generate final connectivity matrix
    final_connectivity = gen_connectivity(
        isom_data, isom_data["Atoms output"], isom_data["Coords output"]
    )
    
    # Calculate connectivity differences
    diff = final_connectivity - init_connectivity
    
    # Remove TS bonds from consideration
    diff = _remove_ts_bonds_from_diff(diff, isom_data, file)
    
    # Check if any connectivity changed
    isomerized = np.any(diff)
    
    return isomerized


def _ensure_radii_defined(atom_types_conn):
    """Ensure all atoms have defined VDW and covalent radii.
    
    Assigns default radius of 1.0 Angstrom for atoms without defined radii.
    
    Args:
        atom_types_conn (list): Element symbols for atoms
    """
    for atom in atom_types_conn:
        if atom not in bondi:
            bondi[atom] = 1.0
        if atom not in rcov:
            rcov[atom] = 1.0
    
    return bondi, rcov


def _calculate_bond_distance(coords_i, coords_j):
    """Calculate Euclidean distance between two atoms.
    
    Args:
        coords_i (array): Coordinates of first atom
        coords_j (array): Coordinates of second atom
    
    Returns:
        float: Distance in Angstroms
    """
    return np.linalg.norm(np.array(coords_i) - np.array(coords_j))


def _is_bonded(elem_i, elem_j, dist_ij, vdwfrac, covfrac, atom_types_conn):
    """Determine if two atoms are bonded based on distance thresholds.
    
    Atoms are considered bonded if distance is less than either:
    - vdwfrac * (VDW_i + VDW_j), or
    - covfrac * (cov_i + cov_j)
    
    Args:
        elem_i (str): Element symbol of first atom
        elem_j (str): Element symbol of second atom
        dist_ij (float): Distance between atoms
        vdwfrac (float): VDW radii fraction threshold
        covfrac (float): Covalent radii fraction threshold
        atom_types_conn (list): Element symbols for each atom

    Returns:
        bool: True if atoms are bonded
    """
    # Ensure all atoms have defined radii
    bondi, rcov = _ensure_radii_defined(atom_types_conn)

    vdw_ij = bondi[elem_i] + bondi[elem_j]
    rcov_ij = rcov[elem_i] + rcov[elem_j]
    
    return (dist_ij / vdw_ij < vdwfrac or 
            dist_ij / rcov_ij < covfrac)


def gen_connectivity(isom_data, atom_types_conn, COORDINATES_conn):
    """Generate connectivity matrix using VDW and covalent radii.
    
    Creates symmetric connectivity matrix where 1 indicates bonded atoms.
    Bonds are detected when distance is below either VDW or covalent
    radii thresholds (multiplied by user-defined fractions).
    
    Parameters
    ----------
    isom_data : dict
        Contains 'VdW radii fraction' and 'Covalent radii fraction' parameters
    atom_types_conn : list
        Element symbols for each atom
    COORDINATES_conn : list
        3D coordinates for each atom
    
    Returns
    -------
    np.ndarray
        Symmetric connectivity matrix (upper triangular populated)
    """
    # Extract threshold parameters
    vdwfrac = float(isom_data["VdW radii fraction"])
    covfrac = float(isom_data["Covalent radii fraction"])
    
    # Initialize connectivity matrix
    conn_mat = np.zeros((len(atom_types_conn), len(atom_types_conn)))
    
    # Fill upper triangle of connectivity matrix
    for i, elem_i in enumerate(atom_types_conn):
        for j, elem_j in enumerate(atom_types_conn):
            if j > i:
                dist_ij = _calculate_bond_distance(
                    COORDINATES_conn[i], 
                    COORDINATES_conn[j]
                )
                
                if _is_bonded(elem_i, elem_j, dist_ij, vdwfrac, covfrac, atom_types_conn):
                    conn_mat[i][j] = 1
    
    return conn_mat


def _detect_qm_program(outlines):
    """Detect QM program and version from output file.
    
    Args:
        outlines (list): Lines from output file
    
    Returns:
        dict: Metadata dictionary with QM program and run date
    """
    metadata = {}
    
    for i, line in enumerate(outlines):
        # Gaussian detection
        if line.strip() == "Cite this work as:":
            qm_program = outlines[i + 1]
            metadata["QM program"] = qm_program[1:-2]
            
            for j in range(i, i + 60):
                if "**********" in outlines[j]:
                    run_date = outlines[j + 2].strip()
                    metadata["run date"] = run_date
                    break
            break
        
        # ORCA detection
        elif "* O   R   C   A *" in line:
            for j in range(i, i + 100):
                if "Program Version" in outlines[j].strip():
                    version_program = "ORCA version " + outlines[j].split()[2]
                    metadata["QM program"] = version_program
                    break
            break
    
    return metadata


def _extract_gaussian_keywords(outlines, i):
    """Extract keywords line from Gaussian output.
    
    Args:
        outlines (list): Lines from output file
        i (int): Starting line index
    
    Returns:
        str: Complete keywords line
    """
    keywords_line = ""
    for j in range(i, i + 10):
        if "----------" in outlines[j]:
            break
        else:
            keywords_line += outlines[j].rstrip("\n")[1:]
    return keywords_line[2:]


def _parse_gaussian_keywords(keywords_line, cclib_data):
    """Parse Gaussian keywords for solvation, dispersion, and calculation type.
    
    Args:
        keywords_line (str): Gaussian keywords line
        cclib_data (dict): cclib data dictionary to update
    """
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
        elif keyword == 'gen' or 'gen' in keyword.split('/'):
            cclib_data["metadata"]["basis set"] = 'gen'
        elif keyword == 'genecp' or 'genecp' in keyword.split('/'):
            cclib_data["metadata"]["basis set"] = 'genecp'
    
    if calcfc_found and ts_found:
        calc_type = "transition_state"
    
    cclib_data["metadata"]["solvation"] = qm_solv
    cclib_data["metadata"]["dispersion model"] = qm_disp
    cclib_data["metadata"]["ground or transition state"] = calc_type


def _extract_gaussian_metadata(cclib_data, outlines):
    """Extract Gaussian-specific metadata from output file.
    
    Args:
        cclib_data (dict): cclib data dictionary to update
        outlines (list): Lines from output file
    """
    cclib_data["rotational"] = {}
    
    for i, line in enumerate(outlines):
        # Extract memory
        if "%mem" in line:
            mem = line.strip().split("=")[-1]
            cclib_data["metadata"]["memory"] = mem
        
        # Extract number of processors
        elif "%nprocs" in line:
            nprocs = int(line.strip().split("=")[-1])
            cclib_data["metadata"]["processors"] = nprocs
        
        # Extract keywords line
        elif "#" in line and "keywords line" not in cclib_data["metadata"]:
            keywords_line = _extract_gaussian_keywords(outlines, i)
            cclib_data["metadata"]["keywords line"] = keywords_line
            _parse_gaussian_keywords(keywords_line, cclib_data)
        
        # Basis set name
        elif line[1:15] == "Standard basis":
            cclib_data["metadata"]["basis set"] = line.split()[2]
        
        # Functional
        if not hasattr(cclib_data, "BOMD") and line[1:9] == "SCF Done":
            t1 = line.split()[2]
            if t1 == "E(RHF)":
                cclib_data["metadata"]["functional"] = "HF"
            else:
                cclib_data["metadata"]["functional"] = t1[
                    t1.index("(") + 2 : t1.rindex(")")
                ]
        
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
        
        if "functional" in cclib_data["metadata"] and "grid type" in cclib_data["metadata"]:
            break
    
    # Track convergence during Freq calcs
    if "optdone" in cclib_data:
        if cclib_data["optdone"] in [True, 'true']:
            cclib_data["opt times converged"] = 1


def _extract_gaussian_properties(cclib_data, outlines):
    """Extract energy, S**2, and rotational data from Gaussian output.
    
    Args:
        cclib_data (dict): cclib data dictionary to update
        outlines (list): Lines from output file
    """
    zero_point_corr = 0.0  # Initialize for G4 calculations
    
    for i in reversed(range(0, len(outlines) - 30)):
        # TD-DFT calculations
        if "E(TD-HF/TD-DFT)" in outlines[i]:
            td_e = float(outlines[i].strip().split()[-1])
            if not hasattr(cclib_data, "energy"):
                cclib_data["energy"] = {}
            cclib_data["energy"]["TD energy"] = cclib.parser.utils.convertor(
                td_e, "hartree", "eV"
            )
        
        # G4 calculations
        elif outlines[i].strip().startswith("E(ZPE)="):
            zero_point_corr = float(outlines[i].strip().split()[1])
        elif outlines[i].strip().startswith("G4(0 K)"):
            G4_energy = float(outlines[i].strip().split()[2])
            G4_energy -= zero_point_corr
            if not hasattr(cclib_data, "energy"):
                cclib_data["energy"] = {}
            cclib_data["energy"]["G4 energy"] = cclib.parser.utils.convertor(
                G4_energy, "hartree", "eV"
            )
        
        # ONIOM calculations
        elif "ONIOM: extrapolated energy" in outlines[i].strip():
            oniom_e = float(outlines[i].strip().split()[4])
            if not hasattr(cclib_data, "energy"):
                cclib_data["energy"] = {}
            cclib_data["energy"]["ONIOM energy"] = cclib.parser.utils.convertor(
                oniom_e, "hartree", "eV"
            )
        
        # S**2 values
        elif "S**2 before annihilation" in outlines[i]:
            cclib_data["S2 after annihilation"] = float(
                outlines[i].strip().split()[-1]
            )
            cclib_data["S2 before annihilation"] = float(
                outlines[i].strip().split()[-3][:-1]
            )
        
        # Point group
        elif "Full point group" in outlines[i]:
            point_group = outlines[i].strip().split()[3]
            cclib_data["rotational"]["symmetry point group"] = point_group
            break
        
        # Convergence
        elif "Stationary point found" in outlines[i]:
            cclib_data["opt times converged"] = 2
        
        # Symmetry number
        elif "Rotational symmetry number" in outlines[i]:
            symmno = int(outlines[i].strip().split()[3].split(".")[0])
            cclib_data["rotational"]["symmetry number"] = symmno
        
        # Rotational constants
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
            cclib_data["rotational"]["rotational constants"] = roconst
        
        # Rotational temperatures
        elif outlines[i].find("Rotational temperature ") > -1:
            rotemp = [float(outlines[i].strip().split()[3])]
            cclib_data["rotational"]["rotational temperatures"] = rotemp
        
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
                        float(outlines[i].strip().replace(":", " ").split()[4]),
                        float(outlines[i].strip().replace(":", " ").split()[5]),
                    ]
            cclib_data["rotational"]["rotational temperatures"] = rotemp


def _extract_orca_keywords(outlines, i):
    """Extract keywords line from ORCA output.
    
    Args:
        outlines (list): Lines from output file
        i (int): Starting line index
    
    Returns:
        str: Complete keywords line
    """
    keywords_line = ""
    for j in range(i, i + 100):
        if "*" in outlines[j]:
            break
        else:
            keywords_line += outlines[j][6:]
    return keywords_line[1:].rstrip("\n")


def _parse_orca_keywords(keywords_line, cclib_data):
    """Parse ORCA keywords for calculation type and processors.
    
    Args:
        keywords_line (str): ORCA keywords line
        cclib_data (dict): cclib data dictionary to update
    """
    calc_type = "ground_state"
    for keyword in keywords_line.split():
        if keyword.lower() in ["optts", 'neb-ts']:
            calc_type = "transition_state"
            break
        if keyword.lower()[0:3] == 'pal':
            cclib_data["metadata"]["processors"] = keyword[3]
    cclib_data["metadata"]["ground or transition state"] = calc_type


def _extract_orca_metadata(cclib_data, outlines):
    """Extract ORCA-specific metadata from output file.
    
    Args:
        cclib_data (dict): cclib data dictionary to update
        outlines (list): Lines from output file
    """
    # Track convergence
    if "optdone" in cclib_data:
        if cclib_data["optdone"] in [True, 'true']:
            cclib_data["opt times converged"] = 1
    
    # Extract final energy
    for i in reversed(range(0, len(outlines))):
        if outlines[i][:25] == "FINAL SINGLE POINT ENERGY":
            orca_e = float(outlines[i].split()[-1])
            if not hasattr(cclib_data, "energy"):
                cclib_data["energy"] = {}
            cclib_data["energy"]["final single point energy"] = cclib.parser.utils.convertor(
                orca_e, "hartree", "eV"
            )
            break
    
    # Extract input parameters
    for i, line in enumerate(outlines):
        # Number of processors
        if "%pal" in line:
            pal_line = ''
            for j in range(i, i + 3):
                if outlines[j][0] not in ['%', '!'] or "%pal" in outlines[j]:
                    pal_line += outlines[j].rstrip("\n")[5:]
            if 'nprocs' in pal_line:
                nprocs = pal_line.strip().split()[2]
                cclib_data["metadata"]["processors"] = nprocs
        
        # Memory
        elif '%maxcore' in line:
            mem = int(line.strip().split()[3])
            cclib_data["metadata"]["memory"] = f'{mem}MB'
        
        # Input line
        elif "!" in line:
            keywords_line = _extract_orca_keywords(outlines, i)
            cclib_data["metadata"]["keywords line"] = keywords_line
            _parse_orca_keywords(keywords_line, cclib_data)
        
        elif 'END OF INPUT' in line:
            break


def _convert_and_save_json(cclib_data, file):
    """Convert NumPy arrays to lists and save as JSON.
    
    Args:
        cclib_data (dict): cclib data dictionary
        file (str): Output file path
    
    Returns:
        dict: Converted cclib_data
    """
    name_path = os.path.basename(Path(file))
    dir_path = os.path.dirname(Path(file))
    
    def _convert_ndarrays(obj):
        """Recursively convert all numpy arrays to lists."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: _convert_ndarrays(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_convert_ndarrays(v) for v in obj]
        elif isinstance(obj, tuple):
            return tuple(_convert_ndarrays(v) for v in obj)
        else:
            return obj
    
    cclib_data = _convert_ndarrays(cclib_data)
    
    with open(f'{dir_path}/{name_path.split(".")[0]}.json', "w") as outfile:
        json.dump(cclib_data, outfile, indent=1)
    
    return cclib_data


def get_json_data(self, file, cclib_data):
    """Extract metadata and GoodVibes data for JSON file.
    
    Parses QM output files (Gaussian, ORCA) to extract metadata,
    energy values, rotational constants, and other properties.
    
    Parameters
    ----------
    self : object
        Object with args.w_dir_main attribute
    file : str
        Path to QM output file
    cclib_data : dict
        cclib parsed data dictionary
    
    Returns
    -------
    dict
        Updated cclib_data with metadata and properties
    """
    outlines = read_file(os.getcwd(), self.args.w_dir_main, file)
    
    # Detect QM program
    metadata = _detect_qm_program(outlines)
    if metadata:
        cclib_data["metadata"] = metadata
    
    # Process based on QM program
    if cclib_data["metadata"]["QM program"].lower().find("gaussian") > -1:
        _extract_gaussian_metadata(cclib_data, outlines)
        _extract_gaussian_properties(cclib_data, outlines)
    
    elif cclib_data["metadata"]["QM program"].lower().find("orca") > -1:
        _extract_orca_metadata(cclib_data, outlines)
    
    # Save to JSON if data exists
    if cclib_data != {}:
        cclib_data = _convert_and_save_json(cclib_data, file)
    
    return cclib_data

def _extract_scf_energy(cclib_data):
    """Extract SCF energy from cclib data.
    
    Args:
        cclib_data (dict): cclib parsed data
    
    Returns:
        float: Energy in hartree
    """
    # in eV, converted to hartree using the conversion factor from cclib
    E_dup = cclib_data["scfenergies"][-1]
    E_dup = cclib.parser.utils.convertor(E_dup, "eV", "hartree")
    return E_dup


def _extract_thermodynamic_data(cclib_data, E_dup, errortype):
    """Extract enthalpy and free energy from cclib data.
    
    For single atoms or SP calculations, uses SCF energy as H and G.
    
    Args:
        cclib_data (dict): cclib parsed data
        E_dup (float): SCF energy in hartree
        errortype (str): Current error type
    
    Returns:
        tuple: (H_dup, G_dup, errortype)
    """
    try:
        H_dup = cclib_data["enthalpy"]
        G_dup = cclib_data["freeenergy"]
    except (AttributeError, KeyError):
        if cclib_data["natom"] == 1:
            if cclib_data["metadata"]["keywords line"].find("freq") == -1:
                errortype = "sp_calc"
                cclib_data["metadata"]["ground or transition state"] = "SP calculation"
            H_dup = E_dup
            G_dup = E_dup
        else:
            raise  # Re-raise for multi-atom systems
    
    return H_dup, G_dup, errortype


def _extract_rotational_constants(cclib_data):
    """Extract rotational constants from cclib data.
    
    Validates that exactly 3 rotational constants exist.
    
    Args:
        cclib_data (dict): cclib parsed data
    
    Returns:
        list or None: Rotational constants or None if invalid
    """
    try:
        ro_dup = cclib_data["rotational"]["rotational constants"]
        if len(ro_dup) != 3:
            ro_dup = None
    except (KeyError, TypeError):
        ro_dup = None
    
    return ro_dup


def get_cclib_params(cclib_data, errortype):
    """Retrieve energy and rotational constant information from cclib data.
    
    Extracts SCF energy, enthalpy, free energy, and rotational constants.
    Handles special cases like single atoms and SP calculations.
    
    Parameters
    ----------
    cclib_data : dict
        Parsed cclib data dictionary
    errortype : str
        Current error type classification
    
    Returns
    -------
    tuple
        (E_dup, H_dup, G_dup, ro_dup, errortype) where:
        - E_dup (float): SCF energy in hartree
        - H_dup (float): Enthalpy in hartree
        - G_dup (float): Free energy in hartree
        - ro_dup (list or None): Rotational constants (GHz) or None
        - errortype (str): Updated error type
    """
    # Extract SCF energy
    E_dup = _extract_scf_energy(cclib_data)
    
    # Extract thermodynamic data
    H_dup, G_dup, errortype = _extract_thermodynamic_data(cclib_data, E_dup, errortype)
    
    # Extract rotational constants
    ro_dup = _extract_rotational_constants(cclib_data)
    
    return E_dup, H_dup, G_dup, ro_dup, errortype
