"""
QDESCP (Quantum Descriptor Calculator) Utility Module

This module provides comprehensive functionality for calculating and analyzing quantum
chemical descriptors for molecular systems. It integrates various computational
chemistry tools and methods to generate meaningful molecular and atomic descriptors.

Key Features:
    - Calculation of quantum chemical descriptors using xTB methods (with MORFEUS)
    - Boltzmann averaging of molecular properties
    - RDKit-based molecular descriptors
    - Processing of NMR chemical shifts
    - MORFEUS-based geometric and electronic analysis
    - Automated pattern recognition in molecular structures
"""

from __future__ import annotations
from typing import Dict, List, Optional, Union, Any, Tuple
import json
import sys
import os
import re
import ast
import warnings
from pathlib import Path

# Scientific computing imports
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS, Descriptors

# MORFEUS package imports
from morfeus import (
    SASA, Dispersion, BuriedVolume, ConeAngle, 
    SolidAngle, Pyramidalization, read_xyz, 
    read_geometry, XTB
)
from morfeus.data import HARTREE_TO_KCAL
from aqme.utils import load_sdf, periodic_table

# Suppress warnings
warnings.filterwarnings('ignore')

# Type aliases
NDArray = np.ndarray
MoleculeType = Union[rdkit.Chem.rdchem.Mol, None]

# Physical constants and unit conversions
GAS_CONSTANT = 8.3144621  # Gas constant in J/(K⋅mol)
J_TO_AU = 4.184 * 627.509541 * 1000.0  # Joules to atomic units conversion
TEMPERATURE = 298.15  # Default temperature in Kelvin
HARTREE_TO_EV = 27.2114  # Hartree to eV conversion

def convert_ndarrays(data: Dict[str, Any]) -> None:
    """
    Recursively convert numpy arrays to lists for JSON serialization.
    
    This function modifies the input dictionary in-place, converting any numpy
    arrays to lists so they can be serialized to JSON. It recursively processes
    nested dictionaries.
    
    Args:
        data: Dictionary containing values to convert
        
    Notes:
        - Modifies the dictionary in-place
        - Handles nested dictionaries recursively
        - Only converts numpy arrays, preserves other types
    """
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            data[key] = value.tolist()
        elif isinstance(value, dict):
            convert_ndarrays(value)

def get_boltz(energy):
    """
    Calculates the Boltzmann weights for a list of energies.
    
    Parameters:
    energy (list of floats): List of energies for which Boltzmann weights will be calculated.

    Returns:
        List of normalized Boltzmann weights corresponding to the input energies.
        The weights sum to 1.0.

    Notes:
        - Energies are shifted to make the minimum energy zero to prevent
          numerical underflow in exponential calculations
        - Uses the formula: w_i = exp(-E_i/kT) / sum(exp(-E_j/kT))
    """
    if not energy:
        return []
        
    # Shift energies to prevent numerical underflow
    shifted_energies = np.array(energy) - min(energy)
    
    # Calculate Boltzmann factors
    boltz_factors = np.exp(-shifted_energies * J_TO_AU / (GAS_CONSTANT * TEMPERATURE))
    
    # Normalize to get weights
    weights = boltz_factors / np.sum(boltz_factors)
    
    return weights.tolist()

def get_boltz_props_nmr(
    json_files: List[str],
    name: str,
    boltz_dir: str,
    self,
    atom_props: List[str],
    nmr_atoms: Optional[List[int]] = None,
    nmr_slope: Optional[List[float]] = None,
    nmr_intercept: Optional[List[float]] = None,
    nmr_experim: Optional[str] = None
) -> Dict[str, Any]:
    """
    Process NMR properties and calculate Boltzmann-weighted averages.

    This function handles the calculation and analysis of NMR chemical shifts,
    including comparison with experimental data when available.

    Args:
        json_files: List of paths to JSON files containing conformer data
        name: Base name for output files
        boltz_dir: Directory for Boltzmann-weighted results
        self: Instance containing logging capabilities
        atom_props: List of atomic properties to process
        nmr_atoms: List of atom indices for NMR calculation
        nmr_slope: Slope parameters for NMR shift scaling
        nmr_intercept: Intercept parameters for NMR shift scaling
        nmr_experim: Path to experimental NMR data file (optional)

    Returns:
        Dictionary containing processed NMR data and Boltzmann averages

    Raises:
        FileNotFoundError: If experimental data file cannot be loaded
        KeyError: If required JSON fields are missing
    """
    # Initialize experimental data
    exp_data = None
    if nmr_experim:
        try:
            exp_data = pd.read_csv(nmr_experim)
        except Exception as e:
            error_msg = f'Error loading NMR shifts file: {nmr_experim}.\n{e}'
            self.args.log.write(f'\nx  {error_msg}')
            self.args.log.finalize()
            sys.exit()

    # Calculate Boltzmann weights based on energies
    energy = []
    for k, json_file in enumerate(json_files):
        json_data = read_json(json_file)
        
        # Extract energy for Boltzmann factor calculation
        try:
            energy_value = json_data["optimization"]["scf"]["scf energies"][-1]  # SCF energy
        except KeyError:
            self.args.log.write(f"x  No SCF energy found in JSON file: {json_file}")
            energy_value = 0
        energy.append(energy_value)

        # Calculate NMR chemical shifts
        json_data["properties"]["NMR"]["NMR Chemical Shifts"] = get_chemical_shifts(json_data, nmr_atoms, nmr_slope, nmr_intercept)

        # Merge with experimental data if available
        if nmr_experim is not None:
            list_shift = json_data["properties"]["NMR"]["NMR Chemical Shifts"]
            df = pd.DataFrame(
                list_shift.items(),
                columns=["atom_idx", "conf_{}".format(k + 1)],
            )
            df["atom_idx"] = df["atom_idx"] + 1  # Adjust index to match experimental data
            exp_data = exp_data.merge(df, on=["atom_idx"], how='left')

        # Save updated JSON data back to the file
        with open(json_file, "w", encoding='utf-8') as outfile:
            json.dump(json_data, outfile)

    # Calculate Boltzmann factors
    boltz = get_boltz(energy)

    # Process atomic properties (NMR Chemical Shifts)
    full_json_data = {}
    for i, prop in enumerate(atom_props):
        prop_list = []
        for json_file in json_files:
            json_data = read_json(json_file)
            prop_list.append(list(json_data["properties"]["NMR"][prop].values()))

        # Average properties based on Boltzmann weights
        avg_prop = average_prop_atom_nmr(boltz, prop_list)
        dictavgprop = {}
        for j, key in enumerate(json_data["properties"]["NMR"][prop].keys()):
            dictavgprop[key] = avg_prop[j]
        full_json_data[prop] = dictavgprop

    # Merge with experimental data if available and calculate error
    if nmr_experim is not None:
        list_shift = full_json_data[prop]
        df = pd.DataFrame(
            list_shift.items(),
            columns=["atom_idx", "boltz_avg"],
        )
        df["atom_idx"] = df["atom_idx"].astype(int) + 1
        exp_data = exp_data.merge(df, on=["atom_idx"])
        exp_data["error_boltz"] = abs(
            exp_data["experimental_ppm"] - exp_data["boltz_avg"]
        )

        # Save the experimental data with predictions
        qdescp_nmr = nmr_experim.split(".csv")[0] + "_predicted.csv"
        exp_data.round(2).to_csv(qdescp_nmr, index=False)
        self.args.log.write(f"o  The {os.path.basename(qdescp_nmr)} file containing Boltzmann weighted NMR shifts was successfully created in {self.args.initial_dir}")

    # Save averaged properties to a JSON file
    final_boltz_file = str(boltz_dir) + "/" + name + "_boltz.json"
    with open(final_boltz_file, "w", encoding='utf-8') as outfile:
        json.dump(full_json_data, outfile)

def get_chemical_shifts(
    json_data: Dict[str, Any],
    nmr_atoms: Union[List[int], str],
    nmr_slope: Union[List[float], str],
    nmr_intercept: Union[List[float], str]
) -> Dict[int, float]:
    """
    Calculate scaled NMR chemical shifts from JSON data.

    This function processes NMR tensors from quantum chemistry calculations
    and applies appropriate scaling factors to obtain chemical shifts.

    Args:
        json_data: Dictionary containing quantum chemistry results
        nmr_atoms: List of atomic numbers or string representation
        nmr_slope: List of slope parameters or string representation
        nmr_intercept: List of intercept parameters or string representation

    Returns:
        Dictionary mapping atom indices to their scaled chemical shifts

    Notes:
        - Scaling formula: shift = (intercept - tensor) / (-slope)
        - Only processes atoms specified in nmr_atoms
        - Handles both list and string inputs for parameters
    """
    # Convert string inputs to lists if necessary
    def ensure_list(param):
        return ast.literal_eval(param) if isinstance(param, str) else param

    nmr_atoms = ensure_list(nmr_atoms)
    nmr_slope = ensure_list(nmr_slope)
    nmr_intercept = ensure_list(nmr_intercept)

    # Extract required data
    atomic_numbers = json_data["atoms"]["elements"]["number"]
    nmr_tensors = json_data["properties"]["NMR"]["NMR isotopic tensors"]
    
    # Calculate shifts
    shifts = {}
    for i, (atom, tensor) in enumerate(zip(atomic_numbers, nmr_tensors)):
        if atom in nmr_atoms:
            idx = nmr_atoms.index(atom)
            slope = nmr_slope[idx]
            intercept = nmr_intercept[idx]
            
            # Apply scaling formula
            scaled_shift = (intercept - tensor) / (-slope)
            shifts[i] = scaled_shift
            
    return shifts

def average_prop_atom(
    weights: List[float],
    properties: List[Union[float, List[float], NDArray]]
) -> Union[float, NDArray]:
    """
    Calculate Boltzmann-weighted average of atomic properties.

    Handles both scalar and vector atomic properties, with robust error handling
    for missing or invalid data.

    Args:
        weights: Boltzmann weights for each configuration
        properties: List of atomic properties. Each property can be:
                   - A single number (float/int)
                   - A list of numbers (for multi-dimensional properties)
                   - A numpy array

    Returns:
        For scalar properties: float (rounded to 4 decimal places)
        For vector properties: numpy array (rounded to 4 decimal places)
        Returns np.nan if calculation fails

    Notes:
        - Handles xTB calculation failures
        - Supports both scalar and vector properties
        - Automatically handles None values in lists
        - Returns np.nan for any calculation errors
    """
    try:
        # Verify input lengths match
        if len(properties) != len(weights):
            return np.nan

        # Convert properties to weighted values
        weighted_props = []
        for prop, weight in zip(properties, weights):
            if isinstance(prop, (float, int)):
                weighted_props.append(prop * weight)
            elif isinstance(prop, (list, np.ndarray)):
                # Handle None values in lists
                weighted_values = [
                    0 if val is None else val * weight 
                    for val in prop
                ]
                weighted_props.append(weighted_values)
            else:
                weighted_props.append(np.nan)

        # Calculate final result based on property type
        if not weighted_props:
            return np.nan

        if isinstance(weighted_props[0], (float, int)):
            # Handle scalar properties
            result = sum(x for x in weighted_props if not pd.isna(x))
            return round(result, 4)
        else:
            # Handle vector properties
            try:
                result = np.sum(weighted_props, axis=0)
                return np.round(result, 4)
            except (ValueError, TypeError):
                return np.nan

    except Exception:
        return np.nan

def average_prop_atom_nmr(
    weights: List[float],
    properties: List[List[float]]
) -> Union[NDArray, float]:
    """
    Calculate Boltzmann-weighted average of atomic NMR properties.

    This function handles NMR-specific property averaging, with special handling
    for missing or invalid data.

    Args:
        weights: Boltzmann weights for each configuration
        properties: List of NMR property lists for each configuration

    Returns:
        Numpy array of averaged properties or np.nan if calculation fails

    Notes:
        - Handles NaN and None values
        - Uses vectorized operations for efficiency
        - Returns np.nan if any property list contains invalid values
    """
    try:
        # Check for invalid values
        for prop_list in properties:
            if any(p is None or pd.isna(p) for p in prop_list):
                return np.nan

        # Convert to numpy array for efficient calculation
        prop_array = np.array(properties)
        
        # Calculate weighted average along first axis (configurations)
        weighted_avg = np.sum(prop_array * np.array(weights)[:, np.newaxis], axis=0)
        
        return weighted_avg
        
    except (ValueError, TypeError):
        return np.nan

def average_prop_mol(weights, prop):
    """
    Calculate Boltzmann-weighted average of molecular properties.

    This function computes the weighted average of molecular properties using
    Boltzmann weights from different molecular conformations.

    Parameters:
        weights (list of float): Boltzmann weights for each configuration,
                                should sum to approximately 1.0
        prop (list of float): List of molecular properties for each configuration

    Returns:
    boltz_avg (float): Boltzmann-weighted average of the molecular properties, rounded to 4 decimal places.
    """
    
    # Initialize the Boltzmann-weighted average to 0.0
    boltz_avg = 0.0
    
    # Loop through each property and its corresponding weight
    for i, p in enumerate(prop):
        # If the property is 'NaN', break the loop and return 'NaN'
        if str(p).lower() == 'nan' or p is None:
            boltz_avg = np.nan
            break
        # Otherwise, calculate the weighted sum of properties
        boltz_avg += p * weights[i]
    
    # If the result is a valid number, round it to 4 decimal places
    if str(boltz_avg).lower() != 'nan':
        boltz_avg = round(boltz_avg, 4)

    # Return the Boltzmann-weighted average, rounded to 4 decimal places (or 'NaN' if encountered)
    return boltz_avg

def get_rdkit_properties(
    self,
    full_json_data: Dict[str, Any],
    mol: MoleculeType
) -> Dict[str, Any]:
    """
    Calculate and store molecular descriptors using RDKit.

    Computes a comprehensive set of molecular descriptors available in RDKit,
    handling both modern and legacy RDKit versions.

    Args:
        self: Instance containing logging capabilities
        full_json_data: Dictionary to store calculated descriptors
        mol: RDKit molecule object for analysis

    Returns:
        Updated dictionary containing valid RDKit descriptors

    Notes:
        - Filters out NaN values automatically
        - Falls back to basic descriptors for older RDKit versions
        - Warns if using legacy version with limited descriptors
    """
    if mol is None:
        return full_json_data

    try:
        # Calculate all available descriptors
        descriptors = Descriptors.CalcMolDescriptors(mol)
        
        # Filter and store valid descriptors
        for name, value in descriptors.items():
            if value is not None and not pd.isna(value):
                full_json_data[name] = value

    except AttributeError:
        # Fallback for older RDKit versions
        warning_msg = (
            "WARNING! Install a newer version of RDKit to get all descriptors. "
            "You can use: 'pip install rdkit --upgrade'."
        )
        self.args.log.write(f"x  {warning_msg}")
        
        # Calculate basic descriptor (MolLogP)
        try:
            full_json_data["MolLogP"] = Descriptors.MolLogP(mol)
        except Exception as e:
            self.args.log.write(f"x  Error calculating MolLogP: {e}")

    try:
        descrs = rdkit.Chem.Descriptors.CalcMolDescriptors(mol)
        for descr in descrs:
            if str(descrs[descr]).lower() != 'nan':
                full_json_data[descr] = descrs[descr]
    except AttributeError:
        self.args.log.write(f"x  WARNING! Install a newer version of RDKit to get all the RDKit descriptors in the database with all the descriptors. You can use: 'pip install rdkit --upgrade'.")
        full_json_data["MolLogP"] = rdkit.Chem.Descriptors.MolLogP(mol)

    return full_json_data


def read_json(file: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """
    Load and parse a JSON file with robust error handling.

    This function attempts to read a JSON file, with special handling for
    Windows-specific backslash escaping issues. It includes multiple fallback
    strategies for parsing problematic JSON files.

    Args:
        file: Path to the JSON file to read

    Returns:
        Dictionary containing the parsed JSON data, or None if parsing fails

    Notes:
        - Attempts standard JSON parsing first
        - On failure, tries to fix common Windows path escaping issues
        - Returns None for non-JSON files or if all parsing attempts fail
        - Uses UTF-8 encoding for universal character support
    """
    if not str(file).endswith(".json"):
        return None

    try:
        # First attempt: Standard JSON parsing
        with open(file, "r", encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError:
        try:
            # Second attempt: Fix Windows backslash issues
            with open(file, "r", encoding='utf-8') as f:
                content = f.read()
                # Fix unescaped backslashes not part of known escape sequences
                fixed_content = re.sub(r'\\(?![\\nt"rbfu/])', r'\\\\', content)
                return json.loads(fixed_content)
        except Exception:
            return None
    except Exception:
        return None


def get_matches_idx_n_prefix(
    self,
    smarts_targets: List[str],
    name_initial: str
) -> Dict[str, Dict[str, Union[List[int], List[str], int]]]:
    """
    Locate atom indices and generate atom prefixes for given SMARTS patterns.

    This function identifies specific atoms in a molecule that match given SMARTS
    patterns and generates descriptive prefixes for them. It's used for analyzing
    specific functional groups or atomic patterns in molecules.

    Args:
        self: Instance containing logging capabilities
        smarts_targets: List of SMARTS patterns to search for
        name_initial: Path to the initial molecular structure file

    Returns:
        Dictionary mapping each SMARTS pattern to its analysis results:
        - sorted_indices: List of matched atom indices
        - match_names: List of descriptive atom prefixes
        - n_types: Number of unique atom types in match

    Notes:
        - Handles both CSEARCH-generated and regular SDF files
        - Sanitizes patterns with quotes
        - Warns about missing or multiple matches
        - Creates unique, consistent atom naming across molecules
    """
    pattern_dict: Dict[str, Dict[str, Union[List[int], List[str], int]]] = {}

    if len(smarts_targets) > 0:
        # Create RDKit mol object from input file
        mol = get_mol_assign(self,name_initial)

        # Process each SMARTS pattern
        for pattern in smarts_targets:
            # Sanitize pattern by removing quotes
            if "'" in pattern or '"' in pattern:
                pattern = pattern.replace("'", "").replace('"', "")

            # Find atoms matching the pattern
            matches, idx_set = get_atom_matches(self, pattern, mol)

            # Handle different match scenarios
            if len(matches) == 0:
                self.args.log.write(
                    f"x  WARNING! SMARTS pattern {pattern} not found in "
                    f"{os.path.basename(name_initial)}, atomic descriptors will "
                    "not be generated."
                )
            elif matches[0] == -1:
                self.args.log.write(
                    f"x  WARNING! Mapped atom {pattern} not found in "
                    f"{os.path.basename(name_initial)}, atomic descriptors will "
                    "not be generated."
                )
            elif len(matches) > 1:
                self.args.log.write(
                    f"x  WARNING! More than one {pattern} patterns were found in "
                    f"{os.path.basename(name_initial)}, atomic descriptors will "
                    "not be generated."
                )
            elif len(matches) == 1:
                # Get sorted atom types for consistent ordering
                n_types, sorted_indices = sort_atom_types(matches, mol)

                # Generate descriptive atom names
                match_names = get_prefix_atom_props(
                    sorted_indices,
                    mol,
                    pattern,
                    idx_set
                )

                # Store pattern analysis results
                pattern_dict[pattern] = {
                    'sorted_indices': sorted_indices,
                    'match_names': match_names,
                    'n_types': n_types,
                }

    return pattern_dict

def setup_env(self):
    """Configure environment variables for reproducible calculations.
    
    Returns:
        dict: Environment variables
    """
    env = dict(os.environ)
    env.update({
        "OMP_STACKSIZE": self.args.stacksize,
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", 
        "OPENBLAS_NUM_THREADS": "1",
        'OMP_MAX_ACTIVE_LEVELS': '1',
        "NUMEXPR_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "MKL_DYNAMIC": "FALSE",
        "BLIS_NUM_THREADS": "1"
    })
    return env

def calculate_morfeus_descriptors(
    final_xyz_path: str,
    self,
    charge: int,
    mult: Optional[int],
    smarts_targets: List[str],
    name_initial: str
) -> Dict[str, Any]:
    """
    Calculate comprehensive molecular descriptors using the MORFEUS package.

    This function computes a wide range of molecular and atomic descriptors including:
    - Global geometric properties (SASA, dispersion)
    - Local atomic properties
    - Electronic structure properties via xTB
    - Solvation effects
    - Excited state properties

    Args:
        final_xyz_path: Path to XYZ coordinate file
        self: Instance containing logging capabilities
        charge: Molecular charge
        mult: Spin multiplicity (None for automatic)
        smarts_targets: SMARTS patterns for specific atom selection
        name_initial: Base name for file operations

    Returns:
        Dictionary containing all calculated descriptors

    Notes:
        - Uses different xTB methods (PTB, GFN) for specific properties
        - Includes error handling for each calculation type
        - Calculates excited states only for closed-shell systems
    """
    morfeus_data: Dict[str, Any] = {}

    # Log calculation start
    mol_name = os.path.basename(final_xyz_path).replace('.xyz','')
    self.args.log.write(f"\no  Running MORFEUS for {mol_name} "
                       f"with charge {charge} and multiplicity {mult}")

    # Calculate number of unpaired electrons
    n_unpaired = int(mult) - 1 if mult is not None else 0

    # Load molecular structure
    try:
        elements, coordinates = read_xyz(final_xyz_path)
        elements, coordinates = read_geometry(final_xyz_path)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error loading molecule from "
                           f"file {final_xyz_path}: {e}")
        return morfeus_data

    # Set up environment
    env = setup_env(self)

    # calculate MORFEUS global descriptors that do not come from xTB
    morfeus_data = morfeus_global_descps(self, elements, coordinates, morfeus_data)

    # calculate MORFEUS local descriptors that do not come from xTB
    morfeus_data = morfeus_local_descps(self, elements, coordinates, morfeus_data, smarts_targets, name_initial)

    # xTB calculations through MORFEUS (with 1 proc to be reproducible)
    # calculate PTB method for descriptors that support it
    morfeus_data = morfeus_ptb_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, env)
    
    # calculate GFN2 method for descriptors not included in PTB
    morfeus_data,energy = morfeus_gfn2_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, env)

    # calculate GFN2 method with ALPB H2O for descriptors related to solvation
    morfeus_data = morfeus_solv_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, env)

    # calculate GFN2 method in triplet state to calculate S0 to T1 energy gaps
    if n_unpaired == 0:
        morfeus_data = morfeus_t1_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, energy, env)

    # for atomic properties that fail
    for prop in morfeus_data:
        if isinstance(morfeus_data[prop],list):
            if len(morfeus_data[prop]) == 0:
                morfeus_data[prop] = [None]*len(morfeus_data["Buried volume"])

    return morfeus_data


def morfeus_global_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Calculate global molecular descriptors using MORFEUS methods.

    This function computes molecular-level geometric descriptors including
    solvent-accessible surface area (SASA) and dispersion properties.
    These calculations do not rely on xTB methods.

    Args:
        self: Instance containing logging capabilities
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results

    Returns:
        Updated morfeus_data dictionary with global descriptors added

    Notes:
        - Calculates total SASA for the molecule
        - Computes dispersion area and volume
        - All values rounded to 4 decimal places
        - Returns None for failed calculations
    """
    # Initialize descriptor values
    sasa_area_global = disp_area_global = disp_vol_global = None

    # Calculate Global SASA
    try:
        sasa = SASA(elements, coordinates)
        sasa_area_global = round(sasa.area, 4)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating SASA from Morfeus: {e}"
        )

    # Calculate Global Dispersion properties
    try:
        disp = Dispersion(elements, coordinates)
        disp_area_global = round(disp.area, 4)
        disp_vol_global = round(disp.volume, 4)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating Global Dispersion from Morfeus: {e}"
        )

    # Update dictionary with calculated descriptors
    morfeus_data.update({
        "SASA": sasa_area_global,
        "Dispersion area": disp_area_global,
        "Dispersion volume": disp_vol_global,
    })

    return morfeus_data


def morfeus_local_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any],
    smarts_targets: List[str],
    name_initial: str
) -> Dict[str, Any]:
    """
    Calculate atom-specific descriptors using MORFEUS methods.

    This function computes various geometric properties for specified atoms in
    the molecule. For each property, it generates a list with values for all
    atoms, using NaN for atoms not in the target set.

    Args:
        self: Instance containing logging capabilities and settings
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results
        smarts_targets: SMARTS patterns identifying atoms of interest
        name_initial: Path to initial molecular structure file

    Returns:
        Updated morfeus_data dictionary with local descriptors added

    Notes:
        - Properties calculated only for atoms matching SMARTS patterns
        - Lists have same length as number of atoms, with NaN for non-targets
        - Example: For first atom only, SASA = [2.3452, NaN, NaN, ...]
        - All values rounded to 4 decimal places (2 for buried volume)
    """
    # Initialize all property lists
    local_sasa_areas = local_buried_volumes = local_cone_angles = []
    local_solid_angles = local_Pyramidalization = []
    local_vol_Pyramidalization = local_disp = []

    # Get target atom indices from SMARTS patterns
    atom_matches = []
    pattern_dict = get_matches_idx_n_prefix(self, smarts_targets, name_initial)
    if pattern_dict:
        for pattern in pattern_dict:
            atom_matches.extend(pattern_dict[pattern]['sorted_indices'])

    # Calculate local SASA
    try:
        sasa = SASA(elements, coordinates)
        local_sasa_areas = [
            round(area, 4) for area in sasa.atom_areas.values()
        ]
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating local SASA from Morfeus: {e}"
        )

    # Calculate local dispersion
    try:
        disp = Dispersion(elements, coordinates)
        local_disp = [
            round(p_int, 4) for p_int in disp.atom_p_int.values()
        ]
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating Local Dispersion from Morfeus: {e}"
        )

    # Calculate buried volume for each atom
    try:
        local_buried_volumes = []
        for i in range(1, len(coordinates) + 1):
            if i - 1 in atom_matches:  # MORFEUS uses 1-based indices
                # Create BuriedVolume instance with appropriate radius
                if self.args.vbur_radius == 'auto':
                    bv = BuriedVolume(
                        elements,
                        coordinates,
                        i
                    )
                else:
                    bv = BuriedVolume(
                        elements,
                        coordinates,
                        i,
                        radius=self.args.vbur_radius
                    )
                buried_volume = round(bv.fraction_buried_volume * 100, 2)
                local_buried_volumes.append(buried_volume)
            else:
                local_buried_volumes.append(np.nan)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating local Buried Volume from Morfeus: {e}"
        )

    # Calculate cone angles
    try:
        local_cone_angles = []
        for i in range(1, len(coordinates) + 1):
            if i - 1 in atom_matches:
                try:
                    cone_angle = ConeAngle(elements, coordinates, i)
                    local_cone_angles.append(round(cone_angle.cone_angle, 4))
                except Exception:
                    local_cone_angles.append(np.nan)
            else:
                local_cone_angles.append(np.nan)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating local Cone Angle from Morfeus: {e}"
        )

    # Calculate solid angles
    try:
        local_solid_angles = []
        for i in range(1, len(coordinates) + 1):
            if i - 1 in atom_matches:
                try:
                    solid_angle = SolidAngle(elements, coordinates, i)
                    local_solid_angles.append(round(solid_angle.cone_angle, 4))
                except Exception:
                    local_solid_angles.append(np.nan)
            else:
                local_solid_angles.append(np.nan)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating local Solid Angle from Morfeus: {e}"
        )

    # Calculate pyramidalization
    try:
        local_Pyramidalization = []
        local_vol_Pyramidalization = []
        for i in range(1, len(coordinates) + 1):
            if i - 1 in atom_matches:
                pyr = Pyramidalization(coordinates, i)
                local_Pyramidalization.append(round(pyr.P, 4))
                local_vol_Pyramidalization.append(round(pyr.P_angle, 4))
            else:
                local_Pyramidalization.append(np.nan)
                local_vol_Pyramidalization.append(np.nan)
    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating Pyramidalization from Morfeus: {e}"
        )

    # Update morfeus_data with all local descriptors
    morfeus_data.update({
        "Atom SASA": local_sasa_areas,
        "Atom dispersion": local_disp,
        "Buried volume": local_buried_volumes,
        "Cone angle": local_cone_angles,
        "Solid angle": local_solid_angles,
        "Pyramidalization": local_Pyramidalization,
        "Pyramidaliz. volume": local_vol_Pyramidalization,
    })

    return morfeus_data


def morfeus_ptb_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any],
    charge: int,
    n_unpaired: int,
    env: Dict[str, str]
) -> Dict[str, Any]:
    """
    Calculate quantum chemical descriptors using PTB method.

    This function computes electronic structure properties using the
    PTB method implemented in xTB. It focuses on properties well-suited
    for this level of theory, such as orbital energies and dipole moments.

    Args:
        self: Instance containing logging capabilities
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons

    Returns:
        Updated morfeus_data dictionary with PTB-based descriptors

    Notes:
        - Uses PTB method from xTB for fast calculations
        - Orbital energies reported in eV
        - Dipole moments in Debye
        - Computes both molecular and atomic properties
        - Returns None for failed calculations
    """
    # Initialize electronic properties
    homo = lumo = homo_lumo_gap = dipole_moment = None
    charges = atom_dipole_vect = bond_orders = {}

    try:
        # Initialize PTB calculation
        xtb_ptb = XTB(
            elements,
            coordinates,
            env_variables=env,
            charge=charge,
            n_unpaired=n_unpaired,
            solvent=None,
            method='ptb'
        )

        # Calculate orbital properties
        homo = xtb_ptb.get_homo(unit="eV")
        lumo = xtb_ptb.get_lumo(unit="eV")
        homo_lumo_gap = xtb_ptb.get_homo_lumo_gap(unit="eV")

        # Calculate electric properties
        dipole_moment = xtb_ptb.get_dipole_moment(unit="debye")
        atom_dipole_vect = xtb_ptb.get_atom_dipoles()
        atom_dipole_modules = [
            np.linalg.norm(atom_dipole_vect[i])
            for i in atom_dipole_vect.keys()
        ]

        # Calculate atomic and bonding properties
        charges = list(xtb_ptb.get_charges(model="Mulliken").values())
        bond_orders = list(xtb_ptb.get_bond_orders().values())

    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating extra Morfeus descriptors with PTB: {e}"
        )

    # Update dictionary with calculated properties
    morfeus_data.update({
        "HOMO": homo,
        "LUMO": lumo,
        "HOMO-LUMO gap": homo_lumo_gap,
        "Dipole module": dipole_moment,
        "Dipole moment": atom_dipole_modules,
        "Partial charge": charges,
        "Bond orders": bond_orders,
    })

    return morfeus_data

def morfeus_gfn2_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any],
    charge: int,
    n_unpaired: int,
    env: Dict[str, str]
) -> Tuple[Dict[str, Any], float]:
    """
    Calculate electronic and structural descriptors using GFN-xTB method.

    This function computes various molecular and atomic descriptors not available
    in the PTB method, including electronic properties, reactivity indices,
    and polarizability.

    Args:
        self: Instance containing logging capabilities and settings
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons

    Returns:
        Tuple containing:
            - Updated morfeus_data dictionary with GFN2 descriptors
            - Total energy of the system

    Notes:
        - Uses GFN2-xTB method with optional solvent
        - Includes various reactivity indices (Fukui, HOMO/LUMO, etc.)
        - Calculates both global and local (atomic) properties
        - References for methods:
            - Normalized electrophilicity: DOI: 10.1021/jp034707u
            - Local electrophilicity: DOI: 10.1002/qua.25706
    """
    # Initialize all properties as None
    ip = ea = hardness = softness = chemical_potential = None
    polarizability = atom_polarizabilities = fod_pop = None
    fod = fukui_plus = fukui_minus = fukui_radical = fukui_dual = None
    electrophilicity = nucleofugality = None
    local_electrophil = normal_electrophil = normal_nucleophil = None
    electrofugality = energy = fermi_level = covCN = None

    try:
        # Initialize GFN2-xTB calculation
        xtb2 = XTB(
            elements,
            coordinates,
            env_variables=env,
            charge=charge,
            n_unpaired=n_unpaired,
            solvent=getattr(self.args, "qdescp_solvent", None),
            method=2
        )

        # Calculate global electronic properties
        ip = xtb2.get_ip(corrected=True)
        ea = xtb2.get_ea(corrected=True)
        energy = xtb2.get_energy()
        fermi_level = xtb2.get_fermi_level()

        # Calculate global reactivity descriptors
        electrophilicity = xtb2.get_global_descriptor("electrophilicity", corrected=True)
        nucleofugality = xtb2.get_global_descriptor("nucleofugality", corrected=True)
        electrofugality = xtb2.get_global_descriptor("electrofugality", corrected=True)
        hardness = xtb2.get_hardness()
        softness = xtb2.get_softness()
        chemical_potential = xtb2.get_chemical_potential()

        # Calculate polarizability and FOD properties
        polarizability = xtb2.get_molecular_polarizability()
        atom_polarizabilities = list(xtb2.get_atom_polarizabilities().values())
        fod_pop = list(xtb2.get_fod_population().values())
        fod = xtb2.get_nfod()

        # Calculate Fukui indices and coordination numbers
        fukui_plus = list(xtb2.get_fukui("plus").values())
        fukui_minus = list(xtb2.get_fukui("minus").values())
        fukui_radical = list(xtb2.get_fukui("radical").values())
        fukui_dual = list(xtb2.get_fukui("dual").values())
        covCN = list(xtb2.get_covcn().values())

        # Calculate local reactivity descriptors
        normal_electrophil = list(electrophilicity * np.array(fukui_plus))
        normal_nucleophil = list(-ip * np.array(fukui_minus))
        local_electrophil = list(xtb2.get_fukui("local_electrophilicity").values())

    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating extra Morfeus descriptors with GFN2: {e}"
        )

    # Update morfeus_data with all calculated properties
    morfeus_data.update({
        "total energy": energy,
        "Fermi-level": fermi_level,
        "Coord. numbers": covCN,
        "IP": ip,
        "EA": ea,
        "Hardness": hardness,
        "Softness": softness,
        "Chem. potential": chemical_potential,
        "Polarizability": polarizability,
        "Atom Polarizability": atom_polarizabilities,
        "Atom FOD": fod_pop,
        "Total FOD": fod,
        "Fukui+": fukui_plus,
        "Fukui-": fukui_minus,
        "Fukui_rad": fukui_radical,
        "Fukui dual": fukui_dual,
        "Electrophil.": local_electrophil,
        "Normaliz. electrophil.": normal_electrophil,
        "Normaliz. nucleophil.": normal_nucleophil,
        "Electrophilicity": electrophilicity,
        "Nucleofugality": nucleofugality,
        "Electrofugality": electrofugality,
    })

    return morfeus_data, energy


def morfeus_solv_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any],
    charge: int,
    n_unpaired: int,
    env: Dict[str, str]
) -> Dict[str, Any]:
    """
    Calculate solvation descriptors using GFN-xTB with ALPB water model.

    This function computes various solvation-related properties including
    solvation free energy and hydrogen bonding contributions.

    Args:
        self: Instance containing logging capabilities
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons

    Returns:
        Updated morfeus_data dictionary with solvation descriptors

    Notes:
        - Uses ALPB implicit solvation model for water
        - All energies converted to kcal/mol
        - Includes both global and atom-specific properties
        - Returns None for failed calculations
    """
    # Initialize solvation properties
    g_solv = g_solv_hb = atom_hb_terms = None

    try:
        # Setup xTB calculation with water solvation
        xtb2_solv = XTB(
            elements,
            coordinates,
            env_variables=env,
            charge=charge,
            n_unpaired=n_unpaired,
            solvent="h2o",
            method=2
        )
        
        # Calculate solvation energies and convert to kcal/mol
        g_solv = round(xtb2_solv.get_solvation_energy() * HARTREE_TO_KCAL, 2)
        g_solv_hb = round(xtb2_solv.get_solvation_h_bond_correction() * HARTREE_TO_KCAL, 2)
        
        # Calculate atomic H-bond contributions
        atom_hb_terms = [
            round(v * HARTREE_TO_KCAL, 2)
            for v in xtb2_solv.get_atomic_h_bond_corrections().values()
        ]

    except Exception as e:
        self.args.log.write(
            f"x  WARNING! Error calculating solvation descriptors in MORFEUS with GFN2: {e}"
        )

    # Update morfeus_data with solvation properties
    morfeus_data.update({
        "G solv. in H2O": g_solv,
        "G of H-bonds H2O": g_solv_hb,
        "H bond H2O": atom_hb_terms,
    })

    return morfeus_data


def morfeus_t1_descps(
    self,
    elements: List[str],
    coordinates: NDArray,
    morfeus_data: Dict[str, Any],
    charge: int,
    n_unpaired: int,
    energy: float,
    env: Dict[str, str]
) -> Dict[str, Any]:
    """
    Calculate singlet-triplet energy gap using GFN-xTB.

    This function computes the energy difference between the singlet ground
    state and first triplet excited state using GFN2-xTB method.

    Args:
        self: Instance containing logging capabilities
        elements: List of atomic symbols
        coordinates: Array of atomic coordinates
        morfeus_data: Dictionary to store calculation results
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons in ground state
        energy: Ground state energy (S0)

    Returns:
        Updated morfeus_data dictionary with S0-T1 gap added

    Notes:
        - Uses GFN2-xTB method for energy calculations
        - Energy gap is reported in kcal/mol
        - Triplet state uses n_unpaired + 2 electrons
        - Includes solvent effects if specified in args
    """
    try:
        # Calculate triplet state energy
        xtb2t1 = XTB(
            elements, 
            coordinates, 
            env_variables=env, 
            charge=charge, 
            n_unpaired=n_unpaired+2,
            solvent=getattr(self.args, "qdescp_solvent", None),
            method=2
        )
        energy_T1 = xtb2t1.get_energy()
        
        # Calculate and store S0-T1 gap in kcal/mol
        S0_T1_gap = round((energy_T1 - energy) * HARTREE_TO_KCAL, 2)
        morfeus_data["S0-T1 gap"] = S0_T1_gap
        
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating S0-T1 gap: {e}")
        morfeus_data["S0-T1 gap"] = None
        
    return morfeus_data


def full_level_boltz(
    descp_dict: Dict[str, List[str]],
    json_files: List[str],
    energy: List[float],
    smarts_targets: List[str],
    full_json_data: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Calculate Boltzmann-weighted averages for all properties at the full level.

    This function processes both atomic and molecular properties, computing
    weighted averages based on Boltzmann factors from energies.

    Args:
        descp_dict: Dictionary containing property lists by category
        json_files: List of paths to JSON files with property data
        energy: List of energies for Boltzmann weighting
        smarts_targets: List of SMARTS patterns defining atom selections
        full_json_data: Dictionary to store averaged properties

    Returns:
        Updated dictionary containing all Boltzmann-weighted properties

    Notes:
        - Handles both atomic and molecular properties
        - Checks property availability across all molecules
        - Uses Boltzmann weighting for averages
        - Gracefully handles missing properties
    """
    # Calculate Boltzmann weights once
    boltz = get_boltz(energy)

    # Process atomic properties
    atomic_props = False
    for i, prop in enumerate(descp_dict['atom_props']):
        avg_prop = np.nan

        # Check for atomic properties in first iteration
        if i == 0:
            for json_file in json_files:
                json_data = read_json(json_file)
                if any(
                    atom_prop in json_data 
                    for atom_prop in descp_dict['atom_props']
                ):
                    atomic_props = True
                    break

        # Calculate atomic property averages if available
        if atomic_props:
            try:
                prop_list = [
                    read_json(json_file)[prop]
                    for json_file in json_files
                ]
                avg_prop = average_properties(boltz, prop_list)
            except KeyError:
                full_json_data[prop] = np.nan
        
        # Update results with averaged property
        full_json_data = update_full_json_data(
            full_json_data,
            prop,
            avg_prop,
            smarts_targets
        )

    # Process molecular properties
    for prop in descp_dict['mol_props']:
        try:
            # Calculate molecular property averages
            prop_list = [
                read_json(json_file)[prop]
                for json_file in json_files
            ]
            avg_prop = average_properties(
                boltz,
                prop_list,
                is_atom_prop=False
            )
            full_json_data[prop] = avg_prop
        except KeyError:
            full_json_data[prop] = np.nan
    
    return full_json_data


def get_descriptors(level: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Get quantum chemical descriptors categorized by analysis level.

    This function provides predefined sets of molecular and atomic descriptors
    organized by analysis level (de novo, interpret, or full).

    Args:
        level: Analysis level ('denovo', 'interpret', or 'full')

    Returns:
        Nested dictionary containing descriptor lists:
        - Top level: Analysis level
        - Second level: 'mol' and 'atoms' categories
        - Values: Lists of descriptor names

    Notes:
        - De novo: Basic descriptors for initial analysis
        - Interpret: Advanced descriptors for detailed analysis
        - Full: Placeholder for custom descriptor sets
        - Returns empty dict for unknown levels
        
    Example:
        >>> descriptors = get_descriptors('denovo')
        >>> descriptors['mol']  # Get molecular descriptors
        >>> descriptors['atoms']  # Get atomic descriptors
    """
    descriptors = {
        'denovo': {
            'mol': [
                "HOMO-LUMO gap", "HOMO", "LUMO", "IP", "EA",
                "Dipole module", "SASA", "G solv. in H2O",
                "G of H-bonds H2O"
            ],
            'atoms': [
                "Partial charge", "Atom SASA", "Buried volume",
                "H bond H2O", "Fukui+", "Fukui-",
                "Electrophil.", "Normaliz. electrophil.",
                "Normaliz. nucleophil."
            ]
        },
        'interpret': {
            'mol': [
                "Fermi-level", "Polarizability", "Total FOD",
                "Hardness", "Softness", "Dispersion area",
                "Dispersion volume", "Chem. potential",
                "Electrophilicity", "Electrofugality",
                "Nucleofugality", "S0-T1 gap"
            ],
            'atoms': [
                "Fukui_rad", "Fukui dual", "Dipole moment",
                "Coord. numbers", "Atom Polarizability",
                "Atom FOD", "Atom dispersion",
                "Pyramidalization", "Pyramidaliz. volume"
            ]
        },
        'full': {
            'mol': [],
            'atoms': []
        }
    }

    return descriptors.get(level, {})


def find_level_names(
    df_level: pd.DataFrame, 
    level: str
) -> List[str]:
    """
    Select descriptors for different analysis levels.

    This function filters DataFrame columns based on the specified analysis level
    (de novo or interpret), keeping essential descriptors and level-specific ones.

    Args:
        df_level: DataFrame containing molecular descriptors
        level: Analysis level ('denovo' or 'interpret')

    Returns:
        List of column names to keep for the specified analysis level

    Notes:
        - Always includes 'code_name' and 'SMILES' if present
        - De novo level includes basic descriptors
        - Interpret level includes both basic and advanced descriptors
        - Matches column suffixes regardless of prefix
    """
    # Essential descriptors always included
    descriptors_denovo = ['code_name', 'SMILES']
    
    # Get descriptors based on analysis level
    if level == 'denovo':
        descriptors = get_descriptors('denovo')
        atom_suffixes = descriptors['atoms'] + descriptors['mol']
    if level == 'interpret':
        denovo_descr = get_descriptors('denovo')
        interpret_descr = get_descriptors('interpret')
        atom_suffixes = (
            denovo_descr['atoms'] + interpret_descr['atoms'] +
            denovo_descr['mol'] + interpret_descr['mol']
        )

    # Build list of columns to keep
    cols_to_keep = [
        c for c in descriptors_denovo 
        if c in df_level.columns
    ]
    
    # Add columns matching descriptor suffixes
    cols_to_keep.extend(
        col for col in df_level.columns
        if any(col.endswith(suffix) for suffix in atom_suffixes)
    )

    return cols_to_keep


def fix_cols_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize column names in descriptor DataFrame.

    Ensures consistent capitalization and naming conventions for key columns
    including SMILES, code_name, charge, and multiplicity.

    Args:
        df: Input DataFrame with molecular descriptors

    Returns:
        DataFrame with standardized column names

    Notes:
        - Case-insensitive matching for column names
        - Standardizes: SMILES, code_name, charge, mult
        - Preserves all other column names as-is
    """
    name_mapping = {
        'smiles': 'SMILES',
        'code_name': 'code_name',
        'charge': 'charge',
        'mult': 'mult'
    }
    
    # Create mapping of actual column names to standardized names
    rename_dict = {}
    for col in df.columns:
        lower_col = col.lower()
        if lower_col in name_mapping:
            rename_dict[col] = name_mapping[lower_col]
            
    # Apply renaming if any matches found
    if rename_dict:
        df = df.rename(columns=rename_dict)
        
    return df


def remove_atom_descp(
    df: pd.DataFrame,
    atom_props: List[str]
) -> pd.DataFrame:
    """
    Remove atomic descriptors from a dataframe that weren't explicitly specified.

    Args:
        df: Pandas DataFrame containing molecular descriptors
        atom_props: List of atomic property names to keep

    Returns:
        DataFrame with unspecified atomic descriptors removed

    Notes:
        - Modifies column structure but preserves row indexing
        - Resets index after dropping columns
        - Returns original DataFrame if no columns match atom_props
    """
    cols_to_drop = [col for col in df.columns if col in atom_props]
    if cols_to_drop:
        return df.drop(cols_to_drop, axis=1).reset_index(drop=True)
    return df

def assign_prefix_atom_props(
    prefix_list: List[str],
    atom_props: List[str],
    interpret_atoms: List[str],
    denovo_atoms: List[str]
) -> Tuple[List[str], List[str], List[str]]:
    """
    Assign prefixes to different sets of atomic properties.

    This function applies a set of prefixes to three different lists of
    atomic properties, maintaining consistency across property types.

    Args:
        prefix_list: List of prefixes to apply
        atom_props: List of general atomic properties
        interpret_atoms: List of interpretation-level atomic properties
        denovo_atoms: List of de novo-level atomic properties

    Returns:
        Tuple containing:
        - Modified atom_props with prefixes
    """
    return (
        add_prefix(prefix_list, atom_props),
        add_prefix(prefix_list, interpret_atoms),
        add_prefix(prefix_list, denovo_atoms)
    )


def add_prefix(
    prefix_list: List[str], 
    prop_list: List[str]
) -> List[str]:
    """
    Add prefixes to atomic property names.

    Creates a new list of property names by combining each prefix with
    each property name. Useful for distinguishing properties of different
    atom types or positions.

    Args:
        prefix_list: List of prefixes to add
        prop_list: List of property names to modify

    Returns:
        List of property names with prefixes added
    """
    # Create list of prefixed property names for each prefix
    prefixed_props = [
        [f'{prefix}{prop}' for prop in prop_list]
        for prefix in prefix_list
    ]
    
    # Flatten list of lists into single list
    return [
        prop for prefix_group in prefixed_props 
        for prop in prefix_group
    ]


def collect_descp_lists() -> Dict[str, Union[List[str], str]]:
    """
    Organize descriptor names for different analysis levels.

    This function collects and organizes molecular and atomic descriptors
    into appropriate categories for different levels of analysis (de novo,
    interpret, and full).

    Returns:
        Dictionary containing:
        - Lists of molecular descriptors by level
        - Lists of atomic descriptors by level
        - CSV filename templates for each level
        
    Example structure:
        {
            'denovo_mols': ['HOMO', 'LUMO', ...],
            'denovo_atoms': ['charge', 'SASA', ...],
            'interpret_mols': ['polarizability', ...],
            'interpret_atoms': ['fukui', ...],
            'mol_props': [...],
            'atom_props': [...],
            'qdescp_denovo_csv': "filename.csv",
            ...
        }
    """
    # Collect descriptors for each level
    levels = {
        'denovo': get_descriptors('denovo'),
        'interpret': get_descriptors('interpret'),
        'full': get_descriptors('full')
    }
    
    # Organize molecular descriptors
    denovo_mols = levels['denovo']['mol']
    interpret_mols = denovo_mols + levels['interpret']['mol']
    mol_props = interpret_mols + levels['full']['mol']
    
    # Organize atomic descriptors
    denovo_atoms = levels['denovo']['atoms']
    interpret_atoms = denovo_atoms + levels['interpret']['atoms']
    atom_props = interpret_atoms + levels['full']['atoms']
    
    # Return organized dictionary
    return {
        # Molecular descriptors by level
        'denovo_mols': denovo_mols,
        'interpret_mols': interpret_mols,
        'mol_props': mol_props,
        
        # Atomic descriptors by level
        'denovo_atoms': denovo_atoms,
        'interpret_atoms': interpret_atoms,
        'atom_props': atom_props,
        
        # Output file templates
        'qdescp_denovo_csv': "QDESCP_denovo_descriptors.csv",
        'qdescp_interpret_csv': "QDESCP_interpret_descriptors.csv",
        'qdescp_csv': "QDESCP_full_descriptors.csv"
    }

def average_properties(
    boltz_weights: List[float],
    prop_list: List[Union[float, List[float], NDArray]],
    is_atom_prop: bool = True
) -> Union[float, NDArray]:
    """
    Calculate Boltzmann-weighted averages for molecular or atomic properties.

    This function serves as a unified interface for averaging properties,
    delegating to specialized functions based on property type.

    Args:
        boltz_weights: List of Boltzmann weights for each configuration
        prop_list: List of property values (scalar or array) for each configuration
        is_atom_prop: True for atomic properties, False for molecular properties

    Returns:
        Boltzmann-weighted average of the properties.
        Returns scalar for molecular properties, array for atomic properties.
    """
    if is_atom_prop:
        return average_prop_atom(boltz_weights, prop_list)
    return average_prop_mol(boltz_weights, prop_list)


def update_full_json_data(
    full_json_data: Dict[str, Any],
    prop: str,
    avg_prop: Union[float, NDArray],
    smarts_targets: List[str]
) -> None:
    """
    Update JSON data structure with averaged property values.

    This function handles the insertion of averaged molecular properties into
    the JSON data structure, with special handling for array-type properties
    and NaN values.

    Args:
        full_json_data: Dictionary to update with new property values
        prop: Name of the property being added
        avg_prop: Averaged property value(s)
        smarts_targets: List of SMARTS patterns defining property scope

    Notes:
        - Preserves numpy arrays when SMARTS targets are present
        - Converts arrays to lists for JSON compatibility otherwise
        - Handles both scalar and array-type properties
        - Modifies full_json_data in place
    """
    if len(smarts_targets) > 0 or np.isnan(avg_prop).any():
        full_json_data[prop] = avg_prop
    else:
        full_json_data[prop] = avg_prop.tolist()

    return full_json_data

def dict_to_json(
    self,
    filepath: Union[str, Path],
    data: Dict[str, Any],
) -> None:
    """
    Save a dictionary to a JSON file with error handling.

    Args:
        filepath: Path to the output JSON file
        data: Dictionary to serialize to JSON

    Raises:
        IOError: If file cannot be written
        TypeError: If data contains non-serializable objects
        
    Notes:
        - Uses UTF-8 encoding for universal character support
        - Converts numpy arrays to lists automatically
        - Creates parent directories if they don't exist
    """
    filepath = Path(filepath)
    
    # Create parent directories if they don't exist
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Convert numpy arrays to lists for JSON serialization
        converted_data = data.copy()
        convert_ndarrays(converted_data)
        
        # Write to file with proper formatting
        with open(filepath, "w", encoding='utf-8') as outfile:
            json.dump(converted_data, outfile)
            
    except TypeError as e:
        type_error = f'Data contains non-serializable objects: {e}'
        self.args.log.write(type_error)
        raise TypeError(type_error)
    except IOError as e:
        io_error = f"Error writing to file {filepath}: {e}"
        self.args.log.write(io_error)
        raise IOError(io_error)
    except Exception as e:
        exception_error = f"Unexpected error saving JSON: {e}"
        self.args.log.write(exception_error)
        raise Exception(exception_error)

def get_mols_qdescp(qdescp_files: List[str]) -> List[MoleculeType]:
    """
    Extract molecule objects from input files.

    This function reads molecules from either SDF files or files containing
    SMILES strings, handling both formats transparently.

    Args:
        qdescp_files: List of paths to input files (SDF or files with SMILES)

    Returns:
        List of RDKit molecule objects with hydrogens added

    Notes:
        - First attempts to find SMILES strings in the file
        - Falls back to SDF parsing if no SMILES found
        - Automatically adds hydrogens to all molecules
        - Returns empty list if no valid molecules found
    """
    mol_list = []
    
    for file in qdescp_files:
        # Try to read file content
        try:
            with open(file, "r", encoding='utf-8') as f:
                lines = f.readlines()
        except Exception:
            continue
            
        # First try to find SMILES in the file
        smi_exist = False
        for i, line in enumerate(lines):
            if ">  <SMILES>" in line and i + 1 < len(lines):
                try:
                    smi = lines[i + 1].split()[0]
                    mol = Chem.MolFromSmiles(smi)
                    if mol is not None:
                        mol_list.append(Chem.AddHs(mol))
                        smi_exist = True
                        break
                except Exception:
                    continue
                    
        # If no SMILES found, try reading as SDF
        if not smi_exist:
            try:
                mols = load_sdf(file)
                if mols:
                    mol_list.append(mols[0])
            except Exception:
                continue
    
    return mol_list


def get_mol_assign(self,
        name_initial: str) -> MoleculeType:
    """
    Create RDKit molecule object from SDF file, supporting multiple formats.

    This function handles both CSEARCH-generated SDF files (with embedded SMILES)
    and regular SDF files. It attempts to extract SMILES first, then falls back
    to direct SDF parsing.

    Args:
        name_initial: Base name of the SDF file (without extension)

    Returns:
        RDKit molecule object with hydrogen atoms added

    Raises:
        FileNotFoundError: If SDF file doesn't exist
        ValueError: If molecule cannot be parsed from file
        
    Notes:
        - Prefers SMILES representation if available
        - Automatically adds hydrogen atoms
        - Handles both CSEARCH and standard SDF formats
    """
    sdf_path = Path(f'{name_initial}.sdf')
    
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    try:
        # Read SDF file content
        with open(sdf_path, "r", encoding='utf-8') as f:
            lines = f.readlines()

        # Try to find and use SMILES string first
        for i, line in enumerate(lines):
            if ">  <SMILES>" in line and i + 1 < len(lines):
                smiles = lines[i + 1].strip().split()[0]
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    return Chem.AddHs(mol)

        # Fall back to SDF parsing if no SMILES found
        mols = load_sdf(str(sdf_path))
        if not mols:
            val_error = f"x  WARNING! No valid molecules found in {sdf_path}"
            self.args.log.write(val_error)
            raise ValueError(val_error)
        
        return mols[0]  # Return first molecule

    except Exception as e:
        exc_error = f"Error processing SDF file {sdf_path}: {str(e)}"
        self.args.log.write(exc_error)
        raise ValueError(exc_error)

def auto_pattern(
    mol_list: List[Chem.Mol],
    smarts_targets: List[str],
) -> List[str]:
    """
    Automatically detect common SMARTS patterns across molecules.

    Args:
        mol_list: List of RDKit molecule objects
        smarts_targets: Existing SMARTS patterns or names (may be empty)

    Returns:
        Updated list of SMARTS patterns including automatically detected ones

    Notes:
        - Only includes functional groups or SMARTS patterns that appear exactly once per molecule
        - Preserves and validates existing SMARTS patterns from smarts_targets
        - idx_added_per_mol is a list of sets; each set contains atom indices
          (for that molecule) that are already assigned to a matched pattern
    """
    if not mol_list:
        return smarts_targets

    n_mols = len(mol_list)

    # Predefined functional groups (name -> SMARTS)
    func_groups = {
        'O=C[H]': '[CX3H1](=O)',
        'O=COC':  '[CX3](=O)O[#6]',
        'O=C[O-]': '[CX3](=O)[O-]',
        'O=C(N)': 'C(=O)[NX3;H2,H1,H0;!$(NC=O)]',
        'O=CO[H]': 'C(=O)[OX2H1]',
        'C#N': '[CX2]#[NX1]',
        'C=C': 'C=C',
        'C#C': 'C#C',
        'N=O': '[NX2]=[OX1]',
        'N=N': '[NX2]=[NX2]',
        'N#N': '[NX2]#[NX1]',
        'S=O': '[#16X3]=O',
        'O=S=O': 'S(=O)(=O)',
        'O-H': '[OX2H]',
        'N-H': '[NX3H]',
        'S-H': '[SX2H]'
    }

    # Index sets per molecule to track already-assigned atom indices
    idx_added_per_mol: List[set] = [set() for _ in range(n_mols)]

    # --- Step 0: Validate existing SMARTS and collect their atom indices per molecule ---
    validated_targets: List[str] = []
    for target in smarts_targets:
        smarts = func_groups.get(target, target)
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue

        # must appear exactly once in every molecule
        appear_once_in_all = True
        all_matches_per_mol = []
        for mol in mol_list:
            matches = mol.GetSubstructMatches(pattern)
            all_matches_per_mol.append(matches)
            if len(matches) != 1:
                appear_once_in_all = False
                break

        if not appear_once_in_all:
            continue

        # collect matching atom indices per molecule (flatten tuples)
        for i, matches in enumerate(all_matches_per_mol):
            for match in matches:
                for idx in match:
                    idx_added_per_mol[i].add(idx)

        validated_targets.append(target)

    smarts_targets = validated_targets.copy()

    # --- Step 1: Detect predefined functional groups (only if exactly once per mol) ---
    for group_name, smarts in func_groups.items():
        if group_name in smarts_targets:
            continue  # already validated/present

        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue

        appear_once_in_all = True
        all_matches_per_mol = []
        for mol in mol_list:
            matches = mol.GetSubstructMatches(pattern)
            all_matches_per_mol.append(matches)
            if len(matches) != 1:
                appear_once_in_all = False
                break

        if not appear_once_in_all:
            continue

        smarts_targets.append(group_name)
        for i, matches in enumerate(all_matches_per_mol):
            for match in matches:
                for idx in match:
                    idx_added_per_mol[i].add(idx)

    # --- Step 2: Add individual unique atoms (atoms appearing exactly once in each mol) ---
    unique_atoms_per_mol = []
    for mol in mol_list:
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        single_occurrence = {sym for sym in set(symbols) if symbols.count(sym) == 1}
        unique_atoms_per_mol.append(single_occurrence)

    if unique_atoms_per_mol:
        common_unique_atoms = set.intersection(*unique_atoms_per_mol)
    else:
        common_unique_atoms = set()

    for atom_symbol in common_unique_atoms:
        if atom_symbol in smarts_targets:
            continue

        overlaps_any = False
        for i, mol in enumerate(mol_list):
            atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == atom_symbol]
            if len(atom_indices) != 1:
                overlaps_any = True
                break
            idx = atom_indices[0]
            if idx in idx_added_per_mol[i]:
                overlaps_any = True
                break

        if not overlaps_any:
            smarts_targets.append(atom_symbol)
            for i, mol in enumerate(mol_list):
                idx = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == atom_symbol)
                idx_added_per_mol[i].add(idx)

    return smarts_targets

def remove_invalid_smarts(
    self,
    mol_list: List[MoleculeType],
    smarts_targets: List[str],
    min_match_percent: float = 0.75
) -> List[str]:
    """
    Remove SMARTS patterns that don't match consistently across molecules.

    This function validates each SMARTS pattern against a list of molecules
    and removes patterns that:
    1. Don't parse correctly as SMARTS
    2. Match fewer than min_match_percent of molecules
    3. Match multiple times in a single molecule

    Args:
        self: Instance containing logging capabilities
        mol_list: List of RDKit molecule objects to check
        smarts_targets: List of SMARTS patterns to validate
        min_match_percent: Minimum fraction of molecules that must match (default: 0.75)

    Returns:
        List of validated SMARTS patterns

    Notes:
        - Handles both atom mapping numbers and SMARTS patterns
        - Supports both mapped and unmapped molecules
        - Automatically logs warnings for invalid or inconsistent patterns
    """
    def _clean_pattern(pattern: str) -> str:
        """Convert lists to pattern (for direct module calls) and remove quotes from pattern string."""
        return pattern.strip("'").strip('"')

    def _find_mapped_atom_matches(mol: MoleculeType, atom_num: int) -> List[Tuple[int, ...]]:
        """Find matches for a specific atom number in mapped/unmapped molecules."""
        mol_idxs = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
        
        # Handle unmapped molecules (SDF files)
        if len(set(mol_idxs)) == 1 and mol_idxs[0] == 0:
            for i, atom in enumerate(mol.GetAtoms()):
                if i == atom_num - 1:  # Convert from 1-based to 0-based indexing
                    return [(int(atom.GetIdx()),)]
        
        # Handle mapped molecules (SMILES)
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_num:
                return [(int(atom.GetIdx()),)]
        
        return []

    def _find_smarts_matches(mol: MoleculeType, pattern: str) -> List[Tuple[int, ...]]:
        """Find all matches for a SMARTS pattern in a molecule."""
        try:
            return mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        except:
            try:
                return mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
            except:
                return []

    def _check_pattern(pattern: str, mol: MoleculeType) -> bool:
        """Check if pattern matches exactly once in molecule."""
        pattern = _clean_pattern(pattern)
        
        # Handle atom numbers vs SMARTS patterns
        if pattern.isdigit():
            matches = _find_mapped_atom_matches(mol, int(pattern))
        else:
            matches = _find_smarts_matches(mol, pattern)
        
        return len(matches) == 1

    # Track patterns to remove
    patterns_to_remove = set()
    required_matches = len(mol_list) * min_match_percent

    # Validate each pattern
    for pattern in smarts_targets:
        matches_found = sum(1 for mol in mol_list if _check_pattern(pattern, mol))
        
        if matches_found < required_matches:
            patterns_to_remove.add(pattern)
            if pattern in self.args.qdescp_atoms:
                if matches_found == 0:
                    self.args.log.write(
                        f"x  WARNING! SMARTS pattern {pattern} was not specified correctly. "
                        "Make sure to use proper format: \"[C]\" for atoms, \"[C=N]\" for bonds, etc."
                    )
                else:
                    self.args.log.write(
                        f"x  WARNING! SMARTS pattern {pattern} is not present (or matches "
                        f"multiple times) in {min_match_percent*100}% or more of the molecules. "
                        "Atomic descriptors will not be generated for this pattern."
                    )

    # Remove invalid patterns
    valid_patterns = [p for p in smarts_targets if p not in patterns_to_remove]

    # Log automatically detected patterns
    for pattern in valid_patterns:
        if pattern not in self.args.qdescp_atoms:
            self.args.log.write(
                f'o  Pattern "{pattern}" validated, it will be used for atomic descriptor calculations'
            )

    return valid_patterns


def get_atom_matches(
    self,
    pattern: str,
    mol: MoleculeType
) -> Tuple[List[Tuple[int, ...]], Optional[str]]:
    """
    Find atoms or groups matching a pattern in a molecule.

    This function handles both atom mapping numbers and SMARTS patterns
    to identify specific atoms or functional groups in a molecule.

    Args:
        self: Instance containing logging capabilities
        pattern: Either an atom mapping number or SMARTS pattern
        mol: RDKit molecule object to search

    Returns:
        Tuple containing:
        - List of matching atom index tuples
        - Mapped atom number if pattern is a digit, None otherwise

    Notes:
        - Handles both mapped atoms (by number) and SMARTS patterns
        - Returns [] for no matches
        - Provides special handling for mapped vs unmapped molecules
    """
    matches = []
    idx_set = None
    
    # Check if pattern is a number (mapped atom) or SMARTS pattern
    if not str(pattern).isalpha() and str(pattern).isdigit():
        idx_set = pattern

        # for non-mapped mols (i.e. SDF input files)
        mol_idxs = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
        if len(set(mol_idxs)) == 1 and mol_idxs[0] == 0:
            for i,atom in enumerate(mol.GetAtoms()):
                if i == int(pattern)-1: # atoms in SDF starts in index 1, but Python starts in idx 0
                    pattern_idx = int(atom.GetIdx())
                    matches = ((int(pattern_idx),),)
        # for mapped SMILES
        else:
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == int(pattern):
                    pattern_idx = int(atom.GetIdx())
                    matches = ((int(pattern_idx),),)

    else: 
        try:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        except:
            try:
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
            except:
                self.args.log.write(f"x  WARNING! SMARTS pattern was not specified correctly! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")

    return matches, idx_set


def sort_atom_types(
    matches: List[Tuple[int, ...]], 
    mol: MoleculeType
) -> Tuple[int, List[int]]:
    """
    Sort atoms based on their types and connectivity.

    This function processes RDKit substructure matches and sorts atoms based
    on their element type and connectivity to ensure consistent ordering
    across different molecules.

    Args:
        matches: List of atom index tuples from RDKit substructure match
        mol: RDKit molecule object containing the matched atoms

    Returns:
        Tuple containing:
        - Number of unique atom types in the match
        - List of sorted atom indices

    Notes:
        - For single atom type matches: sorts by number of neighbors
        - For multiple atom types: sorts by atomic number
        - Returns empty list if matches or molecule is None
    """
    if not matches or not mol:
        return 0, []
        
    atom_indices = list(matches[0])
    atom_types = []
    
    for atom_idx in atom_indices:
        atom_types.append(mol.GetAtoms()[atom_idx].GetSymbol())

    n_types = len(set(atom_types))
    if n_types == 1:
        # For single atom type, sort by connectivity
        sorted_indices = sorted(
            atom_indices, 
            key=lambda idx: len(mol.GetAtoms()[idx].GetNeighbors())
        )
    else:
        # For multiple atom types, sort by atomic number
        sorted_indices = sorted(
            atom_indices, 
            key=lambda idx: mol.GetAtoms()[idx].GetAtomicNum()
        )

    return n_types, sorted_indices


def get_prefix_atom_props(
    sorted_indices: List[int],
    mol: MoleculeType,
    pattern: str,
    idx_set: Optional[str]
) -> List[str]:
    """
    Generate unique match names for each atom in a functional group.

    This function creates unique identifiers for atoms in a matched functional
    group, handling both mapped atoms and SMARTS patterns.

    Args:
        sorted_indices: List of sorted atom indices from the match
        mol: RDKit molecule object containing the atoms
        pattern: Original SMARTS pattern or atom number
        idx_set: Optional explicit atom mapping number

    Returns:
        List of unique match names for each atom

    Notes:
        - Handles both mapped atoms and SMARTS pattern matches
        - Creates unique names for multiple instances of same atom type
        - Preserves atom ordering from sorted_indices
    """
    match_names = []
    atom_counters = {}

    # Iterate through matched atoms in their sorted order
    for atom_idx in sorted_indices:
        atom_type = mol.GetAtoms()[atom_idx].GetSymbol()

        # Count atoms of this type in the match
        n_atoms_of_type = sum(
            1 for idx in sorted_indices 
            if mol.GetAtoms()[idx].GetSymbol() == atom_type
        )

        # Initialize counter for new atom type
        if atom_type not in atom_counters:
            atom_counters[atom_type] = 1
        # If the pattern is a single atom number, we include that number in the match name
        if str(pattern).isdigit():  # This means the pattern is an atom number
            match_name = f'Atom_{pattern}'
        else:
            # If it's a SMARTS pattern or more than one atom
            if pattern in periodic_table():
                # Case where it's just an atom type without SMARTS
                if n_atoms_of_type == 1:
                    match_name = f'{atom_type}'
                else:
                    match_name = f'{atom_type}_{atom_counters[atom_type]}'
            else:
                if pattern[0] == '#':
                    match_name = f'{atom_type}'
                # Regular SMARTS pattern handling
                elif n_atoms_of_type == 1:
                    if idx_set is None:
                        match_name = f'{pattern}_{atom_type}'
                    else:
                        match_name = f'{pattern}_{atom_type}{idx_set}'
                else:
                    if idx_set is None:
                        match_name = f'{pattern}_{atom_type}_{atom_counters[atom_type]}'
                    else:
                        match_name = f'{pattern}_{atom_type}{idx_set}_{atom_counters[atom_type]}'
        
        atom_counters[atom_type] += 1  # Increment counter for that atom type

        match_names.append(match_name)

    return match_names


def update_atom_props_json(
    sorted_indices: List[int],
    match_names: List[str],
    atom_props: List[str],
    json_data: Dict[str, Any],
    prefixes_atom_prop: List[str],
    pattern: str,
    n_types: int
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Update JSON data with atomic descriptors and derived properties.

    This function assigns atomic descriptors to matched atoms and computes
    aggregate properties (min/max) for functional groups.

    Args:
        sorted_indices: List of sorted atom indices from the match
        match_names: List of unique identifiers for each atom
        atom_props: List of atomic property names to process
        json_data: Dictionary containing property values to update
        prefixes_atom_prop: List of existing property prefixes
        pattern: SMARTS pattern or atom number identifier
        n_types: Number of unique atom types in the match

    Returns:
        Tuple containing:
        - Updated list of property prefixes
        - Updated JSON data dictionary

    Notes:
        - Assigns properties per atom with unique identifiers
        - Computes min/max for groups with same atom type
        - Handles missing values gracefully
        - Preserves existing prefixes
    """
    # Assign properties to individual atoms
    for atom_idx, match_name in zip(sorted_indices, match_names):
        idx_xtb = atom_idx
        for prop in atom_props:
            try:
                # Assign property value or None if not available
                prop_key = f'{match_name}_{prop}'
                json_data[prop_key] = (
                    json_data[prop][idx_xtb] 
                    if json_data[prop] is not None 
                    else None
                )
                
                # Add prefix if new
                prefix = f'{match_name}_'
                if prefix not in prefixes_atom_prop:
                    prefixes_atom_prop.append(prefix)
            except (KeyError, IndexError):
                # Skip missing or invalid properties
                continue

    # Calculate aggregate properties for groups of same atom type
    if len(match_names) > 1 and n_types == 1:
        for prop in atom_props:
            # Collect values for this property across all atoms
            prop_values = [
                json_data[f'{name}_{prop}']
                for name in match_names
            ]

            # Calculate min/max if all values are valid
            if None not in prop_values:
                json_data[f'{pattern}_max_{prop}'] = max(prop_values)
                json_data[f'{pattern}_min_{prop}'] = min(prop_values)
            else:
                json_data[f'{pattern}_max_{prop}'] = None
                json_data[f'{pattern}_min_{prop}'] = None

            # Add new prefixes if needed
            for prefix in [f'{pattern}_max_', f'{pattern}_min_']:
                if prefix not in prefixes_atom_prop:
                    prefixes_atom_prop.append(prefix)

    return prefixes_atom_prop, json_data