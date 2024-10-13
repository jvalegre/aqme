######################################################.
#        This file stores QDESCP functions           #
######################################################.

import json
import sys
import os
import numpy as np
import pandas as pd
import ast
import math
import rdkit
from rdkit.Chem import Descriptors
import warnings
warnings.filterwarnings('ignore')
from morfeus import SASA, Dispersion, BuriedVolume, ConeAngle, SolidAngle, Pyramidalization, read_xyz, read_geometry

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15
Hartree = 27.2114 #eV UNIT CONVERSION from 1 hatree = 27.2114 eV

def convert_ndarrays(data):
    """Convert properties that are of type ndarray to lists to serialize them in JSON"""
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
    weights (list of floats): List of Boltzmann weights corresponding to the input energies.
    """
    # Shift energies so that the minimum energy is zero (to prevent negative exponents from dominating)
    energ = [number - min(energy) for number in energy]

    boltz_sum = 0.0 # Initialize the sum of Boltzmann factors to zero

    # Calculate the sum of Boltzmann factors for all energies
    for e in energ:
        # Apply the Boltzmann factor formula: exp(-energy / (k * T))
        # where J_TO_AU is a conversion factor, GAS_CONSTANT is the gas constant, and T is temperature
        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
    
    weights = [] # Initialize the list to store Boltzmann weights

    # Calculate each individual Boltzmann weight by normalizing with the total sum of Boltzmann factors
    for e in energ:
        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
        weights.append(weight)

    # Return the list of normalized Boltzmann weights
    return weights


def get_boltz_props(json_files, name, boltz_dir, calc_type, self, mol_props, atom_props, smarts_targets, 
                    denovo_mols, denovo_atoms,interpret_mols, interpret_atoms, mol=None):
    """
    Retrieves the properties from json files and gives Boltzmann averaged properties for rdkit, NMR and morfues descriptors.
    """
    # Ensure smarts_targets is a list even if None
    if smarts_targets is None:
        smarts_targets = []

    def average_properties(boltz, prop_list, is_atom_prop=True):
        """Calculate average properties based on Boltzmann weights."""
        return average_prop_atom(boltz, prop_list) if is_atom_prop else average_prop_mol(boltz, prop_list)

    def update_avg_json_data(avg_json_data, prop, avg_prop, smarts_targets):
        """Update avg_json_data with averaged properties."""
        if len(smarts_targets) > 0 or np.isnan(avg_prop).any():
            avg_json_data[prop] = avg_prop
        else:
            avg_json_data[prop] = avg_prop.tolist()

    # Calculate Boltzmann weights
    energy = []
    for _, json_file in enumerate(json_files):
        json_data = read_json(json_file)
        energy.append(json_data["total energy"] if calc_type.lower() == "xtb" else json_data["optimization"]["scf"]["scf energies"][-1])

    boltz = get_boltz(energy)
    avg_json_data = {}
    denovo_json_data = {}
    interpret_json_data = {}
    atomic_props = True

    # Get weighted atomic properties
    for i, prop in enumerate(atom_props):
        if i == 0:  # Filter to ensure all molecules have atomic properties for qdescp_atoms option
            for json_file in json_files:
                json_data = read_json(json_file)
                for atom_prop in atom_props:
                    if atom_prop not in json_data:
                        atomic_props = False
        if atomic_props:
            try:
                prop_list = [read_json(json_file)[prop] for json_file in json_files]
                avg_prop = average_properties(boltz, prop_list)
            except KeyError:
                avg_json_data[prop] = np.nan
        else:
            avg_prop = np.nan
        update_avg_json_data(avg_json_data, prop, avg_prop, smarts_targets)

    # Get weighted molecular properties from XTB
    for prop in mol_props:
        try:
            prop_list = [read_json(json_file)[prop] for json_file in json_files]
            avg_prop = average_properties(boltz, prop_list, is_atom_prop=False)
            avg_json_data[prop] = avg_prop
        except KeyError:
            avg_json_data[prop] = np.nan

    # Get denovo atomic properties
    for i, prop in enumerate(denovo_atoms):
        if atomic_props:
            try:
                prop_list = [read_json(json_file)[prop] for json_file in json_files]
                avg_prop = average_properties(boltz, prop_list)
            except KeyError:
                avg_json_data[prop] = np.nan
        else:
            avg_prop = np.nan
        update_avg_json_data(denovo_json_data, prop, avg_prop, smarts_targets)

    # Get denovo molecular properties
    for prop in denovo_mols:
        try:
            prop_list = [read_json(json_file)[prop] for json_file in json_files]
            avg_prop = average_properties(boltz, prop_list, is_atom_prop=False)
            denovo_json_data[prop] = avg_prop
        except KeyError:
            avg_json_data[prop] = np.nan

    # Get interpret atomic properties
    for i, prop in enumerate(interpret_atoms):
        if atomic_props:
            try:
                prop_list = [read_json(json_file)[prop] for json_file in json_files]
                avg_prop = average_properties(boltz, prop_list)
            except KeyError:
                avg_json_data[prop] = np.nan

        else:
            avg_prop = np.nan
        update_avg_json_data(interpret_json_data, prop, avg_prop, smarts_targets)

    # Get interpret molecular properties
    for prop in interpret_mols:
        try:
            prop_list = [read_json(json_file)[prop] for json_file in json_files]
            avg_prop = average_properties(boltz, prop_list, is_atom_prop=False)
            interpret_json_data[prop] = avg_prop
        except KeyError:
            avg_json_data[prop] = np.nan

    # Calculate RDKit descriptors if molecule is provided
    if mol is not None:
        # Calculate all RDKit properties for avg_json_data
        avg_json_data,_,_ = get_rdkit_properties(self,avg_json_data, mol)
        
        # Calculate selected RDKit properties for denovo_json_data
        _,denovo_rdkit_json_data,_ = get_rdkit_properties(self,{}, mol)
        
        # Merge selected RDKit properties with denovo_json_data
        denovo_json_data.update(denovo_rdkit_json_data)

        # Calculate selected RDKit properties for interpret_json_data
        _,_,interpret_rdkit_json_data = get_rdkit_properties(self,{}, mol)
        
        # Merge selected RDKit properties with interpret_json_data
        interpret_json_data.update(interpret_rdkit_json_data)
    
    convert_ndarrays(avg_json_data)
    convert_ndarrays(denovo_json_data)
    convert_ndarrays(interpret_json_data)

    # Save the averaged properties to a file
    final_boltz_file = os.path.join(boltz_dir, f"{name}_full_boltz.json")
    with open(final_boltz_file, "w") as outfile:
        json.dump(avg_json_data, outfile)
    
    # Save the denovo properties to a second file
    final_denovo_file = os.path.join(boltz_dir, f"{name}_denovo_boltz.json")
    with open(final_denovo_file, "w") as outfile:
        json.dump(denovo_json_data, outfile)

    # Save the interpret properties to a second file
    final_interpret_file = os.path.join(boltz_dir, f"{name}_interpret_boltz.json")
    with open(final_interpret_file, "w") as outfile:
        json.dump(interpret_json_data, outfile)


def get_boltz_props_nmr(json_files,name,boltz_dir,self,atom_props,nmr_atoms=None,nmr_slope=None,nmr_intercept=None,nmr_experim=None):

    """
    Function to process NMR properties and calculate Boltzmann averaged properties.
    """

    # Load experimental NMR data
    if nmr_experim is not None:
        try:
            exp_data = pd.read_csv(nmr_experim)
        except Exception as e:
            self.args.log.write(f'\nx  Error loading experimental NMR shifts file: {nmr_experim}.\n{e}')
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
        with open(json_file, "w") as outfile:
            json.dump(json_data, outfile)

    # Calculate Boltzmann factors
    boltz = get_boltz(energy)

    # Process atomic properties (NMR Chemical Shifts)
    avg_json_data = {}
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
        avg_json_data[prop] = dictavgprop

    # Merge with experimental data if available and calculate error
    if nmr_experim is not None:
        list_shift = avg_json_data[prop]
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
    with open(final_boltz_file, "w") as outfile:
        json.dump(avg_json_data, outfile)

def get_chemical_shifts(json_data, nmr_atoms, nmr_slope, nmr_intercept):
    """
    Retrieves and scales NMR shifts from json files
    """

    if not isinstance(nmr_atoms, list):
        nmr_atoms = ast.literal_eval(nmr_atoms)
    if not isinstance(nmr_slope, list):
        nmr_slope = ast.literal_eval(nmr_slope)
    if not isinstance(nmr_intercept, list):
        nmr_intercept = ast.literal_eval(nmr_intercept)

    atoms = json_data["atoms"]["elements"]["number"]
    tensor = json_data["properties"]["NMR"]["NMR isotopic tensors"]
    shifts = {}
    i = 0
    for atom, ten in zip(atoms, tensor):
        if atom in nmr_atoms:
            # assigning values from arrays
            index = nmr_atoms.index(atom)
            slope_nuc = nmr_slope[index]
            intercept_nuc = nmr_intercept[index]

            scaled_nmr = (intercept_nuc - ten) / (-slope_nuc)
            shifts[i] = scaled_nmr
        else:
            pass
        i += 1

    return shifts

def average_prop_atom(weights, prop):
    """
    Returns the atomic properties averaged using the Boltzmann average.

    Parameters:
    weights (list of floats): Boltzmann weights calculated for each configuration.
    prop (list): List of atomic properties (either floats or lists of values) for each configuration.

    Returns:
    boltz_res: The Boltzmann-weighted average of the atomic properties.
    """
    # List to store the Boltzmann-weighted properties
    boltz_avg = []
    
    # Loop through each property and its corresponding weight
    if len(prop) == len(weights): # if xTB fails a calculation in a conformer, that property won't be included
        for i, p in enumerate(prop):
            # If the property is a single float or int, multiply it by the corresponding weight
            if isinstance(p, (float, int)):
                boltz_avg.append(p * weights[i])
            # If the property is a list (e.g., multi-dimensional property), multiply each element by the weight
            elif isinstance(p, list):
                boltz_avg.append([0 if number is None else number * weights[i] for number in p])
            # If the property type is not recognized, append 0 as a placeholder
            else:
                boltz_avg.append(np.nan)

        # Check if the first element in the list is a list (indicating multi-dimensional properties)
        try:
            if isinstance(boltz_avg[0], list):
                # If so, sum along the axis (i.e., sum element-wise for multi-dimensional properties)
                boltz_res = np.sum(boltz_avg, axis=0)
                # Round each element in the result to 4 decimal places if it's a NumPy array
                boltz_res = np.round(boltz_res, 4)
            else:
                # Otherwise, sum all the weighted values as a scalar and round to 4 decimal places
                boltz_res = round(sum(boltz_avg), 4)
        except (ValueError,IndexError): # in some cases, there are atoms that are missing properties because of xTB calculation errors
            boltz_res = np.nan

    else:
        boltz_res = np.nan
    
    # Return the Boltzmann-weighted result
    return boltz_res

def average_prop_atom_nmr(weights, prop):
    """
    Returns Boltzmann averaged atomic properties
    """

    boltz_avg = []
    for i, p in enumerate(prop):
        if p == 'NaN':
            boltz_avg = 'NaN'
            break
        boltz_avg.append([number * weights[i] for number in p])
    if boltz_avg == 'NaN':
        boltz_res = 'NaN'
    else:
        boltz_res = np.sum(boltz_avg, 0)
    return boltz_res

def average_prop_mol(weights, prop):
    """
    Returns Boltzmann averaged molecular properties, rounded to 4 decimal places.

    Parameters:
    weights (list of floats): Boltzmann weights for each configuration.
    prop (list of floats): List of molecular properties corresponding to each configuration.

    Returns:
    boltz_avg (float): Boltzmann-weighted average of the molecular properties, rounded to 4 decimal places.
    """
    
    # Initialize the Boltzmann-weighted average to 0.0
    boltz_avg = 0.0
    
    # Loop through each property and its corresponding weight
    for i, p in enumerate(prop):
        # If the property is 'NaN', break the loop and return 'NaN'
        if p == 'NaN':
            boltz_avg = 'NaN'
            break
        # Otherwise, calculate the weighted sum of properties
        boltz_avg += p * weights[i]
    
    # If the result is a valid number, round it to 4 decimal places
    if boltz_avg != 'NaN':
        boltz_avg = round(boltz_avg, 4)

    # Return the Boltzmann-weighted average, rounded to 4 decimal places (or 'NaN' if encountered)
    return boltz_avg

def get_rdkit_properties(self,avg_json_data, mol):
    """
    Calculates RDKit molecular descriptors
    """

    #Level: denovo descriptors
    denovo_json_data = {}
    #Level: interpret descriptors
    interpret_json_data = {}

    try:
        #level: full
        descrs = Descriptors.CalcMolDescriptors(mol)
        for descr in descrs:
            if descrs[descr] != np.nan and str(descrs[descr]).lower() != 'nan':
                avg_json_data[descr] = descrs[descr]
        
        # descriptors for the level_ denovo
        denovo_json_data["MolLogP"] = descrs.get("MolLogP", None)
        # descriptors for the level: interpret
        interpret_json_data["MolLogP"] = descrs.get("MolLogP", None)

    # For older versions of RDKit 
    except AttributeError:
        self.args.log.write(f"x  WARNING! Install a newer version of RDKit to get all the RDKit descriptors in the databse with all the descriptors. You can use: 'pip install rdkit --upgrade'.")
        avg_json_data["MolLogP"] = rdkit.Chem.Descriptors.MolLogP(mol)

        # descriptors for the level_ denovo
        denovo_json_data["MolLogP"] = avg_json_data["MolLogP"]
        # descriptors for the level: interpret
        interpret_json_data["MolLogP"] = avg_json_data["MolLogP"]


    return avg_json_data, denovo_json_data, interpret_json_data

def read_gfn1(file,self):
    """
    Read .gfn1 output file created from xTB and return parsed data.

    Parameters:
    file (str): Path to the .gfn1 file.

    Returns:
    dict: Parsed data containing Mulliken charges, CM5 charges, and proportions.
          Returns None if there's an issue with the file.
    """

    # Check if the file exists
    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    # Open and read the file safely
    with open(file, "r") as f:
        data = f.readlines()

    # Ensure the file contains data
    if not data:
        self.args.log.write(f"x  WARNING! The file {file} is empty.")
        return None

    # Initialize variables for data extraction
    start, end = None, None

    # Find the start of the Mulliken/CM5 charges section
    for i, line in enumerate(data):
        if "Mulliken/CM5 charges" in line:
            start = i + 1
            break

    # Find the end of the section (start of Wiberg/Mayer or generalized Born)
    if start is not None:
        for j in range(start, len(data)):
            if "Wiberg/Mayer (AO) data" in data[j] or "generalized Born model" in data[j]:
                end = j - 1
                break
    else:
        self.args.log.write(f"x  WARNING! Mulliken/CM5 charges section not found in {file}.")
        return None

    # Ensure both start and end are found
    if start is None or end is None:
        self.args.log.write(f"x  WARNING! Required sections not found in {file}.")
        return None

    # Extract the relevant data between the start and end lines
    pop_data = data[start:end]

    # Initialize lists to store the parsed data
    mulliken, cm5, s_prop, p_prop, d_prop = [], [], [], [], []

    # Process each line in the extracted data section
    for line in pop_data:
        try:
            item = line.split()
            # Extract and round the required values from each line
            q_mull = round(float(item[-5]), 5)
            q_cm5 = round(float(item[-4]), 5)
            s_prop_ind = round(float(item[-3]), 3)
            p_prop_ind = round(float(item[-2]), 3)
            d_prop_ind = round(float(item[-1]), 3)
            mulliken.append(q_mull)
            cm5.append(q_cm5)
            s_prop.append(s_prop_ind)
            p_prop.append(p_prop_ind)
            d_prop.append(d_prop_ind)
        except (ValueError, IndexError) as e:
            # Handle errors related to parsing the line
            self.args.log.write(f"x  WARNING! Error parsing line in {file}: {line}. Error: {e}")
            return None

    # Store the parsed data in a dictionary and return it
    localgfn1 = {
        "mulliken charges": mulliken,
        "cm5 charges": cm5,
        "s proportion": s_prop,
        "p proportion": p_prop,
        "d proportion": d_prop,
    }

    return localgfn1

def read_wbo(file,self):
    """
    Read wbo output file created from xTB. Return data.
    """

    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    with open(file, "r") as f:
        data = f.readlines()

    bonds, wbos = [], []
    for line in data:
        item = line.split()  # Split the line into components
        bond = [int(item[0]), int(item[1])]  # Extract bond indices
        wbo = round(float(item[2]), 3)  # Extract and round the WBO value
        bonds.append(bond)  # Add bond to list
        wbos.append(wbo)  # Add WBO to list

    return bonds, wbos


def calculate_global_CDFT_descriptors(file,self):
    """
    Read .gfn1 output file created from xTB and calculate CDFT descriptors with FDA approximations part 1.
    """

    # Check if the file exists
    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    with open(file, "r") as f:
        data = f.readlines()
        
    # Initialize variables
    delta_SCC_IP, delta_SCC_EA, electrophilicity_index = None, None, None
    chemical_hardness, chemical_softness = None, None
    chemical_potential, mulliken_electronegativity = None, None
    electrodonating_power_index, electroaccepting_power_index = None, None
    intrinsic_reactivity_index = None
    electrofugality, nucleofugality, nucleophilicity_index, net_electrophilicity = None, None, None, None

    # Extract relevant values from the file
    for line in data:
        if "delta SCC IP (eV):" in line:
            delta_SCC_IP = float(line.split()[-1])
        elif "delta SCC EA (eV):" in line:
            delta_SCC_EA = float(line.split()[-1])
        elif "Global electrophilicity index (eV):" in line:
            electrophilicity_index = float(line.split()[-1])

    # Check if required descriptors were found
    if delta_SCC_IP is not None and delta_SCC_EA is not None and electrophilicity_index is not None:
        # Calculate CDFT descriptors
        chemical_hardness = round((delta_SCC_IP - delta_SCC_EA), 4)
        chemical_potential = round(-(delta_SCC_IP + delta_SCC_EA) / 2, 4)
        mulliken_electronegativity = round(-chemical_potential, 4)
        electrofugality = round(-delta_SCC_EA + electrophilicity_index, 4)
        nucleofugality = round(delta_SCC_IP + electrophilicity_index, 4)

        if chemical_hardness != 0:
            chemical_softness = round(1 / chemical_hardness, 4)
            electrodonating_power_index = round(((delta_SCC_IP + 3 * delta_SCC_EA)**2) / (8 * chemical_hardness), 4)
            electroaccepting_power_index = round(((3 * delta_SCC_IP + delta_SCC_EA)**2) / (8 * chemical_hardness), 4)
            intrinsic_reactivity_index = round((delta_SCC_IP + delta_SCC_EA) / chemical_hardness, 4)

            if electroaccepting_power_index != 0:
                nucleophilicity_index = round(10 / electroaccepting_power_index, 4)

        net_electrophilicity = round((electrodonating_power_index - electroaccepting_power_index), 4)

    else:
        self.args.log.write(f"x  WARNING! delta_SCC_IP, delta_SCC_EA, or electrophilicity_index were not found in the file. Global Conceptual DFT descriptors cannot be fully calculated.")

    # Collect descriptors into a dictionary
    cdft_descriptors = {
        "IP": delta_SCC_IP,
        "EA": delta_SCC_EA,
        "Electrophil. idx": electrophilicity_index,
        "Hardness": chemical_hardness,
        "Softness": chemical_softness,
        "Chem. potential": chemical_potential,
        "Electronegativity": mulliken_electronegativity,
        "Electrodon. power idx": electrodonating_power_index,
        "Electroaccep. power idx": electroaccepting_power_index,
        "Nucleophilicity idx": nucleophilicity_index,
        "Electrofugality": electrofugality,
        "Nucleofugality": nucleofugality,
        "Intrinsic React. idx": intrinsic_reactivity_index,
        "Net Electrophilicity": net_electrophilicity
    }

    return cdft_descriptors


def calculate_global_CDFT_descriptors_part2(file, file_Nminus1, file_Nminus2, file_Nplus1, file_Nplus2, cdft_descriptors,self):
    """
    Read .gfn1 output file created from xTB and calculate CDFT descriptors with FDA approximations part 2
    """
    corr_xtb = 4.8455  # correction from XTB

    def extract_scc_energy(lines, filename):
        """
        Extract SCC energy from the file. If not found, return None and log a warning.
        """
        for line in lines:
            if "SCC energy" in line:
                return float(line.split()[3])
        
        self.args.log.write(f"x  WARNING! Could not find SCC energy value in the file: {filename}")
        return None

    try:
        # Open and read files
        with open(file, "r") as f:
            data = f.readlines()
        with open(file_Nminus1, "r") as f1:
            data1 = f1.readlines()
        with open(file_Nminus2, "r") as f2:
            data2 = f2.readlines()
        with open(file_Nplus1, "r") as f3:
            data3 = f3.readlines()
        with open(file_Nplus2, "r") as f4:
            data4 = f4.readlines()
    except Exception as e:
        self.args.log.write(f"x  WARNING! An error occurred while processing {file}: {e}")
        return None

    # Extract SCC energies, handle cases where None is returned
    scc_energy = extract_scc_energy(data, file)
    scc_energy_Nminus1 = extract_scc_energy(data1, file_Nminus1)
    scc_energy_Nminus2 = extract_scc_energy(data2, file_Nminus2)
    scc_energy_Nplus1 = extract_scc_energy(data3, file_Nplus1)
    scc_energy_Nplus2 = extract_scc_energy(data4, file_Nplus2)

    if None in [scc_energy, scc_energy_Nminus1, scc_energy_Nminus2, scc_energy_Nplus1, scc_energy_Nplus2]:
        self.args.log.write(f"x  WARNING! Missing SCC energy in one or more files. Cannot calculate global CDFT descriptors.")
        return None

    # Convert SCC energies to Hartree
    scc_energy *= Hartree
    scc_energy_Nminus1 *= Hartree
    scc_energy_Nminus2 *= Hartree
    scc_energy_Nplus1 *= Hartree
    scc_energy_Nplus2 *= Hartree

    # Extract required CDFT descriptors
    delta_SCC_IP = cdft_descriptors.get("IP")
    delta_SCC_EA = cdft_descriptors.get("EA")
    chemical_hardness = cdft_descriptors.get("Hardness")

    if None in [delta_SCC_IP, delta_SCC_EA, chemical_hardness]:
        self.args.log.write("x  WARNING! Missing required CDFT descriptors (IP, EA, or Hardness).")
        return None

    # Initialize variables
    Vertical_second_IP, Vertical_second_EA = None, None
    hyper_hardness, Global_hypersoftness = None, None
    Electrophilic_descriptor, w_cubic = None, None

    # Calculations if all SCC energies and descriptors are available:
    # Calculating the following descriptors
        # 1) Second IP
    Vertical_second_IP = round((((scc_energy_Nminus2 - scc_energy_Nminus1) - corr_xtb)), 4)
        # 2) Second EA
    Vertical_second_EA = round((((scc_energy_Nplus1 - scc_energy_Nplus2) + corr_xtb)), 4)
        # 3) Hyperhardnes
    hyper_hardness = round((-((0.5) * (delta_SCC_IP + delta_SCC_EA - Vertical_second_IP - Vertical_second_EA))), 4)

    if chemical_hardness != 0:
        Global_hypersoftness = round((hyper_hardness / ((chemical_hardness) ** 3)), 4)

        # 4) Electrophilic descriptor calculations
    A = ((scc_energy_Nplus1 - scc_energy) + corr_xtb)
    c = (Vertical_second_IP - (2 * delta_SCC_IP) + A) / ((2 * Vertical_second_IP) - delta_SCC_IP - A)
    a = -((delta_SCC_IP + A) / 2) + (((delta_SCC_IP - A) / 2) * c)
    b = ((delta_SCC_IP - A) / 2) - (((delta_SCC_IP + A) / 2) * c)
    Gamma = (-3 * c) * (b - (a * c))
    Eta = 2 * (b - (a * c))
    chi = -a
    Mu = a

    discriminant = Eta ** 2 - (2 * Gamma * Mu)  # Checking the square root
    if discriminant < 0:
        self.args.log.write(f"x  WARNING! Negative discriminant encountered, skipping Electrophilic descriptor calculation.")
    else:
        inter_phi = math.sqrt(discriminant)
        Phi = inter_phi - Eta
        Electrophilic_descriptor = round(((chi * (Phi / Gamma)) - (((Phi / Gamma) ** 2) * ((Eta / 2) + (Phi / 6)))), 4)

        # 5) Cubic electrophilicity index
    Gamma_cubic = 2 * delta_SCC_IP - Vertical_second_IP - delta_SCC_EA
    Eta_cubic = delta_SCC_IP - delta_SCC_EA

    if Eta_cubic != 0:
        Mu_cubic = (1 / 6) * ((-2 * delta_SCC_EA) - (5 * delta_SCC_IP) + Vertical_second_IP)
        w_cubic = round(((Mu_cubic ** 2) / (2 * Eta_cubic)) * (1 + ((Mu_cubic / (3 * (Eta_cubic) ** 2)) * Gamma_cubic)), 4)
    else:
        self.args.log.write(f"x  WARNING! Eta_cubic is zero, skipping cub. electrophilicity idx calculation.")

    # Return the calculated descriptors
    cdft_descriptors2 = {
        "Second IP": Vertical_second_IP,
        "Second EA": Vertical_second_EA,
        "Hyperhardness": hyper_hardness,
        "Hypersoftness": Global_hypersoftness,
        "Electrophilic descrip.": Electrophilic_descriptor,
        "cub. electrophilicity idx": w_cubic
    }

    return cdft_descriptors2

def calculate_local_CDFT_descriptors(file_fukui, cdft_descriptors, cdft_descriptors2,self):
    """
    Read fukui output file created from XTB and calculate local CDFT descriptors.
    """

    with open(file_fukui, "r") as f:
        data = f.readlines()

    f_pos, f_negs, f_rads = [], [], []
    start, end = None, None

    for i, line in enumerate(data):
        if "f(+) " in line:
            start = i + 1
        elif "-------------" in line and start is not None:
            end = i
            break

    if start is not None and end is not None:
        fukui_data = data[start:end]
        for line in fukui_data:
            try:
                f_po, f_neg, f_rad = map(lambda x: round(float(x), 4), line.split()[-3:])
                f_pos.append(f_po)
                f_negs.append(f_neg)
                f_rads.append(f_rad)
            except ValueError:
                continue
    else:
        self.args.log.write("WARNING: Fukui data not found in the file. Please check the '.fukui' file.")
        return None

    if not f_pos or not f_negs or not f_rads:
        self.args.log.write("WARNING: Fukui data lists are empty. Please check the '.fukui' file.")
        return None

    if None in [cdft_descriptors, cdft_descriptors2]:
        self.args.log.write("x  WARNING! Missing required CDFT descriptors (Softness, Hypersoftness, Electrophil. idx or Nucleophilicity idx).")
        return None

    chemical_softness = cdft_descriptors.get("Softness")
    Global_hypersoftness = cdft_descriptors2.get("Hypersoftness")
    electrophilicity_index = cdft_descriptors.get("Electrophil. idx")
    nucleophilicity_index = cdft_descriptors.get("Nucleophilicity idx")

    if None in [chemical_softness, Global_hypersoftness, electrophilicity_index, nucleophilicity_index]:
        self.args.log.write("x  WARNING! Missing required CDFT descriptors (Softness, Hypersoftness, Electrophil. idx or Nucleophilicity idx).")
        return None

    # Calculating local descriptors:
        # 1) dual descrip.
    dual_descriptor = [round(f_po - f_neg, 4) for f_po, f_neg in zip(f_pos, f_negs)]
        # 2) softness+, softness- and softness0
    s_pos = [round(chemical_softness * f_po, 4) for f_po in f_pos]
    s_negs = [round(chemical_softness * f_neg, 4) for f_neg in f_negs]
    s_rads = [round(chemical_softness * f_rad, 4) for f_rad in f_rads]
        # 3) Rel. nucleophilicity
    Relative_nucleophilicity = [round(s_neg / s_po, 4) if s_po != 0 else None for s_neg, s_po in zip(s_negs, s_pos)]
        # 4) Rel. electrophilicity
    Relative_electrophilicity = [round(s_po / s_neg, 4) if s_neg != 0 else None for s_neg, s_po in zip(s_negs, s_pos)]
        # 5) GC Dual Descrip.
    Grand_canonical_dual_descriptor = [round(Global_hypersoftness * dual, 4) for dual in dual_descriptor]
        # 6) softness+, softness- and softness0
    w_pos = [round(electrophilicity_index * f_po, 4) for f_po in f_pos]
    w_negs = [round(electrophilicity_index * f_neg, 4) for f_negs in f_negs]
    w_rads = [round(electrophilicity_index * f_rad, 4) for f_rad in f_rads]
        # 7) softness+, softness- and softness0
    Multiphilic_descriptor = [round(electrophilicity_index * dual, 4) for dual in dual_descriptor]
        # 8) softness+, softness- and softness0
    Nu_pos = [round(nucleophilicity_index * f_po, 4) for f_po in f_pos]
    Nu_negs = [round(nucleophilicity_index * f_neg, 4) for f_neg in f_negs]
    Nu_rads = [round(nucleophilicity_index * f_rad, 4) for f_rad in f_rads]

    localDescriptors = {
        "fukui+": f_pos,
        "fukui-": f_negs,
        "fukui0": f_rads,
        "dual descrip.": dual_descriptor,
        "softness+": s_pos,  # s+
        "softness-": s_negs,  # s-
        "softness0": s_rads,  # srad
        "Rel. nucleophilicity": Relative_nucleophilicity,  # s+/s-
        "Rel. electrophilicity": Relative_electrophilicity,  # s-/s+
        "GC Dual Descrip.": Grand_canonical_dual_descriptor,
        "Electrophil.": w_pos,  # w+
        "Nucleophil.": w_negs,  # w-
        "Radical attack": w_rads,  # wrad
        "Mult. descrip.": Multiphilic_descriptor,
        "Nu_Electrophil.": Nu_pos,  # Nu+
        "Nu_Nucleophil.": Nu_negs,  # Nu-
        "Nu_Radical attack": Nu_rads  # Nurad
    }

    return localDescriptors

def read_xtb(file,self):
    """
    Read xtb.out file and return a dictionary of extracted properties.
    """

    # Check if the file exists
    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    with open(file, "r") as f:
        data = f.readlines()

    # Initialize variables
    energy, homo_lumo, homo, lumo = np.nan, np.nan, np.nan, np.nan
    dipole_module, Fermi_level, transition_dipole_moment = np.nan, np.nan, np.nan
    total_charge, total_SASA = np.nan, np.nan
    total_C6AA, total_C8AA, total_alpha = np.nan, np.nan, np.nan
    atoms, numbers, chrgs = [], [], []
    covCN, C6AA, alpha = [], [], []

    # Parsing file data
    for i, line in enumerate(data):
        if "SUMMARY" in line:
            energy = float(data[i + 2].split()[3])
        elif "total charge" in line:
            total_charge = int(float(data[i].split()[3]))
        elif "(HOMO)" in line:
            if data[i].split()[3] != "(HOMO)":
                homo = round(float(data[i].split()[3]), 4)
                homo_occ = round(float(data[i].split()[1]), 4)
            else:
                homo = round(float(data[i].split()[2]), 4)
                homo_occ = 0
        elif "(LUMO)" in line:
            if data[i].split()[3] != "(LUMO)":
                lumo = round(float(data[i].split()[3]), 4)
                lumo_occ = round(float(data[i].split()[1]), 4)
            else:
                lumo = round(float(data[i].split()[2]), 4)
                lumo_occ = 0
        elif "molecular dipole:" in line:
            dipole_module = float(data[i + 3].split()[-1])
        elif "transition dipole moment" in line:
            transition_dipole_moment = float(data[i + 2].split()[-1])
        elif "Fermi-level" in line:
            Fermi_level = float(data[i].split()[-2])

    homo_lumo = round(float(lumo - homo), 4)

    # Getting atomic properties related to charges, dispersion, etc.
    start, end = 0, 0
    for j in range(len(data)):
        if "#   Z          covCN" in data[j]:
            start = j + 1
            break
    for k in range(start, len(data)):
        if "Mol. " in data[k]:
            end = k - 1
            total_C6AA = float(data[k].split()[-1])
            total_C8AA = float(data[k + 1].split()[-1])
            total_alpha = float(data[k + 2].split()[-1])
            break

    chrg_data = data[start:end]
    for line in chrg_data:
        item = line.split()
        numbers.append(int(item[0]))
        atoms.append(item[2])
        covCN.append(float(item[3]))
        chrgs.append(float(item[4]))
        C6AA.append(float(item[5]))
        alpha.append(float(item[6]))

    properties_dict = {
        "Total energy": energy,
        "Total charge": total_charge,
        "HOMO-LUMO gap": homo_lumo,
        "HOMO": homo,
        "LUMO": lumo,
        "atoms": atoms,
        "numbers": numbers,
        "charges": chrgs, 
        "Dipole module": dipole_module,
        "Fermi-level": Fermi_level,
        "Trans. dipole moment": transition_dipole_moment,
        "Coord. numbers": covCN,
        "Disp. coeff. C6": C6AA,
        "Polariz. alpha": alpha,
        "HOMO occup.": homo_occ,
        "LUMO occup.": lumo_occ,
        "Total SASA": total_SASA, 
        "Total disp. C6": total_C6AA,
        "Total disp. C8": total_C8AA,
        "Total polariz. alpha": total_alpha, 
    }

    return properties_dict


def read_fod(file,self):
    """
    Read xtb.fod files. Return FOD-related properties.
    """

    # Check if the file exists
    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    # Try to open the file and read its contents
    with open(file, "r") as f:
        data = f.readlines()  # Read all lines from the file


    # Initialize variables for storing FOD-related data
    total_fod = None  # Will store the total FOD value
    start_fod, end_fod = None, None  # Will mark the start and end of the FOD data section

    # Look for the line containing 'Loewdin FODpop' to identify the start of the FOD data
    for j, line in enumerate(data):
        if "Loewdin FODpop" in line:
            try:
                # Extract the total FOD value from two lines above 'Loewdin FODpop'
                total_fod = float(data[j - 2].split()[-1])
                # Mark the start of FOD data, which is the next line after 'Loewdin FODpop'
                start_fod = j + 1
            except (IndexError, ValueError) as e:
                # Handle potential errors from accessing invalid indices or incorrect value types
                self.args.log.write(f"x  WARNING! Error extracting total FOD: {e}")
                return None
            break  # Stop the loop once 'Loewdin FODpop' is found

    # If the 'Loewdin FODpop' section was not found, return None
    if start_fod is None:
        self.args.log.write(f"x  WARNING! 'Loewdin FODpop' not found in the file {file}.")
        return None

    # Look for the line that contains 'Wiberg/Mayer' to mark the end of the FOD data section
    for k in range(start_fod, len(data)):
        if "Wiberg/Mayer" in data[k]:
            # The FOD data ends just before the 'Wiberg/Mayer' section
            end_fod = k - 1
            break

    # If the end of the FOD section is not found, return None
    if end_fod is None:
        self.args.log.write(f"x  WARNING! 'Wiberg/Mayer' section not found in the file {file}.")
        return None

    # Extract and process the lines between start_fod and end_fod, which contain the FOD data
    fod_data = data[start_fod:end_fod]

    # Initialize lists to store FOD-related properties
    fod, s_prop_fod, p_prop_fod, d_prop_fod = [], [], [], []

    # Loop through each line of the FOD data and extract the relevant properties
    for line in fod_data:
        try:
            item = line.split()  # Split the line into individual components
            fod.append(float(item[1]))  # Extract the FOD value (index 1)
            s_prop_fod.append(float(item[2]))  # Extract the s proportion (index 2)
            p_prop_fod.append(float(item[3]))  # Extract the p proportion (index 3)
            d_prop_fod.append(float(item[4]))  # Extract the d proportion (index 4)
        except (IndexError, ValueError) as e:
            # Handle cases where the line is not properly formatted or contains invalid values
            self.args.log.write(f"x  WARNING! Error processing FOD data in line '{line.strip()}': {e}")
            return None

    # Create a dictionary to hold all the extracted FOD properties
    properties_FOD = {
        "Total FOD": total_fod,  # The total FOD value extracted from above
        "FOD": fod,  # List of FOD values for individual atoms
        "FOD s proportion": s_prop_fod,  # List of s proportion values
        "FOD p proportion": p_prop_fod,  # List of p proportion values
        "FOD d proportion": d_prop_fod,  # List of d proportion values
    }

    # Return the dictionary containing all FOD-related properties
    return properties_FOD


def read_json(file):
    """
    Takes json files and parses data into pandas table. Returns data.
    """

    if file.find(".json") > -1:
        f = open(file, "r")  # Opening JSON file
        data = json.loads(f.read())  # read file
        f.close()
        return data
    else:
        pass


def read_solv(file_solv):
    '''
    Retrieve properties from the single-point in solvent
    '''
    
    with open(file_solv, "r") as f:
        data = f.readlines()

        # Get molecular properties related to solvation (in kcal/mol)
        hartree_to_kcal = 627.509
        g_solv, g_elec, g_sasa, g_hb, g_shift = np.nan,np.nan,np.nan,np.nan,np.nan
        for _,line in enumerate(data):
            if '-> Gsolv' in line:
                g_solv = float(line.split()[3])*hartree_to_kcal
            elif '-> Gelec' in line:
                g_elec = float(line.split()[3])*hartree_to_kcal
            elif '-> Gsasa' in line:
                g_sasa = float(line.split()[3])*hartree_to_kcal
            elif '-> Ghb' in line:
                g_hb = float(line.split()[3])*hartree_to_kcal
            elif '-> Gshift' in line:
                g_shift = float(line.split()[3])*hartree_to_kcal

        # getting atomic properties related to solvation
        born_rad, atom_sasa, h_bond = [], [], []
        start_solv, end_solv = 0, 0
        for j in range(len(data)):
            if "#   Z     Born rad" in data[j]:
                start_solv = j + 1
                break
        for k in range(start_solv, len(data)):
            if "total SASA " in data[k]:
                end_solv = k - 1
                break

        solv_data = data[start_solv:end_solv]
        for line in solv_data:
            item = line.split()
            born_rad.append(float(item[3]))
            atom_sasa.append(float(item[4]))
            try:
                h_bond.append(float(item[5]))
            except IndexError:
                h_bond.append(0.0)

    properties_dict = {
        "G solv. in H2O": g_solv,
        "G of H-bonds H2O": g_hb,
        "G solv. elec.": g_elec,
        "G solv. SASA": g_sasa,
        "G solv. shift": g_shift,
        "Born radii": born_rad, 
        "Atomic SASAs": atom_sasa,
        "H bond H2O": h_bond, 
    }

    return properties_dict


def calculate_global_morfeus_descriptors(final_xyz_path,self):
    """
    Calculate global descriptors using the MORFEUS package. Return them in a structured dictionary.
    """

    # Initialize the descriptor variables as None
    sasa_area_global = None
    disp_area_global = None
    disp_vol_global = None

    # Try to load the molecular structure from the XYZ file
    try:
        elements, coordinates = read_xyz(final_xyz_path)
        elements, coordinates = read_geometry(final_xyz_path)
    except Exception as e:
        # Handle the error if the molecule cannot be loaded
        self.args.log.write(f"x  WARNING! Error loading molecule from file {final_xyz_path}: {e}")
        return {
            "Global SASA": sasa_area_global,
            "Disp. Area": disp_area_global,
            "Disp. Vol.": disp_vol_global,
        }

    # Calculate Global SASA (Solvent Accessible Surface Area)
    try:
        sasa = SASA(elements, coordinates)
        sasa_area_global = round(sasa.area, 4)  # Assign the calculated SASA value
    except Exception as e:
        # Handle any errors during SASA calculation
        self.args.log.write(f"x  WARNING! Error calculating SASA from Morfeus: {e}")

    # Calculate Global Dispersion Area and Volume
    try:
        disp = Dispersion(elements, coordinates)
        disp_area_global = round(disp.area, 4)  # Assign the calculated dispersion area
        disp_vol_global = round(disp.volume, 4)  # Assign the calculated dispersion volume
    except Exception as e:
        # Handle any errors during Dispersion calculation
        self.args.log.write(f"x  WARNING! Error calculating Global Dispersion from Morfeus: {e}")

    # Create the final dictionary with the descriptors, even if some values are None
    global_properties_morfeus = {
        "Global SASA": sasa_area_global,
        "Disp. Area": disp_area_global,
        "Disp. Vol.": disp_vol_global,
    }

    return global_properties_morfeus

def calculate_local_morfeus_descriptors(final_xyz_path,self):
    """
    Calculate local descriptors using the MORFEUS package. Return them in a structured dictionary.
    """

    # Initialize variables as None
    local_sasa_areas = None
    local_buried_volumes = None
    local_cone_angles = None
    local_solid_angles = None
    local_Pyramidalization = None
    local_vol_Pyramidalization = None
    local_disp = None

    # Try to load the molecular structure from the XYZ file
    try:
        elements, coordinates = read_xyz(final_xyz_path)
        elements, coordinates = read_geometry(final_xyz_path)
    except Exception as e:
        # If there's an issue loading the molecule, log the error and return an empty dictionary
        self.args.log.write(f"x  WARNING! Error loading molecule from file {final_xyz_path}: {e}")
        return {
            "SASA": local_sasa_areas,
            "Buried volume": local_buried_volumes,
            "Cone angle": local_cone_angles,
            "Solid angle": local_solid_angles,
            "Pyramidaliz. P": local_Pyramidalization,
            "Pyramidaliz. Vol": local_vol_Pyramidalization,
            "Dispersion": local_disp
        }

    # Calculate local SASA
    try:
        sasa = SASA(elements, coordinates)
        local_sasa_areas = [round(area, 4) for area in sasa.atom_areas.values()]
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local SASA from Morfeus: {e}")

    # Local buried Volume
    try:
        local_buried_volumes = []
        for i in range(len(coordinates)):
            bv = BuriedVolume(elements, coordinates, i)
            buried_volume = round(bv.fraction_buried_volume, 4)
            local_buried_volumes.append(buried_volume)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Buried Volume from Morfeus: {e}")

    # Local ConeAngle
    try:
        local_cone_angles = []
        for i in range(len(coordinates)):
            try:
                cone_angle = ConeAngle(elements, coordinates, i)
                local_cone_angles.append(round(cone_angle.cone_angle, 4))
            except Exception as e:
                local_cone_angles.append(None)  # Append None if there is an error for specific atom
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Cone Angle from Morfeus: {e}")

    # Local Solid Angle
    try:
        local_solid_angles = []
        for i in range(len(coordinates)):
            try:
                solid_angle = SolidAngle(elements, coordinates, i)
                local_solid_angles.append(round(solid_angle.cone_angle, 4))
            except Exception as e:
                local_solid_angles.append(None)  # Append None if there is an error for specific atom
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Solid Angle from Morfeus: {e}")

    # Pyramidalization
    try:
        local_Pyramidalization, local_vol_Pyramidalization = [], []
        for i in range(len(coordinates)):
            pyr = Pyramidalization(coordinates, i)
            local_Pyramidalization.append(round(pyr.P, 4))
            local_vol_Pyramidalization.append(round(pyr.P_angle, 4))
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating Pyramidalization from Morfeus: {e}")

    # Local Dispersion
    try:
        disp = Dispersion(elements, coordinates)
        local_disp = [round(p_int, 4) for p_int in disp.atom_p_int.values()]
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating Local Dispersion from Morfeus: {e}")

    # Return the dictionary with local descriptors, even if some values are None
    local_properties_morfeus = {
        "SASA": local_sasa_areas,
        "Buried volume": local_buried_volumes,
        "Cone angle": local_cone_angles,
        "Solid angle": local_solid_angles,
        "Pyramidaliz. P": local_Pyramidalization,
        "Pyramidaliz. Vol": local_vol_Pyramidalization,
        "Dispersion": local_disp
    }

    return local_properties_morfeus


def get_descriptors(level):
    """
    Returns descriptors for a given level from XTB and Morfeus
    """
    descriptors = {
        'denovo': {
            'mol': ["HOMO-LUMO gap", "HOMO", "LUMO", "IP", "EA", "Dipole module", "Total charge", "Global SASA", "G solv. in H2O", "G of H-bonds H2O"],
            'atoms': ["cm5 charges", "Electrophil.", "Nucleophil.", "Radical attack", "SASA", "Buried volume", "Cone angle", "H bond H2O"]
        },
        'interpret': {
            'mol': ["Fermi-level", "Total polariz. alpha", "Total FOD", "Electrophil. idx", "Hardness", "Softness", "Electronegativity",
                    "Nucleophilicity idx", "Second IP", "Second EA", "Disp. Area", "Disp. Vol.", 
                    "HOMO occup.", "LUMO occup."],
            'atoms': ["s proportion", "p proportion", "d proportion", "Coord. numbers",
                      "Polariz. alpha", "FOD", "FOD s proportion", "FOD p proportion", "FOD d proportion",
                      "Solid angle", "Pyramidaliz. P", "Pyramidaliz. Vol", "Dispersion"]
        },
        'full': {
            'mol': ["Total energy", "Total disp. C6", "Total disp. C8", "Chem. potential", "Electrodon. power idx",
                    "Electroaccep. power idx", 
                    "Electrofugality", "Nucleofugality", "Intrinsic React. idx", "Net Electrophilicity", 
                    "Hyperhardness", "Hypersoftness", "Electrophilic descrip.", "cub. electrophilicity idx",
                    "G solv. elec.", "G solv. SASA", "G solv. shift"],
            'atoms': ["fukui+", "fukui-", "fukui0", "dual descrip.", "softness+", "softness-", "softness0", 
                      "Rel. nucleophilicity", "Rel. electrophilicity", "GC Dual Descrip.", "Mult. descrip.", 
                      "Nu_Electrophil.", "Nu_Nucleophil.", "Nu_Radical attack", "Disp. coeff. C6", "Born radii"]
        }
    }

    return descriptors.get(level, {})


def fix_cols_names(df):
    '''
    Set code_name and SMILES using the right format
    '''
    
    for col in df.columns:
        if col.lower() == 'smiles':
            df = df.rename(columns={col: 'SMILES'})
        if col.lower() == 'code_name':
            df = df.rename(columns={col: 'code_name'})
    
    return df


def remove_atom_descp(df,atom_props):
    '''
    Remove atomic descriptors from a dataframe if they were not specified
    '''
    
    for col in df.columns:
        if col in atom_props:
            df = df.drop([col], axis=1).reset_index(drop=True)

    return df


def load_file_formats():
    '''
    Load formats for xTB calculations (loaded in QDESCP and in pytest).
    '''
    
    file_formats = {'_opt.out': 'Optimization',
                    '.out': 'Single-point',
                    '.fod': 'FOD',
                    '.fukui': 'Fukui',
                    '.gfn1': 'GFN1',
                    '.Nminus1': 'Nminus1',
                    '.Nminus2': 'Nminus2',
                    '.Nplus1': 'Nplus1',
                    '.Nplus2': 'Nplus2',
                    '.wbo': 'WBO',
                    '.solv': 'Solvation in H2O'}

    return file_formats