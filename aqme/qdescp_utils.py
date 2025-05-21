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
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Descriptors
import warnings
warnings.filterwarnings('ignore')
from morfeus import SASA, Dispersion, BuriedVolume, ConeAngle, SolidAngle, Pyramidalization, read_xyz, read_geometry
from aqme.utils import load_sdf, periodic_table

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
    with open(final_boltz_file, "w") as outfile:
        json.dump(full_json_data, outfile)

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
        if str(p).lower() == 'nan' or p is None:
            boltz_avg = np.nan
            break
        boltz_avg.append([number * weights[i] for number in p])
    if str(p).lower() == 'nan':
        boltz_res = np.nan
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

def get_rdkit_properties(self,full_json_data, mol):
    """
    Calculates RDKit molecular descriptors
    """

    try:
        #level: full
        descrs = Descriptors.CalcMolDescriptors(mol)
        for descr in descrs:
            if descrs[descr] != np.nan and str(descrs[descr]).lower() != 'nan':
                full_json_data[descr] = descrs[descr]

    # For older versions of RDKit 
    except AttributeError:
        self.args.log.write(f"x  WARNING! Install a newer version of RDKit to get all the RDKit descriptors in the databse with all the descriptors. You can use: 'pip install rdkit --upgrade'.")
        full_json_data["MolLogP"] = rdkit.Chem.Descriptors.MolLogP(mol)

    return full_json_data


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
    with open(file, "r", encoding='utf-8') as f:
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
            q_mull = round(float(item[-5]), 3)
            q_cm5 = round(float(item[-4]), 3)
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
        "mulliken charge": mulliken,
        "cm5 charge": cm5,
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

    with open(file, "r", encoding='utf-8') as f:
        data = f.readlines()

    bonds, wbos = [], []
    for line in data:
        item = line.split()  # Split the line into components
        bond = [int(item[0]), int(item[1])]  # Extract bond indices
        wbo = round(float(item[2]), 3)  # Extract and round the WBO value
        bonds.append(bond)  # Add bond to list
        wbos.append(wbo)  # Add WBO to list

    return bonds, wbos

def calculate_global_CDFT_descriptors(file, file_Nminus1, file_Nminus2, file_Nplus1, file_Nplus2,self):
    """
    Read .gfn1 output file created from xTB and calculate CDFT descriptors with FDA approximations part 2
    """
    corr_xtb = 4.8455  # correction from xTB

    def extract_scc_energy(lines, filename):
        """
        Extract SCC energy from the file. If not found, return None and log a warning.
        """
        for line in lines:
            if "SCC energy" in line:
                return float(line.split()[3])
        
        self.args.log.write(f"x  WARNING! Could not find SCC energy value in the file: {os.path.basename(filename)}. Some CDFT-based descriptors will be missing.")
        return None

    try:
        # Open and read files
        with open(file, "r", encoding='utf-8') as f:
            data = f.readlines()
        with open(file_Nminus1, "r", encoding='utf-8') as f1:
            data1 = f1.readlines()
        with open(file_Nminus2, "r", encoding='utf-8') as f2:
            data2 = f2.readlines()
        with open(file_Nplus1, "r", encoding='utf-8') as f3:
            data3 = f3.readlines()
        with open(file_Nplus2, "r", encoding='utf-8') as f4:
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

    # Convert SCC energies to Hartree
    if scc_energy is not None:
        scc_energy *= Hartree
    if scc_energy_Nminus1 is not None:
        scc_energy_Nminus1 *= Hartree
    if scc_energy_Nminus2 is not None:
        scc_energy_Nminus2 *= Hartree
    if scc_energy_Nplus1 is not None:
        scc_energy_Nplus1 *= Hartree
    if scc_energy_Nplus2 is not None:
        scc_energy_Nplus2 *= Hartree

    # Initialize variables
    delta_SCC_IP, delta_SCC_EA, electrophilicity_index = None, None, None
    chemical_hardness, chemical_softness = None, None
    chemical_potential, mulliken_electronegativity = None, None
    electrodonating_power_index, electroaccepting_power_index = None, None
    intrinsic_reactivity_index = None
    electrofugality, nucleofugality, nucleophilicity_index, net_electrophilicity = None, None, None, None
    Vertical_second_IP, Vertical_second_EA = None, None
    hyper_hardness, Global_hypersoftness = None, None
    Electrophilic_descriptor, w_cubic = None, None

    # Calculate Global CDFT descriptors
    if None not in [scc_energy_Nminus1,scc_energy]:
        delta_SCC_IP = round(((scc_energy_Nminus1 - corr_xtb) - scc_energy),4)
    if None not in [scc_energy_Nplus1,scc_energy]:
        delta_SCC_EA = round((scc_energy - (scc_energy_Nplus1 + corr_xtb)),4)
    if None not in [delta_SCC_IP,delta_SCC_EA]:
        chemical_hardness = round((delta_SCC_IP - delta_SCC_EA), 4)
        chemical_potential = round(-(delta_SCC_IP + delta_SCC_EA) / 2, 4)
        electrophilicity_index = (chemical_potential**2)/(2*chemical_hardness)
        mulliken_electronegativity = round(-chemical_potential, 4)
        electrofugality = round(-delta_SCC_EA + electrophilicity_index, 4)
        nucleofugality = round(delta_SCC_IP + electrophilicity_index, 4)
        electrodonating_maximum_electron_flow = round((-(chemical_potential/chemical_hardness)),4)
        electrodonating_chemical_potential = round(((1/4)*((-3*delta_SCC_IP) - delta_SCC_EA)),4)
        electrodonating_maximum_electron_flow = round((-(electrodonating_chemical_potential/chemical_hardness)),4)
    if None not in [scc_energy_Nminus2,scc_energy_Nminus1]:
        Vertical_second_IP = round((((scc_energy_Nminus2 - scc_energy_Nminus1) - corr_xtb)), 4)
    if None not in [scc_energy_Nplus1,scc_energy_Nplus2]:
        Vertical_second_EA = round((((scc_energy_Nplus1 - scc_energy_Nplus2) + corr_xtb)), 4)
    if None not in [delta_SCC_IP,delta_SCC_EA,Vertical_second_IP,Vertical_second_EA]:  
        hyper_hardness = round((-((0.5) * (delta_SCC_IP + delta_SCC_EA - Vertical_second_IP - Vertical_second_EA))), 4)

    if chemical_hardness is not None and chemical_hardness != 0:
        chemical_softness = round(1 / chemical_hardness, 4)
        electrodonating_power_index = round(((delta_SCC_IP + 3 * delta_SCC_EA)**2) / (8 * chemical_hardness), 4)
        electroaccepting_power_index = round(((3 * delta_SCC_IP + delta_SCC_EA)**2) / (8 * chemical_hardness), 4)
        intrinsic_reactivity_index = round((delta_SCC_IP + delta_SCC_EA) / chemical_hardness, 4)
        if hyper_hardness is not None:
            Global_hypersoftness = round((hyper_hardness / ((chemical_hardness) ** 3)), 4)

        if electroaccepting_power_index != 0:
            nucleophilicity_index = round(10 / electroaccepting_power_index, 4)

    if None not in [electrodonating_power_index,electroaccepting_power_index]:
        net_electrophilicity = round((electrodonating_power_index - electroaccepting_power_index), 4)

    # For electrophilic descriptor calculations
    if None not in [scc_energy_Nplus1,scc_energy,Vertical_second_IP,delta_SCC_IP]:
        A = ((scc_energy_Nplus1 - scc_energy) + corr_xtb)
        c = (Vertical_second_IP - (2 * delta_SCC_IP) + A) / ((2 * Vertical_second_IP) - delta_SCC_IP - A)
        a = -((delta_SCC_IP + A) / 2) + (((delta_SCC_IP - A) / 2) * c)
        b = ((delta_SCC_IP - A) / 2) - (((delta_SCC_IP + A) / 2) * c)
        Gamma = (-3 * c) * (b - (a * c))
        Eta = 2 * (b - (a * c))
        chi = -a
        Mu = a

        discriminant = Eta ** 2 - (2 * Gamma * Mu)  # Checking the square root
        if discriminant >= 0:
            inter_phi = math.sqrt(discriminant)
            Phi = inter_phi - Eta
            Electrophilic_descriptor = round(((chi * (Phi / Gamma)) - (((Phi / Gamma) ** 2) * ((Eta / 2) + (Phi / 6)))), 4)

    # For cubic electrophilicity index
    if None not in [delta_SCC_IP,Vertical_second_IP,delta_SCC_EA]:
        Gamma_cubic = 2 * delta_SCC_IP - Vertical_second_IP - delta_SCC_EA
        Eta_cubic = delta_SCC_IP - delta_SCC_EA

    if Eta_cubic != 0:
        Mu_cubic = (1 / 6) * ((-2 * delta_SCC_EA) - (5 * delta_SCC_IP) + Vertical_second_IP)
        w_cubic = round(((Mu_cubic ** 2) / (2 * Eta_cubic)) * (1 + ((Mu_cubic / (3 * (Eta_cubic) ** 2)) * Gamma_cubic)), 4)

    # Return the calculated descriptors
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
        "Net Electrophilicity": net_electrophilicity,
        "Second IP": Vertical_second_IP,
        "Second EA": Vertical_second_EA,
        "Hyperhardness": hyper_hardness,
        "Hypersoftness": Global_hypersoftness,
        "Electrophilic descrip.": Electrophilic_descriptor,
        "cub. electrophilicity idx": w_cubic,
        "max. electron flow": electrodonating_maximum_electron_flow,
        "Electrodon. Chem. potential": electrodonating_chemical_potential,
        "Electrodon. max. electron flow": electrodonating_maximum_electron_flow
    }

    return cdft_descriptors

def calculate_local_CDFT_descriptors(file_fukui, cdft_descriptors,self):
    """
    Read fukui output file created from XTB and calculate local CDFT descriptors.
    """

    with open(file_fukui, "r", encoding='utf-8') as f:
        data = f.readlines()

    # Initialize variables
    f_pos, f_negs, f_rads = [], [], []
    dual_descriptor,s_pos,s_negs,s_rads = None,None,None,None
    Relative_nucleophilicity,Relative_electrophilicity,Grand_canonical_dual_descriptor = None,None,None
    w_pos,w_negs,w_rads,Multiphilic_descriptor = None,None,None,None
    Nu_pos,Nu_negs,Nu_rads = None,None,None

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

    if f_pos == [] or f_negs == [] or f_rads == []:
        self.args.log.write("x  WARNING! Fukui indices did not generate, please check the '.fukui' file.")

    chemical_softness = cdft_descriptors.get("Softness")
    Global_hypersoftness = cdft_descriptors.get("Hypersoftness")
    electrophilicity_index = cdft_descriptors.get("Electrophil. idx")
    nucleophilicity_index = cdft_descriptors.get("Nucleophilicity idx")

    # Calculating local descriptors
    if chemical_softness is not None and f_pos != [] and f_negs != [] and f_rads != []:
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

    if Global_hypersoftness is not None and dual_descriptor is not None:
            # 5) GC Dual Descrip.
        Grand_canonical_dual_descriptor = [round(Global_hypersoftness * dual, 4) for dual in dual_descriptor]

    if electrophilicity_index is not None and f_pos != [] and f_negs != [] and f_rads != []:
            # 6) softness+, softness- and softness0
        w_pos = [round(electrophilicity_index * f_po, 4) for f_po in f_pos]
        w_negs = [round(electrophilicity_index * f_neg, 4) for f_neg in f_negs]
        w_rads = [round(electrophilicity_index * f_rad, 4) for f_rad in f_rads]
        if dual_descriptor is not None:
                # 7) softness+, softness- and softness0
            Multiphilic_descriptor = [round(electrophilicity_index * dual, 4) for dual in dual_descriptor]

    if nucleophilicity_index is not None and f_pos != [] and f_negs != [] and f_rads != []:
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

    with open(file, "r", encoding='utf-8') as f:
        data = f.readlines()

    # Initialize variables
    energy, homo_lumo, homo, lumo = np.nan, np.nan, np.nan, np.nan
    dipole_module, Fermi_level = np.nan, np.nan
    total_charge = np.nan
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
        chrgs.append(round(float(item[4]),3))
        C6AA.append(float(item[5]))
        alpha.append(float(item[6]))

    properties_dict = {
        "Total energy": energy,
        "Total charge": total_charge,
        "HOMO-LUMO gap_GFN": homo_lumo,
        "HOMO_GFN": homo,
        "LUMO_GFN": lumo,
        "atoms": atoms,
        "numbers": numbers,
        "Partial charge_GFN": chrgs, 
        "Dipole module_GFN": dipole_module,
        "Fermi-level": Fermi_level,
        "Coord. numbers": covCN,
        "Disp. coeff. C6": C6AA,
        "Polariz. alpha": alpha,
        "HOMO occup.": homo_occ,
        "LUMO occup.": lumo_occ,
        "Total disp. C6": total_C6AA,
        "Total disp. C8": total_C8AA,
        "Total polariz. alpha": total_alpha, 
    }

    return properties_dict


def read_ptb(file,self):
    """
    Read xtb.ptb file and return a dictionary of extracted properties.
    """

    # Check if the file exists
    if not os.path.exists(file):
        self.args.log.write(f"x  WARNING! The file {file} does not exist.")
        return None

    with open(file, "r", encoding='utf-8') as f:
        data = f.readlines()

    # Initialize variables
    homo_lumo, homo, lumo = np.nan, np.nan, np.nan
    dipole_module = np.nan
    atom_dipoles, chrgs = [], []

    # Parsing file data
    for i, line in enumerate(data):
        if "(HOMO)" in line:
            if data[i].split()[3] != "(HOMO)":
                homo = round(float(data[i].split()[3]), 4)
            else:
                homo = round(float(data[i].split()[2]), 4)
        elif "(LUMO)" in line:
            if data[i].split()[3] != "(LUMO)":
                lumo = round(float(data[i].split()[3]), 4)
            else:
                lumo = round(float(data[i].split()[2]), 4)
        elif "Total dipole moment" in line:
            dipole_module = float(data[i + 1].split()[-1])

    homo_lumo = round(float(lumo - homo), 4)

    ptb_json = str(os.path.dirname(file)) + "/xtbout_ptb.json"
    if os.path.exists(ptb_json):
        # this part fixes a bug in xTB v1.7.1 when creating the json files
        with open(ptb_json, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        # Remove empty lines at the end of the file
        with open(ptb_json, 'w') as file:
            for line in lines:
                if line.rstrip('\n') != ',':
                    file.write(line)

        json_data = read_json(ptb_json)
        chrgs = json_data['partial charges']
        for dip_vector in json_data['atomic dipole moments']:
            atom_dipoles.append(math.sqrt(sum(pow(element, 2) for element in np.array(dip_vector))))
        os.remove(ptb_json)

    properties_dict = {
        "HOMO-LUMO gap": homo_lumo,
        "HOMO": homo,
        "LUMO": lumo,
        "Partial charge": chrgs,
        "Dipole module": dipole_module,
        "Dipole moment": atom_dipoles,
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
    with open(file, "r", encoding='utf-8') as f:
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
    Loads JSON content from a file and returns it as a Python dictionary.
    This function replaces single backslashes with double backslashes before parsing.
    Returns None if the file cannot be opened or parsed.
    """
    
    if file.find(".json") > -1:
        try:
            with open(file, "r", encoding='utf-8') as f:
                content = f.read()
                fixed_content = content.replace('\\', '\\\\')
                try:
                    data = json.loads(fixed_content)
                    return data
                except Exception:
                    return None
        except Exception:
            return None
    else:
        return None


def read_solv(file_solv):
    '''
    Retrieve properties from the single-point in solvent
    '''
    
    with open(file_solv, "r", encoding='utf-8') as f:
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
        "H bond with H2O": h_bond, 
    }

    return properties_dict


def read_triplet(file_triplet,singlet_e):
    '''
    Retrieve properties from the single-point in triplet
    '''

    hartree_to_kcal = 627.509
    triplet_e, transition_dipole_moment, singlet_triplet_gap = np.nan,np.nan,np.nan

    if os.path.exists(file_triplet):
        with open(file_triplet, "r", encoding='utf-8') as f:
            data = f.readlines()

            # Get molecular properties related to solvation (in kcal/mol)
            for i,line in enumerate(data):
                if "transition dipole moment" in line:
                    transition_dipole_moment = float(data[i + 2].split()[-1])
                elif "SUMMARY" in line:
                    triplet_e = float(data[i + 2].split()[3])

        singlet_triplet_gap = (triplet_e-singlet_e)*hartree_to_kcal

    properties_triplet = {
        "S0-T1 gap": singlet_triplet_gap,
        "Trans. dipole moment": transition_dipole_moment
        }

    return properties_triplet


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
            "Dispersion area": disp_area_global,
            "Dispersion volume": disp_vol_global,
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
        "Dispersion area": disp_area_global,
        "Dispersion volume": disp_vol_global,
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
            "Pyramidalization": local_Pyramidalization,
            "Pyramidaliz. volume": local_vol_Pyramidalization,
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
        for i in range(1,len(coordinates)+1):
            bv = BuriedVolume(elements, coordinates, i)
            buried_volume = round(bv.fraction_buried_volume, 4)
            local_buried_volumes.append(buried_volume)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Buried Volume from Morfeus: {e}")

    # Local ConeAngle
    try:
        local_cone_angles = []
        for i in range(1,len(coordinates)+1):
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
        for i in range(1,len(coordinates)+1):
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
        for i in range(1,len(coordinates)+1):
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
        "Pyramidalization": local_Pyramidalization,
        "Pyramidaliz. volume": local_vol_Pyramidalization,
        "Dispersion": local_disp
    }

    # for properties that failed in all calculations
    for prop in local_properties_morfeus:
        if local_properties_morfeus[prop] == []:
            local_properties_morfeus[prop] = [None]*len(local_properties_morfeus["Buried volume"])

    return local_properties_morfeus


def full_level_boltz(descp_dict,json_files,energy,smarts_targets,full_json_data):
    '''
    Get all the Boltzmann weighted properties (full level)
    '''

    # Calculate Boltzmann weights
    boltz = get_boltz(energy)

    # Get weighted atomic properties
    atomic_props = False
    for i, prop in enumerate(descp_dict['atom_props']):
        if i == 0:  # Filter to ensure all molecules have atomic properties for qdescp_atoms option
            for json_file in json_files:
                json_data = read_json(json_file)
                for atom_prop in descp_dict['atom_props']:
                    if atom_prop in json_data:
                        atomic_props = True
                        break
        if atomic_props:
            try:
                prop_list = [read_json(json_file)[prop] for json_file in json_files]
                avg_prop = average_properties(boltz, prop_list)
            except KeyError:
                full_json_data[prop] = np.nan
            
        else:
            avg_prop = np.nan
        update_full_json_data(full_json_data, prop, avg_prop, smarts_targets)

    # Get weighted molecular properties from XTB
    for prop in descp_dict['mol_props']:
        try:
            prop_list = [read_json(json_file)[prop] for json_file in json_files]
            avg_prop = average_properties(boltz, prop_list, is_atom_prop=False)
            full_json_data[prop] = avg_prop
        except KeyError:
            full_json_data[prop] = np.nan
    
    return full_json_data,atomic_props


def get_descriptors(level):
    """
    Returns descriptors for a given level from XTB and Morfeus
    """
    descriptors = {
        'denovo': {
            'mol': ["HOMO-LUMO gap", "HOMO", "LUMO", "IP", "EA", "Dipole module", "Total charge", "Global SASA", "G solv. in H2O", "G of H-bonds H2O"],
            'atoms': ["Partial charge", "Dipole moment", "Electrophil.", "Nucleophil.", "Radical attack", "SASA", "Buried volume", "H bond with H2O"]
        },
        'interpret': {
            'mol': ["Fermi-level", "Total polariz. alpha", "Total FOD", "Hardness", "Softness", "Electronegativity",
                    "Electrophil. idx", "Nucleophilicity idx", "Second IP", "Second EA", "S0-T1 gap"],
            'atoms': ["s proportion", "p proportion", "d proportion", "Coord. numbers",
                      "Polariz. alpha", "FOD", "Dispersion", "Pyramidalization", "Pyramidaliz. volume"]
        },
        'full': {
            'mol': ["HOMO-LUMO gap_GFN", "HOMO_GFN", "LUMO_GFN", "Dipole module_GFN", "Total energy", "Total disp. C6", "Total disp. C8", "Dispersion area", "Dispersion volume", "Chem. potential", 
                    "Electrodon. power idx", "Electroaccep. power idx", "max. electron flow", "Electrodon. Chem. potential", "Electrodon. max. electron flow",
                    "Electrofugality", "Nucleofugality", "Intrinsic React. idx", "Net Electrophilicity", 
                    "Hyperhardness", "Hypersoftness", "Electrophilic descrip.", "cub. electrophilicity idx",
                    "G solv. elec.", "G solv. SASA", "G solv. shift", "Trans. dipole moment",
                    "HOMO occup.", "LUMO occup."],
            'atoms': ["Partial charge_GFN", "fukui+", "fukui-", "fukui0", "dual descrip.", "softness+", "softness-", "softness0", 
                      "Rel. nucleophilicity", "Rel. electrophilicity", "GC Dual Descrip.", "Mult. descrip.", 
                      "Nu_Electrophil.", "Nu_Nucleophil.", "Nu_Radical attack", "Disp. coeff. C6", "Born radii",
                      "Cone angle", "Solid angle", "FOD s proportion", "FOD p proportion", "FOD d proportion"]
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
        if col.lower() == 'charge': # to fix Charge
            df = df.rename(columns={col: 'charge'})
        if col.lower() == 'mult': # to fix Mult
            df = df.rename(columns={col: 'mult'})    
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
                    '.ptb': 'PTB',
                    '.out': 'Single-point',
                    '.fod': 'FOD',
                    '.fukui': 'Fukui',
                    '.gfn1': 'GFN1',
                    '.Nminus1': 'Nminus1',
                    '.Nminus2': 'Nminus2',
                    '.Nplus1': 'Nplus1',
                    '.Nplus2': 'Nplus2',
                    '.wbo': 'WBO',
                    '.solv': 'Solvation in H2O',
                    '.stgap': 'S0-T1 (or T1-S0) gap',                    
                    }

    return file_formats


def assign_prefix_atom_props(prefix_list,atom_props,interpret_atoms,denovo_atoms):
    '''
    Assign prefixes to atomic properties based on the atoms used
    '''

    atom_props = add_prefix(prefix_list, atom_props)
    interpret_atoms = add_prefix(prefix_list, interpret_atoms)
    denovo_atoms = add_prefix(prefix_list, denovo_atoms)

    return atom_props,interpret_atoms,denovo_atoms


def add_prefix(prefix_list, prop_list):
    '''
    Add prefix to each atomic property
    '''
    
    new_atom_props = []

    for prefix in prefix_list:
        new_atom_props.append([f'{prefix}{ele}' for ele in prop_list])

    for i,new_name in enumerate(new_atom_props):
        if i == 0:
            prop_list = new_name
        else:
            prop_list = prop_list + new_name
            
    return prop_list


def collect_descp_lists():
    '''
    Retrieve all the descriptors names used in the different levels of descriptor interpretation
    '''
    denovo_descriptors = get_descriptors('denovo')
    interpret_descriptors = get_descriptors('interpret')
    full_descriptors = get_descriptors('full')

    # Descriptors in every level
    denovo_mols = denovo_descriptors['mol']
    denovo_atoms = denovo_descriptors['atoms']

    interpret_mols = denovo_mols + interpret_descriptors['mol']
    interpret_atoms = denovo_atoms + interpret_descriptors['atoms']

    mol_props = interpret_mols + full_descriptors['mol'] 
    atom_props =  interpret_atoms + full_descriptors['atoms']

    descp_dict = {'denovo_mols': denovo_mols,
                  'denovo_atoms': denovo_atoms,
                  'qdescp_denovo_csv': "QDESCP_denovo_descriptors.csv",
                  'interpret_mols': interpret_mols,
                  'interpret_atoms': interpret_atoms,
                  'qdescp_interpret_csv': "QDESCP_interpret_descriptors.csv",
                  'mol_props': mol_props,
                  'atom_props': atom_props,
                  'qdescp_csv': "QDESCP_full_descriptors.csv",
                  }

    return descp_dict


def average_properties(boltz, prop_list, is_atom_prop=True):
    """Calculate average properties based on Boltzmann weights."""
    return average_prop_atom(boltz, prop_list) if is_atom_prop else average_prop_mol(boltz, prop_list)


def update_full_json_data(full_json_data, prop, avg_prop, smarts_targets):
    """Update full_json_data with averaged properties."""
    if len(smarts_targets) > 0 or np.isnan(avg_prop).any():
        full_json_data[prop] = avg_prop
    else:
        full_json_data[prop] = avg_prop.tolist()


def dict_to_json(name, dict_data):
    '''
    Saves a dictionary as a JSON file
    '''
    
    with open(name, "w") as outfile:
        json.dump(dict_data, outfile)


def get_mols_qdescp(qdescp_files):
    '''
    Obtaining mols from input files
    '''
    
    mol_list = []
    for file in qdescp_files:
        with open(file, "r", encoding='utf-8') as F:
            lines = F.readlines()
            smi_exist = False
            for i, line in enumerate(lines):
                if ">  <SMILES>" in line:
                    smi = lines[i + 1].split()[0]
                    mol_indiv = Chem.AddHs(Chem.MolFromSmiles(smi))
                    mol_list.append(mol_indiv)
                    smi_exist = True
                    break
            if not smi_exist:
                mols = load_sdf(file)
                mol_indiv = mols[0]
                mol_list.append(mol_indiv)
    
    return mol_list


def get_mol_assign(name_initial):
    '''
    Create mol from SMILES in SDF files generated by CSEARCH or from regular SDF files
    '''

    sdf_file = f'{name_initial}.sdf'
    with open(sdf_file, "r", encoding='utf-8') as F:
        lines = F.readlines()

    smi_exist = False
    for i, line in enumerate(lines):
        if ">  <SMILES>" in line:
            smi = lines[i + 1].split()[0]
            mol = Chem.AddHs(Chem.MolFromSmiles(smi))
            smi_exist = True
    if not smi_exist:
        mols = load_sdf(sdf_file)
        mol = mols[0]

    return mol


def auto_pattern(mol_list,smarts_targets):
    '''
    Obtaing SMARTS patterns from the input files automatically if no patterns are provided
    '''

    mcs = rdFMCS.FindMCS(mol_list)
    if mcs is not None:
        common_substructure = Chem.MolFromSmarts(mcs.smartsString)
        # Filter out non-metal atoms
        atom_types = []
        for atom in common_substructure.GetAtoms():
            atom_types.append(atom.GetSymbol())
        for atom_type in set(atom_types):
            if atom_types.count(atom_type) == 1:
                smarts_targets.append(f'{atom_type}')

    return smarts_targets


def remove_invalid_smarts(self,mol_list,smarts_targets):
    '''
    Delete a SMARTS pattern if it is not compatible with more than 75% of the sdf files
    '''

    patterns_remove,matches = [],[]
    for pattern in smarts_targets:
        if "'" in pattern or '"' in pattern:
            pattern = pattern.replace("'",'').replace('"','')
        num_matches = len(mol_list)
        for mol_indiv in mol_list:
            try:
                # we differentiate if the pattern is a number for mapped atom or we are looking for smarts pattern in the molecule
                if not str(pattern).isalpha() and str(pattern).isdigit():
                    # for non-mapped mols (i.e. SDF input files)
                    mol_idxs = [atom.GetAtomMapNum() for atom in mol_indiv.GetAtoms()]
                    if len(set(mol_idxs)) == 1 and mol_idxs[0] == 0:
                        for i,atom in enumerate(mol_indiv.GetAtoms()):
                            if i == int(pattern)-1: # atoms in SDF starts in index 1, but Python starts in idx 0
                                pattern_idx = int(atom.GetIdx())
                                matches = ((int(pattern_idx),),)
                    # for mapped SMILES
                    else:
                        for atom in mol_indiv.GetAtoms():
                            if atom.GetAtomMapNum() == int(pattern):
                                pattern_idx = int(atom.GetIdx())
                                matches = ((int(pattern_idx),),)
                else:
                    matches = mol_indiv.GetSubstructMatches(Chem.MolFromSmarts(pattern))
            except:
                try: # I tried to make this except more specific for Boost.Python.ArgumentError, but apparently it's not as simple as it looks
                    matches = mol_indiv.GetSubstructMatches(Chem.MolFromSmarts(f'[{pattern}]'))
                except:
                    if pattern in self.args.qdescp_atoms:
                        self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} was not specified correctly and its corresponding atomic descriptors will not be generated! Make sure the qdescp_atoms option uses this format: \"[C]\" for atoms, \"[C=N]\" for bonds, and so on.")
                    patterns_remove.append(pattern)
                    break
            if len(matches) != 1:
                num_matches -= 1
            if (num_matches / len(mol_list)) < 0.75:
                patterns_remove.append(pattern)
                if pattern in self.args.qdescp_atoms:
                    self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} is not present (or there is more than one match) in 75% or more of the molecules. Atomic descriptors will not be generated for this pattern.")
                break

    # remove invalid patterns
    for pattern in patterns_remove:
        smarts_targets.remove(pattern)

    # print patterns recognized automatically
    for smarts_target in smarts_targets:
        if smarts_target not in self.args.qdescp_atoms:
            self.args.log.write(f"\no  Pattern {smarts_target} shared in all input files. Using it for atomic descriptor calculations")

    return smarts_targets


def get_atom_matches(self,pattern,mol):
    '''
    Find the target atoms or groups in the mol object
    '''

    matches = []
    idx_set = None
    
    # we differentiate if it is a number for mapped atom or we are looking for smarts pattern in the molecule
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


def sort_atom_types(matches,mol):
    '''
    Get atom types and sort them to keep the same atom order among different molecules
    '''
    
    atom_indices = list(matches[0])
    atom_types = []
    for atom_idx in atom_indices:
        atom_types.append(mol.GetAtoms()[atom_idx].GetSymbol())

    n_types = len(set(atom_types))
    if n_types == 1:
        sorted_indices = sorted(atom_indices, key=lambda idx: len(mol.GetAtoms()[idx].GetNeighbors()))
    elif n_types > 1:
        sorted_indices = sorted(atom_indices, key=lambda idx: mol.GetAtoms()[idx].GetAtomicNum())

    return n_types, sorted_indices


def get_prefix_atom_props(sorted_indices,mol,pattern,smarts_targets,idx_set):
    '''
    Generate unique match names for each atom type in the functional group
    '''
    
    match_names = []
    atom_counters = {}

    # Separate atoms when functional groups are used
    for atom_idx in sorted_indices:
        atom_type = mol.GetAtoms()[atom_idx].GetSymbol()

        # Check if there's more than one atom of the same type
        n_atoms_of_type = sum(1 for idx in sorted_indices if mol.GetAtoms()[idx].GetSymbol() == atom_type)

        # Initialize counter if it doesn't exist for this atom type
        if atom_type not in atom_counters:
            atom_counters[atom_type] = 1
        # If the pattern is a single atom number, we include that number in the match name
        if str(pattern).isdigit():  # This means the pattern is an atom number
            match_name = f'Atom_{pattern}'
        else:
            # If it's a SMARTS pattern or more than one atom
            if len(smarts_targets) == 1 and pattern in periodic_table():
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


def update_atom_props_json(sorted_indices,match_names,atom_props,json_data,prefixes_atom_prop,pattern,n_types):
    '''
    Assign atomic descriptors to each identified atom and update database for final JSON file
    '''

    for atom_idx, match_name in zip(sorted_indices, match_names):
        idx_xtb = atom_idx
        for prop in atom_props:
            try:
                if json_data[prop] is not None:
                    json_data[f'{match_name}_{prop}'] = json_data[prop][idx_xtb]
                else:
                    json_data[f'{match_name}_{prop}'] = None
                if f'{match_name}_' not in prefixes_atom_prop:
                    prefixes_atom_prop.append(f'{match_name}_')
            except (KeyError,IndexError): # prevents missing values
                pass

    # Adding max and min values for functional groups with the same two atoms
    if len(match_names) > 1 and n_types == 1:
        for prop in atom_props:
            prop_values = []
            for prop_name in match_names:
                prop_values.append(json_data[f'{prop_name}_{prop}'])
            if None not in prop_values:
                json_data[f'{pattern}_max_{prop}'] = max(prop_values)
                json_data[f'{pattern}_min_{prop}'] = min(prop_values)
            else:
                json_data[f'{pattern}_max_{prop}'] = None
                json_data[f'{pattern}_min_{prop}'] = None
            if f'{pattern}_max_{prop}' not in prefixes_atom_prop:
                prefixes_atom_prop.append(f'{pattern}_max_')
            if f'{pattern}_min_{prop}' not in prefixes_atom_prop:
                prefixes_atom_prop.append(f'{pattern}_min_')

    return prefixes_atom_prop, json_data