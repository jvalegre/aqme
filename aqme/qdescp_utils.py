######################################################.
#        This file stores QDESCP functions           #
######################################################.

import json
import sys
import os
import re
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
from morfeus import SASA, Dispersion, BuriedVolume, ConeAngle, SolidAngle, Pyramidalization, read_xyz, read_geometry, XTB
from morfeus.data import HARTREE_TO_KCAL
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


def read_json(file):
    """
    Loads JSON content from a file and returns it as a Python dictionary.
    On Windows, if parsing fails due to invalid backslash escapes, tries to fix by doubling backslashes.
    Returns None if the file cannot be opened or parsed.
    """
    if file.endswith(".json"):
        try:
            with open(file, "r", encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
                try:
                    with open(file, "r", encoding='utf-8') as f:
                        content = f.read()
                        fixed_content = re.sub(r'\\(?![\\nt"rbfu/])', r'\\\\', content)
                        return json.loads(fixed_content)
                except Exception:
                    return None
        except Exception:
            return None
    else:
        return None


def get_matches_idx_n_prefix(self,smarts_targets,name_initial):
    """
    Locate the indices (ie. [1], [1,4], etc) of the atoms involved in the patterns and their atom prefixes (ie. C, CO_C, etc)
    """

    pattern_dict = {}
    if len(smarts_targets) > 0:
        # create mol from SMILES in SDF files generated by CSEARCH or from regular SDF files
        mol = get_mol_assign(name_initial)

        # find the target atoms or groups in the mol object
        for pattern in smarts_targets:
            if "'" in pattern or '"' in pattern:
                pattern = pattern.replace("'",'').replace('"','')
            matches, idx_set = get_atom_matches(self,pattern,mol)
            if len(matches) == 0:
                self.args.log.write(f"x  WARNING! SMARTS pattern {pattern} not found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
            elif matches[0] == -1:
                self.args.log.write(f"x  WARNING! Mapped atom {pattern} not found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
            elif len(matches) > 1:
                self.args.log.write(f"x  WARNING! More than one {pattern} patterns were found in {os.path.basename(name_initial)}, atomic descriptors will not be generated.")
            elif len(matches) == 1:
                # get atom types and sort them to keep the same atom order among different molecules
                n_types, sorted_indices = sort_atom_types(matches,mol)

                # Generate unique match names for each atom type in the functional group
                match_names = get_prefix_atom_props(sorted_indices,mol,pattern,smarts_targets,idx_set)


                pattern_dict[pattern] = {'sorted_indices': sorted_indices,
                                            'match_names': match_names,
                                            'n_types': n_types,
                }
    return pattern_dict


def calculate_morfeus_descriptors(final_xyz_path,self,charge,mult,smarts_targets,name_initial):
    """
    Calculate local descriptors using the MORFEUS package. Return them in a structured dictionary.
    """

    # Electronic Morfeus descriptors using XTB
    morfeus_data = {}

    self.args.log.write(f"\no  Running MORFEUS for {os.path.basename(final_xyz_path).replace('.xyz','')} with charge {charge} and multiplicity {mult}")
    n_unpaired = int(mult) - 1 if mult is not None else 0

    # Try to load the molecular structure from the XYZ file
    try:
        elements, coordinates = read_xyz(final_xyz_path)
        elements, coordinates = read_geometry(final_xyz_path)
    except Exception as e:
        # If there's an issue loading the molecule, log the error and return an empty dictionary
        self.args.log.write(f"x  WARNING! Error loading molecule from file {final_xyz_path}: {e}")

    # calculate MORFEUS global descriptors that do not come from xTB
    morfeus_data = morfeus_global_descps(self, elements, coordinates, morfeus_data)

    # calculate MORFEUS local descriptors that do not come from xTB
    morfeus_data = morfeus_local_descps(self, elements, coordinates, morfeus_data, smarts_targets, name_initial)

    # xTB calculations through MORFEUS (with 1 proc to be reproducible)
    # calculate PTB method for descriptors that support it
    morfeus_data = morfeus_ptb_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired)
    
    # calculate GFN2 method for descriptors not included in PTB
    morfeus_data,energy = morfeus_gfn2_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired)

    # calculate GFN2 method with ALPB H2O for descriptors related to solvation
    morfeus_data = morfeus_solv_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired)

    # calculate GFN2 method in triplet state to calculate S0 to T1 energy gaps
    if n_unpaired == 0:
        morfeus_data = morfeus_t1_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, energy)

    # for properties that failed
    for prop in morfeus_data:
        if morfeus_data[prop] == []:
            morfeus_data[prop] = [None]*len(morfeus_data["Buried volume"])

    return morfeus_data


def morfeus_global_descps(self, elements, coordinates, morfeus_data):
    """
    Calculate MORFEUS local descriptors that do not come from xTB
    """

    # Initialize the descriptor variables as None
    sasa_area_global = disp_area_global = disp_vol_global = None

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

    # Update morfeus_data
    morfeus_data.update({
        "SASA": sasa_area_global,
        "Dispersion area": disp_area_global,
        "Dispersion volume": disp_vol_global,
    })

    return morfeus_data


def morfeus_local_descps(self, elements, coordinates, morfeus_data,smarts_targets,name_initial):
    """
    Calculate MORFEUS local descriptors that do not come from xTB. In all cases,
    the list of values of a property (ie SASA) has the same length as the number of atoms, but the
    code only calculates the properties for those atoms of interest. For example, if only
    the first atom of a molecule was defined in --qdescp_atoms, the list of SASA values
    will be [2.3452, NaN, NaN, NaN, NaN, NaN...]
    """
    
    # initialize
    local_sasa_areas = local_buried_volumes = local_cone_angles = local_solid_angles = []
    local_Pyramidalization = local_vol_Pyramidalization = local_disp = []

    # check atoms to calculate
    atom_matches = []
    pattern_dict = get_matches_idx_n_prefix(self,smarts_targets,name_initial)
    if len(pattern_dict.keys()) > 0:
            for pattern in pattern_dict:
                atom_matches += pattern_dict[pattern]['sorted_indices']

    # Calculate local SASA
    try:
        sasa = SASA(elements, coordinates)
        local_sasa_areas = [round(area, 4) for area in sasa.atom_areas.values()]
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local SASA from Morfeus: {e}")

    # Local Dispersion
    try:
        disp = Dispersion(elements, coordinates)
        local_disp = [round(p_int, 4) for p_int in disp.atom_p_int.values()]
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating Local Dispersion from Morfeus: {e}")

    # Local buried Volume
    try:
        local_buried_volumes = []
        for i in range(1,len(coordinates)+1):
            if i-1 in atom_matches: # MORFEUS starts indices at 1
                bv = BuriedVolume(elements, coordinates, i)
                buried_volume = round(bv.fraction_buried_volume*100, 2)
                local_buried_volumes.append(buried_volume)
            else:
                local_buried_volumes.append(np.nan)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Buried Volume from Morfeus: {e}")

    # Local ConeAngle
    try:
        local_cone_angles = []
        for i in range(1,len(coordinates)+1):
            if i-1 in atom_matches: # MORFEUS starts indices at 1
                try:
                    cone_angle = ConeAngle(elements, coordinates, i)
                    local_cone_angles.append(round(cone_angle.cone_angle, 4))
                except Exception as e:
                    local_cone_angles.append(np.nan)  # Append None if there is an error for specific atom
            else:
                local_cone_angles.append(np.nan)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Cone Angle from Morfeus: {e}")

    # Local Solid Angle
    try:
        local_solid_angles = []
        for i in range(1,len(coordinates)+1):
            if i-1 in atom_matches: # MORFEUS starts indices at 1
                try:
                    solid_angle = SolidAngle(elements, coordinates, i)
                    local_solid_angles.append(round(solid_angle.cone_angle, 4))
                except Exception as e:
                    local_solid_angles.append(np.nan)  # Append None if there is an error for specific atom
            else:
                local_solid_angles.append(np.nan)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating local Solid Angle from Morfeus: {e}")

    # Pyramidalization
    try:
        local_Pyramidalization, local_vol_Pyramidalization = [], []
        for i in range(1,len(coordinates)+1):
            if i-1 in atom_matches: # MORFEUS starts indices at 1
                pyr = Pyramidalization(coordinates, i)
                local_Pyramidalization.append(round(pyr.P, 4))
                local_vol_Pyramidalization.append(round(pyr.P_angle, 4))
            else:
                local_Pyramidalization.append(np.nan)
                local_vol_Pyramidalization.append(np.nan)
    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating Pyramidalization from Morfeus: {e}")

    # Update morfeus_data
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


def morfeus_ptb_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired):
    """
    Calculate PTB method for descriptors that support it
    """
    
    # initialize defaults
    homo = lumo = homo_lumo_gap = dipole_moment = None
    charges = atom_dipole_vect = bond_orders = {}

    try:
        xtb_ptb = XTB(elements, coordinates, n_processes=0, charge=charge, n_unpaired=n_unpaired, solvent=None, method='ptb')
        homo = xtb_ptb.get_homo(unit="eV")
        lumo = xtb_ptb.get_lumo(unit="eV")
        homo_lumo_gap = xtb_ptb.get_homo_lumo_gap(unit="eV")
        dipole_moment = xtb_ptb.get_dipole_moment(unit="debye")
        atom_dipole_vect = xtb_ptb.get_atom_dipoles()
        atom_dipole_modules = [np.linalg.norm(atom_dipole_vect[i]) for i in atom_dipole_vect.keys()]
        charges = list(xtb_ptb.get_charges(model="Mulliken").values())
        bond_orders = list(xtb_ptb.get_bond_orders().values())

    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating extra Morfeus descriptors with PTB: {e}")

    # update morfeus_data
    morfeus_data.update({
        "HOMO": homo,
        "LUMO": lumo,
        "HOMO-LUMO gap": homo_lumo_gap,
        "Dipole module": dipole_moment,
        "Atom dipole moment": atom_dipole_modules,
        "Partial charge": charges,
        "Bond orders": bond_orders,
    })

    return morfeus_data


def morfeus_gfn2_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired):
    """
    Calculate GFN2 method for descriptors not included in PTB
    """
    
    # initialize defaults
    ip = ea = hardness = softness = chemical_potential = None
    polarizability = atom_polarizabilities = fod_pop = None
    fod = fukui_plus = fukui_minus = fukui_radical = None
    electrophilicity = nucleofugality = None
    electrofugality = None

    try:
        xtb2 = XTB(elements, coordinates, n_processes=0, charge=charge, n_unpaired=n_unpaired, solvent=getattr(self.args, "qdescp_solvent", None), method=2)
        ip = xtb2.get_ip(corrected=True)
        ea = xtb2.get_ea(corrected=True)
        electrophilicity = xtb2.get_global_descriptor("electrophilicity", corrected=True)
        nucleofugality = xtb2.get_global_descriptor("nucleofugality", corrected=True)
        electrofugality = xtb2.get_global_descriptor("electrofugality", corrected=True)
        hardness = xtb2.get_hardness()
        softness = xtb2.get_softness()
        chemical_potential = xtb2.get_chemical_potential()
        polarizability = xtb2.get_molecular_polarizability()
        atom_polarizabilities = list(xtb2.get_atom_polarizabilities().values())
        fod_pop = list(xtb2.get_fod_population().values())
        fod = xtb2.get_nfod()
        fukui_plus = list(xtb2.get_fukui("plus").values())
        fukui_minus = list(xtb2.get_fukui("minus").values())
        fukui_radical = list(xtb2.get_fukui("radical").values())
        fukui_dual = list(xtb2.get_fukui("dual").values())
        local_electrophil = list(xtb2.get_fukui("local_electrophilicity").values())
        energy = xtb2.get_energy()
        fermi_level = xtb2.get_fermi_level()
        covCN = list(xtb2.get_covcn().values())

    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating extra Morfeus descriptors with GFN2: {e}")

    # update final JSON
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
        "Atom electrophilicity" : local_electrophil,
        "Electrophilicity": electrophilicity,
        "Nucleofugality": nucleofugality,
        "Electrofugality": electrofugality,
    })

    return morfeus_data,energy


def morfeus_solv_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired):
    """
    Calculate GFN2 method with ALPB H2O for descriptors related to solvation
    """
    
    # initialize
    g_solv = g_solv_hb = atom_hb_terms = None

    try:
        xtb2_solv = XTB(elements, coordinates, n_processes=0, charge=charge, n_unpaired=n_unpaired, solvent="h2o", method=2)
        g_solv = round(xtb2_solv.get_solvation_energy() * HARTREE_TO_KCAL, 2)
        g_solv_hb = round(xtb2_solv.get_solvation_h_bond_correction() * HARTREE_TO_KCAL, 2)
        atom_hb_terms = [round(v * HARTREE_TO_KCAL,2) for v in xtb2_solv.get_atomic_h_bond_corrections().values()]

    except Exception as e:
        self.args.log.write(f"x  WARNING! Error calculating solvation descriptors in MORFEUS with GFN2: {e}")

    # update final JSON
    morfeus_data.update({
        "G solv. in H2O": g_solv,
        "G of H-bonds H2O": g_solv_hb,
        "Atom H bond H2O": atom_hb_terms,
    })

    return morfeus_data


def morfeus_t1_descps(self, elements, coordinates, morfeus_data, charge, n_unpaired, energy):
    """
    Calculate GFN2 method in triplet state to calculate S0 to T1 energy gaps
    """
    
    # get S0-T1 gap
    xtb2t1 = XTB(elements, coordinates, n_processes=0, charge=charge, n_unpaired=n_unpaired+2, solvent=getattr(self.args, "qdescp_solvent", None), method=2)
    energy_T1 = xtb2t1.get_energy()
    S0_T1_gap = round((energy_T1 - energy) * HARTREE_TO_KCAL,2) # in kcal/mol
    morfeus_data["S0-T1 gap"] = S0_T1_gap

    return morfeus_data


def full_level_boltz(descp_dict,json_files,energy,smarts_targets,full_json_data):
    '''
    Get all the Boltzmann weighted properties (full level)
    '''

    # Calculate Boltzmann weights
    boltz = get_boltz(energy)

    # Get weighted atomic properties
    atomic_props = False
    for i, prop in enumerate(descp_dict['atom_props']):
        avg_prop = np.nan
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
    
    return full_json_data


def get_descriptors(level):
    """
    Returns descriptors for a given level from xTB and Morfeus
    """
    descriptors = {
        'denovo': {
            'mol': ["HOMO-LUMO gap", "HOMO", "LUMO", "IP", "EA", "Dipole module", "SASA", "G solv. in H2O", "G of H-bonds H2O"],
            'atoms': ["Partial charge", "Atom SASA", "Buried volume", "Atom H bond H2O", "Fukui+", "Fukui-", "Atom electrophilicity",]
        },
        'interpret': {
            'mol': ["Fermi-level", "Polarizability", "Total FOD", "Hardness", "Softness", "Dispersion area", "Dispersion volume", "Chem. potential",
                    "Electrophilicity", "Electrofugality", "Nucleofugality", "S0-T1 gap"],
            'atoms': ["Fukui_rad", "Fukui dual", "Atom dipole moment", "Coord. numbers",
                      "Atom Polarizability", "Atom FOD", "Atom dispersion", "Pyramidalization", "Pyramidaliz. volume"]
        },
        'full': {
            'mol': [],
            'atoms': []
        }
    }

    return descriptors.get(level, {})


def find_level_names(df_level, level):
    """
    Select which descriptors go into each level (de novo and interpret)
    """

    # fixed descriptors
    descriptors_denovo = ['code_name', 'SMILES']

    # other descriptors, in the same order that they will be in the CSV databases!
    if level == 'denovo':
        atom_suffixes = get_descriptors('denovo')['atoms'] + get_descriptors('denovo')['mol']
    if level == 'interpret': # interpret descriptors
        atom_suffixes = get_descriptors('denovo')['atoms'] + get_descriptors('interpret')['atoms'] + get_descriptors('denovo')['mol'] + get_descriptors('interpret')['mol']

    # build list of columns to keep
    cols_to_keep = []

    # always keep the fixed descriptors if present
    cols_to_keep.extend([c for c in descriptors_denovo if c in df_level.columns])

    # match prefixed atom columns
    for col in df_level.columns:
        for suffix in atom_suffixes:
            if col.endswith(suffix):  # matches regardless of prefix
                cols_to_keep.append(col)
                break  # avoid duplicates if multiple suffixes match

    return cols_to_keep


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
    
    with open(name, "w", encoding='utf-8') as outfile:
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

    # Adding a maximum time limit to avoid very long calculations
    params = rdFMCS.MCSParameters()
    params.Timeout = 30  # seconds maximum
    mcs = rdFMCS.FindMCS(mol_list, parameters=params)
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
        try:
            smarts_targets.remove(pattern)
        except ValueError:
            pass

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