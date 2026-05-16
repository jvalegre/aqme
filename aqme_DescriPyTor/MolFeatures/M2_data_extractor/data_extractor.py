import pandas as pd
import numpy as np
import os
import sys
import math
from enum import Enum
import igraph as ig
from typing import *
import warnings
from morfeus import Sterimol
import networkx as nx
import numpy.typing as npt
# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    # Now you can import from the parent directory
    from gaussian_handler import feather_file_handler, save_to_feather
    # from MolAlign.renumbering import batch_renumbering
    from utils import visualize
    from utils.help_functions import *
    # from utils.rdkit_utils import xyz_list_to_mols, visualize_mols, mols_to_fingerprint_df 
    from extractor_utils.sterimol_utils import *
    from extractor_utils.bond_angle_length_utils import *
    from extractor_utils.dipole_utils import *
    from extractor_utils.vibrations_utils import *
    
except:
    from .gaussian_handler import feather_file_handler, save_to_feather
    # from ..MolAlign.renumbering import batch_renumbering
    from ..utils import visualize
    from ..utils.help_functions import *
    # from ..utils.rdkit_utils import xyz_list_to_mols, visualize_mols, mols_to_fingerprint_df 
    from .extractor_utils.sterimol_utils import *
    from .extractor_utils.bond_angle_length_utils import *
    from .extractor_utils.dipole_utils import *
    from .extractor_utils.vibrations_utils import *
    

warnings.filterwarnings("ignore", category=RuntimeWarning)



def show_highly_correlated_pairs(df, corr_thresh=0.8):
    """
    Display all feature pairs with |correlation| above corr_thresh,
    and highlight those with perfect correlation (|corr|=1).
    """
    corr = df.corr()
    pairs = []
    features = corr.columns
    for i in range(len(features)):
        for j in range(i+1, len(features)):
            cval = corr.iloc[i, j]
            if abs(cval) > corr_thresh:
                is_perfect = abs(cval) == 1.0
                pairs.append((features[i], features[j], cval))
    # Sort: perfect correlations first, then by absolute correlation
    pairs_sorted = sorted(pairs, key=lambda x: -abs(x[2]))
    if pairs_sorted:
        import pandas as pd
        table = pd.DataFrame(pairs_sorted, columns=["Feature 1", "Feature 2", "Correlation"])
    
        return table
    else:
      
        return None


from typing import Dict, Iterable, Optional, Union, Literal


import numpy as np
from typing import Dict, Tuple

def _validate_and_build_maps(n: int, renumbering_dict: Dict[int, int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Validate a 1-based renumbering dictionary and build:
      - old_to_new (0-based): index i (old) → new index
      - new_to_old (0-based): index j (new) → old index

    Key improvements:
      • Removes self-maps.
      • Automatically detects and closes open renumbering chains.
      • Prevents duplicate targets (enforces bijection).
      • Fills missing indices as identity.
    """
    import numpy as np

    # --- unwrap nested dicts (if {mol_name: {..}}) ---
    if len(renumbering_dict) == 1 and isinstance(next(iter(renumbering_dict.values())), dict):
        mol_name, inner = next(iter(renumbering_dict.items()))
        print(f"[DEBUG] Nested renumbering detected for molecule: {mol_name}")
        renumbering_dict = inner

    print(f"\n[DEBUG] Raw renumbering_dict (1-based): {renumbering_dict}")

    # --- remove self-maps ---
    renumbering_dict = {k: v for k, v in renumbering_dict.items() if k != v}
    print(f"[DEBUG] Removed self-maps: {renumbering_dict}")

    # --- validate numeric range ---
    valid = {}
    for k, v in renumbering_dict.items():
        if not (isinstance(k, int) and isinstance(v, int)):
            print(f"[WARN] Skipping non-int map {k}->{v}")
            continue
        if not (1 <= k <= n and 1 <= v <= n):
            print(f"[WARN] Out-of-range map {k}->{v} (n={n})")
            continue
        valid[k] = v
    renumbering_dict = valid

    # --- detect and close open chains ---
    visited = set()
    closed_pairs = {}
    for start in list(renumbering_dict.keys()):
        if start in visited:
            continue
        chain = [start]
        current = start
        while current in renumbering_dict and renumbering_dict[current] not in chain:
            nxt = renumbering_dict[current]
            chain.append(nxt)
            current = nxt
        for node in chain:
            visited.add(node)
        if len(chain) > 1:
            head, tail = chain[0], chain[-1]
            if tail not in renumbering_dict and head != tail:
                print(f"[FIX] Closing open chain {chain} with {tail}->{head}")
                renumbering_dict[tail] = head
            # pairwise connect the chain symmetrically
            for i in range(len(chain) - 1):
                closed_pairs[chain[i]] = chain[i + 1]
            if chain[-1] not in closed_pairs:
                closed_pairs[chain[-1]] = chain[0]
    renumbering_dict = closed_pairs or renumbering_dict
    print(f"[DEBUG] After chain closure: {renumbering_dict}")

    # --- enforce bijection (unique targets) ---
    seen = set()
    clean = {}
    for k, v in renumbering_dict.items():
        if v in seen:
            print(f"[WARN] Duplicate target {v}, skipping {k}->{v}")
            continue
        seen.add(v)
        clean[k] = v
    renumbering_dict = clean
    print(f"[DEBUG] Cleaned to unique targets: {renumbering_dict}")

    # --- fill missing indices as identity ---
    for i in range(1, n + 1):
        if i not in renumbering_dict and i not in renumbering_dict.values():
            renumbering_dict[i] = i
    print(f"[DEBUG] Added identity for unused indices")

    # --- build arrays ---
    old_to_new = np.arange(n, dtype=int)
    for old1, new1 in renumbering_dict.items():
        old_to_new[old1 - 1] = new1 - 1

    # --- ensure bijection (no missing/duplicate targets) ---
    unique, counts = np.unique(old_to_new, return_counts=True)
    dupes = unique[counts > 1]
    missing = np.setdiff1d(np.arange(n), unique)
    if len(dupes) or len(missing):
        print(f"[WARN] Fixing non-bijective mapping → duplicates={dupes.tolist()} missing={missing.tolist()}")
        free = list(missing)
        for i, val in enumerate(old_to_new):
            if val in dupes and free:
                replacement = free.pop(0)
                print(f"[FIX] old_to_new[{i}]={val} replaced by {replacement}")
                old_to_new[i] = replacement

    # --- build inverse safely ---
    new_to_old = np.full(n, -1, dtype=int)
    for old0, new0 in enumerate(old_to_new):
        if new_to_old[new0] != -1:
            print(f"[WARN] Two atoms map to same new index {new0+1}; keeping first")
            continue
        new_to_old[new0] = old0
    mask = new_to_old < 0
    if np.any(mask):
        new_to_old[mask] = np.arange(n)[mask]

    print(f"[DEBUG] old_to_new (0-based): {old_to_new}")
    print(f"[DEBUG] new_to_old (0-based): {new_to_old}\n")

    return old_to_new, new_to_old



def _reindex_like(obj: Union[pd.Series, pd.DataFrame], order: np.ndarray) -> Union[pd.Series, pd.DataFrame]:
    """Reindex a Series/DataFrame (row-wise) by an integer order vector (0-based)."""
    if hasattr(obj, "iloc"):
        out = obj.iloc[order].copy()
        out.index = range(len(out))  # normalize to 0..n-1
        return out
    raise TypeError(f"Unsupported type for reindexing: {type(obj)}")


def _remap_bonds_df(bonds_df: pd.DataFrame, old_to_new: np.ndarray, assume_0_based: bool = True) -> pd.DataFrame:
    """
    Remap bonds (two integer columns) using old_to_new (0-based map).
    If your bonds are 1-based, set assume_0_based=False.
    """
    out = bonds_df.copy()
    cols = out.columns[:2]
    if not assume_0_based:
        # convert 1-based to 0-based, map, then back to 1-based
        out[cols] = out[cols].astype(int) - 1
        out[cols] = out[cols].applymap(lambda idx: int(old_to_new[idx]))
        out[cols] = out[cols] + 1
    else:
        out[cols] = out[cols].astype(int)
        out[cols] = out[cols].applymap(lambda idx: int(old_to_new[idx]))
    return out



class Molecule:
    def __init__(self, molecule_feather_filename, parameter_list=None, new_xyz_df=None , threshold: float = 1.82):
        """
        Initialize a Molecule object with structural and computational data.
        
        Args:
            molecule_feather_filename (str): Path to Feather file containing molecular data
            parameter_list (list, optional): Precomputed parameters list. Defaults to None.
            new_xyz_df (DataFrame, optional): Alternative coordinates DataFrame. Defaults to None.
        """
        self.threshold = threshold  # Distance threshold for connectivity
        # Initialize core properties first
        self._initialize_core_properties(molecule_feather_filename, new_xyz_df)
        
        # Handle parameters and file loading
        self._load_parameters(molecule_feather_filename, parameter_list)
        
        # Initialize derived structural properties
        self._initialize_structural_properties()
        
        # Initialize computational results
        self._initialize_computational_data()


    def _initialize_core_properties(self, filename, new_xyz_df):
        """Set basic molecule identification and path properties"""
        self.molecule_name = os.path.splitext(os.path.basename(filename))[0]
        self.molecule_path = os.path.abspath(os.path.dirname(filename))
        self.new_index_order = new_xyz_df.index.tolist() if new_xyz_df is not None else None

    def _load_parameters(self, filename, parameter_list):
        """Handle parameter loading from file or existing list"""
        if parameter_list is None:
            original_dir = os.getcwd()
            try:
                os.chdir(self.molecule_path)
                self.parameter_list = feather_file_handler(filename)
            finally:
                os.chdir(original_dir)
        else:
            self.parameter_list = parameter_list

    def _initialize_structural_properties(self):
        """Set up structural components"""
        # XYZ coordinates
        try:
            self.xyz_df = self.parameter_list[0]['standard_orientation_df']
            self.coordinates_array = np.array(self.xyz_df[['x', 'y', 'z']].astype(float))
            
            # Connectivity and atom typing
            self.bonds_df = extract_connectivity(self.xyz_df, threshold_distance=self.threshold)
            self.atype_list = nob_atype(self.xyz_df, self.bonds_df)
            
            
        except Exception as e:
            print(f"Warning: Error initializing structural properties for {getattr(self, 'molecule_name', 'Unknown')}: {e}")
            # Set default values for missing properties
            self.xyz_df = pd.DataFrame(columns=['atom', 'x', 'y', 'z'])
            self.coordinates_array = np.array([])
            self.bonds_df = pd.DataFrame(columns=['atom1', 'atom2'])
            self.atype_list = []

    def _initialize_computational_data(self):
        """Initialize computed properties from calculations"""
        try:
            self.gauss_dipole_df = self.parameter_list[0].get('dipole_df', pd.DataFrame())
            self.polarizability_df = self.parameter_list[0].get('pol_df', pd.DataFrame())
            
            # General information
            self.info_df = self.parameter_list[0].get('info_df', pd.DataFrame())
            
            # Charge and vibrational data
            self.charge_dict = self.parameter_list[2] if len(self.parameter_list) > 2 else {}
            self.vibration_dict = self.parameter_list[1] if len(self.parameter_list) > 1 else {}
            self.vibration_mode_dict = build_vibration_mode_dict(self.vibration_dict, self.info_df)
            self.dfs = self.vibration_dict_to_dfs()
            # print('xxx')
            # print(self.dfs)
            # Energy values
            self.energy_value = self.parameter_list[0].get('energy_df', pd.DataFrame()) if len(self.parameter_list) > 2 else pd.DataFrame()
        
        except Exception as e:
            print(f"Warning: Error initializing computational data for {getattr(self, 'molecule_name', 'Unknown')}: {e}")
            # Set default values for missing properties
            if not hasattr(self, 'gauss_dipole_df'): self.gauss_dipole_df = pd.DataFrame()
            if not hasattr(self, 'polarizability_df'): self.polarizability_df = pd.DataFrame()
            if not hasattr(self, 'charge_dict'): self.charge_dict = {}
            if not hasattr(self, 'vibration_dict'): self.vibration_dict = {}
            if not hasattr(self, 'info_df'): self.info_df = pd.DataFrame()
            if not hasattr(self, 'energy_value'): self.energy_value = pd.DataFrame()




    def renumber_atoms(
        self,
        renumbering_dict: Dict[int, int],
        *,
        rebuild_connectivity: bool = False,
        bonds_are_0_based: bool = True,
        verbose: bool = True,
        save_feather: bool = True, out_file: Optional[str] = None
    ) -> None:
        """
        Renumber atoms according to a 1-based mapping and update relevant data.
        - Robust validation (range, duplicates)
        - Vectorized reindexing for xyz_df and charge_dict entries
        - Vibration keys 'vibration_atom_X' remapped to new indices
        - Connectivity can be preserved by remapping bonds_df, or recomputed

        Parameters
        ----------
        renumbering_dict : dict
            1-based dictionary: {old_atom_num: new_atom_num}.
            Unspecified atoms keep their original numbers (identity).
        rebuild_connectivity : bool, default False
            If True: recompute connectivity via `extract_connectivity(self.xyz_df)`.
            If False: preserve original `self.bonds_df` by remapping indices.
        bonds_are_0_based : bool, default True
            Set False if your `self.bonds_df` stores atom numbers 1-based.
        verbose : bool, default True
            Print a short success message and summary.

        Notes
        -----
        - Stores the following helpers:
            self.old_to_new (np.ndarray, 0-based)
            self.new_to_old (np.ndarray, 0-based permutation)
            self.new_index_order (list[int], 1-based new numbers in old-order for backward compat)
        - Assumes self.xyz_df has columns ['x','y','z'] (floats) and one row per atom.
        """
        if not hasattr(self, "xyz_df"):
            raise AttributeError("self.xyz_df is required before renumbering.")

        n = len(self.xyz_df)
        print(f"\n[DEBUG] Starting renumbering for molecule '{self.molecule_name}' with {n} atoms. {renumbering_dict}")
        try:
            print('debug renumbering_dict:', renumbering_dict)
            old_to_new, new_to_old = _validate_and_build_maps(n, renumbering_dict)
        except Exception as e:
            # print error and return without making changes so it continues to work in batch mode
            print(f"[ERROR] Renumbering validation failed for molecule '{self.molecule_name}': {e}")
            return
        # Save maps for downstream use
        self.old_to_new = old_to_new              # old idx -> new idx (0-based)
        self.new_to_old = new_to_old              # new idx -> old idx (0-based)
        self.new_index_order = (old_to_new + 1).tolist()  # 1-based, for compatibility

        # --- Reorder XYZ (vectorized) ---
        old_xyz_df = self.xyz_df.copy()
        print(f"[DEBUG] old_xyz_df before renumbering:\n{old_xyz_df}")
        self.xyz_df = old_xyz_df.iloc[new_to_old].reset_index(drop=True)
        print(f"[DEBUG] new_xyz_df after renumbering:\n{self.xyz_df}")
        # Ensure numeric coordinate dtypes
        for c in ("x", "y", "z"):
            if c in self.xyz_df.columns:
                self.xyz_df[c] = pd.to_numeric(self.xyz_df[c], errors="coerce")
        self.coordinates_array = self.xyz_df[["x", "y", "z"]].to_numpy(dtype=float)

        # --- Connectivity & atom types ---
        if hasattr(self, "bonds_df") and isinstance(self.bonds_df, pd.DataFrame) and not rebuild_connectivity:
            # Preserve original bonds by remapping indices
            self.bonds_df = _remap_bonds_df(self.bonds_df, old_to_new=self.old_to_new, assume_0_based=False)
        else:
            # Recompute from geometry (may differ from original)
            if "extract_connectivity" not in globals() and not hasattr(self, "extract_connectivity"):
                raise NameError("extract_connectivity is not available to recompute bonds. Set rebuild_connectivity=False or provide it.")
            # Support bound method or global function
            _extract = getattr(self, "extract_connectivity", globals().get("extract_connectivity"))
            self.bonds_df = _extract(self.xyz_df)

        # Atom types (depends on your implementation)
        if "nob_atype" not in globals() and not hasattr(self, "nob_atype"):
            # Skip silently if not available
            self.atype_list = getattr(self, "atype_list", None)
        else:
            _atype = getattr(self, "nob_atype", globals().get("nob_atype"))
            self.atype_list = _atype(self.xyz_df, self.bonds_df)

        # --- Renumber charges (vectorized) ---
        if hasattr(self, "charge_dict") and self.charge_dict:
            new_charge_dict: Dict[str, Any] = {}
            for charge_type, charge_data in self.charge_dict.items():
                if isinstance(charge_data, dict):
                    new_sub = {}
                    for subkey, subdata in charge_data.items():
                        if hasattr(subdata, "iloc"):
                            new_sub[subkey] = _reindex_like(subdata, new_to_old)
                        else:
                            raise TypeError(f"Unsupported nested type in charge_dict[{charge_type}][{subkey}]: {type(subdata)}")
                    new_charge_dict[charge_type] = new_sub
                elif hasattr(charge_data, "iloc"):
                    new_charge_dict[charge_type] = _reindex_like(charge_data, new_to_old)
                else:
                    raise TypeError(f"Unsupported type in charge_dict[{charge_type}]: {type(charge_data)}")
            self.charge_dict = new_charge_dict

        # --- Renumber vibration dict keys like 'vibration_atom_{i}' ---
        if hasattr(self, "vibration_dict") and self.vibration_dict:
            print(f"[DEBUG] Starting vibration_dict renumbering with {len(self.vibration_dict)} entries")
            print(f"[DEBUG] Number of atoms (n): {n}")
            print(f"[DEBUG] old_to_new : {list(old_to_new) if hasattr(old_to_new, '__len__') else old_to_new}")

            new_vib_dict: Dict[str, Any] = {}
            for key, val in self.vibration_dict.items():
                if key.startswith("vibration_atom_"):
                    try:
                        old_atom_1b = int(key.split("_")[-1])
                        print(f"[DEBUG] Key '{key}' parsed old_atom_1b={old_atom_1b}")
                    except ValueError:
                        print(f"[WARN] Key '{key}' has unexpected format, keeping as-is")
                        new_vib_dict[key] = val
                        continue

                    # Range check
                    if not (1 <= old_atom_1b <= n):
                        print(f"[WARN] old_atom_1b={old_atom_1b} out of valid range 1–{n}, keeping as-is")
                        new_vib_dict[key] = val
                        continue

                    # Map 1-based old → 0-based → 1-based new
                    try:
                        new_atom_1b = int(old_to_new[old_atom_1b - 1] + 1)
                    except Exception as e:
                        print(f"[ERROR] Failed mapping old_atom_1b={old_atom_1b}: {e}")
                        new_vib_dict[key] = val
                        continue

                    new_key = f"vibration_atom_{new_atom_1b}"
                    print(f"[DEBUG] Mapping {key} → {new_key}")

                    if new_key in new_vib_dict:
                        print(f"[ERROR] Collision detected! {key} and another map to {new_key}")
                        raise ValueError(f"Collision in vibration_dict keys after renumbering: {new_key}")

                    new_vib_dict[new_key] = val
                else:
                    # Non-vibration_atom_ keys copied as-is
                    new_vib_dict[key] = val
                    print(f"[DEBUG] Key '{key}' copied as-is (not vibration_atom_)")

            print(f"[DEBUG] Renumbering finished: {len(new_vib_dict)} keys total "
                f"({len(self.vibration_dict)} before)")
            self.vibration_dict = new_vib_dict

        # Note: self.vibration_mode_dict may need similar treatment if it uses atom indices in keys
        if hasattr(self, "vibration_mode_dict") and self.vibration_mode_dict:
            new_mode_dict: Dict[str, Any] = {}
            for key, val in self.vibration_mode_dict.items():
                # key is frequency, np.float value is array of shape (n, 3)
                if isinstance(val, np.ndarray) and val.shape[0] == n and val.shape[1] == 3:
                    new_array = val[old_to_new, :]
                    new_mode_dict[key] = new_array
                else:
                    new_mode_dict[key] = val  # leave untouched if unexpected format
                    continue
                    
            self.vibration_mode_dict = new_mode_dict
        if verbose:
            moved = {k: v for k, v in renumbering_dict.items() if k != v}
            msg = "Atoms successfully renumbered."
            if moved:
                msg += f" Moved: {moved} (1-based)."
            print(msg)
        
        if save_feather:
            if out_file is None:
                ren_dir = os.path.join(self.molecule_path, 'renumbered')
                os.makedirs(ren_dir,  exist_ok=True)  # just os.makedirs(ren_dir, exist_ok=True) is fine
                out_file = os.path.join(ren_dir, f"{self.molecule_name}_renumbered.feather")
            save_to_feather(
                df_dict={
                    "xyz": self.xyz_df,
                    "dipole": self.gauss_dipole_df,
                    "polarizability": self.polarizability_df,
                    "info": self.info_df
                },
                vib_dict=self.vibration_dict,
                charge_dict=self.charge_dict,
                ev_df=self.energy_value,
                out_file=out_file
            )
            if verbose:
                print(f"Renumbered data saved to {out_file}\n Saved to path: {os.getcwd()}")



    def get_molecule_descriptors_df(self) -> pd.DataFrame:
        """
        Returns a DataFrame with various descriptors calculated for the molecule.

        Returns:
            pd.DataFrame: A DataFrame with the calculated descriptors.
        """
        descriptors_df = pd.concat(
            [self.energy_value, self.gauss_dipole_df, self.get_sterimol()], axis=1)
        descriptors_df = descriptors_df.sort_values(by=['energy_value'])
        return descriptors_df


    def write_xyz_file(self) -> None:
        """
        Writes the molecule's coordinates to an XYZ file.
        """
        data_to_xyz(self.xyz_df, self.molecule_name+'.xyz')
    
    def write_csv_files(self) -> None:
        """
        Writes all class variables to CSV files with the variable name as the file 
        """
        for var_name, var_value in vars(self).items():
            if isinstance(var_value, pd.DataFrame):
                var_value.to_csv(f"{var_name}.csv", index=False)


    def vibration_dict_to_dfs(self, debug: bool = False):
        """
        Converts the vibration dictionary to a DataFrame for easier handling.

        Args:
            debug (bool): If True, prints debugging information.

        Returns:
            list[pd.DataFrame]: A list of DataFrames grouped by frequency.
        """
        
        self.vib_df_list = []
        self.bonds_df_vib = self.bonds_df

        if debug:
            print(f"[DEBUG] Starting vibration_dict_to_dfs")
            print(f"[DEBUG] vibration_dict contains {len(self.vibration_dict)} entries")

        for key, value in self.vibration_dict.items():
            if debug:
                print(f"\n[DEBUG] Processing key: {key}")
                print(f"[DEBUG] Raw value shape: {value.shape}")

            # --- get atom index from key ---
            idx = int(key.split('_')[-1]) - 1

            # --- filter out NaNs ---
            before_shape = value.shape
            value = value[~np.isnan(value).any(axis=1)]
            after_shape = value.shape
            if debug:
                print(f"[DEBUG] Removed NaNs: {before_shape} -> {after_shape}")

            # --- construct new array ---
            new_array = np.array([value[i] for i in range(value.shape[0])])
            if debug:
                print(f"[DEBUG] new_array shape: {new_array.shape}")

            # --- build DataFrame ---
            df = pd.DataFrame(new_array, columns=['x', 'y', 'z'])
            df['frequency'] = self.info_df['Frequency'].values

            indice = [self.xyz_df['atom'].iloc[idx]] * len(value)
            df.insert(0, 'atom', indice)

            if debug:
                print(f"[DEBUG] DataFrame built for {key}: shape {df.shape}")
                print(df.head())

            # --- update attributes ---
            self.vib_df_list.append(df)
            self.vibration_dict[key] = new_array

            # --- rebuild combined xyz_df ---
            xyz_df = pd.DataFrame()
            for vib_df in self.vib_df_list:
                xyz_df = pd.concat([xyz_df, vib_df[['atom','x','y','z','frequency']]], axis=0)

            self.dfs = [group.reset_index(drop=True) for _, group in xyz_df.groupby(xyz_df.iloc[:, -1])]

        if debug:
            print(f"\n[DEBUG] Finished building dfs: {len(self.dfs)} groups")
            for i, df in enumerate(self.dfs[:3]):  # preview only first 3 groups
                print(f"[DEBUG] Group {i} shape: {df.shape}")
                print(df.head())

        return self.dfs


    def visualize_molecule(self, vector=None) -> None:
        """
        Visualizes the molecule using the `visualize` module.
        """
        if vector is not None:
            visualize.show_single_molecule(molecule_name=self.molecule_name, xyz_df=self.xyz_df, dipole_df=vector,origin=[0,0,0])
        else:
            visualize.show_single_molecule(molecule_name=self.molecule_name, xyz_df=self.xyz_df, dipole_df=self.gauss_dipole_df,origin=[0,0,0])


    def process_sterimol_atom_group(self, atoms, radii, sub_structure=True, drop_atoms=None,visualize_bool=None, mode='all') -> pd.DataFrame:

        # connected = get_molecule_connections(self.bonds_df, atoms[0], atoms[1], mode=mode)
        connected = None

        return get_sterimol_df(self.xyz_df, self.bonds_df, atoms, connected, radii, sub_structure=sub_structure, drop_atoms=drop_atoms, visualize_bool=visualize_bool, mode=mode)

    def get_sterimol(self, base_atoms: Union[None, Tuple[int, int]] = None, radii: str = 'CPK', sub_structure=True, drop_atoms=None, visualize_bool=None, mode='all') -> pd.DataFrame:
        """
        Returns a DataFrame with the Sterimol parameters calculated based on the specified base atoms and radii.

        Args:
            base_atoms (Union[None, Tuple[int, int]], optional): The indices of the base atoms to use for the Sterimol calculation. Defaults to None.
            radii (str, optional): The radii to use for the Sterimol calculation. Defaults to 'bondi'.

        Returns:
            pd.DataFrame: A DataFrame with the Sterimol parameters.
            
            to add
            - only_sub- sterimol of only one part.
            - drop some atoms.
        """
        if base_atoms is None:
            base_atoms = get_sterimol_indices(self.xyz_df, self.bonds_df)
        
        if isinstance(base_atoms[0], list):
            # If base_atoms is a list of lists, process each group individually and concatenate the results
            sterimol_list = [self.process_sterimol_atom_group(atoms, radii, sub_structure=sub_structure, drop_atoms=drop_atoms,visualize_bool=visualize_bool, mode=mode) for atoms in base_atoms]
            sterimol_df = pd.concat(sterimol_list, axis=0)

        else:
            # If base_atoms is a single group, just process that group
            sterimol_df = self.process_sterimol_atom_group(base_atoms, radii,sub_structure=sub_structure, drop_atoms=drop_atoms,visualize_bool=visualize_bool, mode=mode)
        return sterimol_df


    def swap_atom_pair(self, pair_indices: Tuple[int, int]) -> pd.DataFrame:
        """
        Swaps the positions of two atoms in the molecule and returns a new DataFrame with the updated coordinates.

        Args:
            pair_indices (Tuple[int, int]): The indices of the atoms to swap.

        Returns:
            pd.DataFrame: A new DataFrame with the updated coordinates.
        """
        pairs = adjust_indices(pair_indices)
        xyz_df = self.xyz_df
        temp = xyz_df.iloc[pairs[0]].copy()
        xyz_df.iloc[pairs[0]] = self.coordinates_array[pairs[1]]
        xyz_df.iloc[pairs[1]] = temp
        return xyz_df

    
    def get_charge_df(self, atoms_indices: List[int], 
                  type: Union[str, List[str]] = 'all'
                 ) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
        """
        Returns a DataFrame with the charges for the specified atoms for the requested type(s).

        Args:
            atoms_indices (List[int]): Indices of atoms to include.
            type (str or List[str], optional): Charge type to extract (e.g., 'nbo', 'hirshfeld', 'cm5').
                                            Use 'all' or a list of types to extract multiple types.
                                            Defaults to 'nbo'.

        Returns:
            pd.DataFrame or Dict[str, pd.DataFrame]: A DataFrame for a single type, or a dictionary 
            of DataFrames if multiple types are requested. If 'all' is specified, any missing/invalid 
            types are skipped instead of raising an error.
        """
        atoms_indices = adjust_indices(atoms_indices)

        def extract_charge(df: pd.DataFrame) -> pd.DataFrame:
            # Extract specified rows, rename indices, and transpose
            return df.iloc[atoms_indices].rename(index=lambda x: f'atom_{x + 1}').T

        # Handle different input cases for 'type'
        if isinstance(type, str):
            if type.lower() == 'all':
                # Instead of raising an error on missing types,
                # we iterate through every type in self.charge_dict
                # and skip any that fail.
                results = {}
                for t, df in self.charge_dict.items():
                    try:
                        results[t] = extract_charge(df)
                    except Exception as e:
                        pass
                return results

            else:
                # Single specific type

                return extract_charge(self.charge_dict[type])

        elif isinstance(type, list):
            # Multiple types, specified by name
            result = {}
            for t in type:
         
                result[t] = extract_charge(self.charge_dict[t])
            return result

        else:
            raise TypeError("Parameter 'type' must be a string or a list of strings.")


    def get_charge_diff_df(self, diff_indices: Union[List[List[int]], List[int]], 
                       type: Union[str, List[str]] = 'all'
                      ) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
        """
        Returns a DataFrame with the charge differences for specified atom pairs for the requested type(s).

        Args:
            diff_indices (List[List[int]] or List[int]):
                - For multiple pairs: a list of lists (e.g., [[atom1, atom2], [atom3, atom4]]).
                - For a single pair: a list (e.g., [atom1, atom2]).
            type (str or List[str], optional): Charge type to extract (e.g., 'nbo', 'hirshfeld', 'cm5').
                                            Use 'all' or a list of types to extract multiple types.
                                            Defaults to 'nbo'.

        Returns:
            pd.DataFrame or Dict[str, pd.DataFrame]: A DataFrame for a single type, or a dictionary of DataFrames if multiple types are requested.
            When using 'all' (or a list), any type that fails to compute will be skipped.
        """
        def compute_diff(df: pd.DataFrame) -> pd.DataFrame:
            # If the first element of diff_indices is a list, assume multiple pairs are provided.
            try:
                
                if isinstance(diff_indices[0], list):
                    diff_list = []
                    for atoms in diff_indices:
                        atoms = adjust_indices(atoms)
                        # Calculate difference between the two specified atoms.
                        diff = pd.DataFrame(
                            [df.iloc[atoms[0]].values - df.iloc[atoms[1]].values],
                            columns=[f'diff_{atoms[0]+1}-{atoms[1]+1}']
                        )
                        diff_list.append(diff)
                        # print(f"Difference for atoms {atoms[0]+1} and {atoms[1]+1}: {diff.values}")
                    return pd.concat(diff_list, axis=1)
                else:
                    diff_indices_adj = adjust_indices(diff_indices)
                    diff_df=pd.DataFrame(
                        [df.iloc[diff_indices_adj[0]].values - df.iloc[diff_indices_adj[1]].values],
                        columns=[f'diff_{diff_indices_adj[0]+1}-{diff_indices_adj[1]+1}']
                    )
        
                    return diff_df  
            except Exception as e:
                raise ValueError(f"Error in computing charge difference: {e}")

        # Handle different input cases for 'type'
        if isinstance(type, str):
            if type.lower() == 'all':
                results = {}
                for t, df in self.charge_dict.items():
                    try:
                        results[t] = compute_diff(df)
                    except Exception as e:
                        pass

                return results
            else:

                return compute_diff(self.charge_dict[type])
        elif isinstance(type, list):
            result = {}
            for t in type:

                try:
                    result[t] = compute_diff(self.charge_dict[t])
                except Exception as e:
                    pass

            return result
        else:
            raise TypeError("Parameter 'type' must be a string or a list of strings.")

    def get_coordinates_mean_point(self, atoms_indices: List[int]) -> np.ndarray:
        """
        Returns the mean point of the specified atoms.

        Args:
            atoms_indices (List[int]): The indices of the atoms to calculate the mean point for.

        Returns:
            np.ndarray: The mean point of the specified atoms.
        """
        atoms_indices = adjust_indices(atoms_indices)
        return np.mean(self.coordinates_array[atoms_indices], axis=0)
 
    def get_coordination_transformation_df(self, base_atoms_indices: List[int],origin=None) -> pd.DataFrame:
        """
        Returns a new DataFrame with the coordinates transformed based on the specified base atoms.

        Args:
            base_atoms_indices (List[int]): The indices of the base atoms to use for the transformation.

        Returns:
            pd.DataFrame: A new DataFrame with the transformed coordinates.
        """
        new_coordinates_df = preform_coordination_transformation(self.xyz_df, base_atoms_indices,origin)
        return new_coordinates_df

    ## not working after renumbering for some reason
    
    
    def get_npa_df_single(self, atoms: List[int], sub_atoms: Optional[List[int]] = None, type='nbo') -> pd.DataFrame:
        """
        Calculates the NPA charges for a single group of atoms.
        
        Parameters
        ----------
        atoms : List[int]
            The indices of the base atoms used to perform the coordination transformation.
        sub_atoms : Optional[List[int]]
            If provided, only these atom indices (rows) will be taken from the transformed coordinates.
            Otherwise, all atoms (rows) are used.
        type : str, optional
            Type of charge dictionary to use (default 'nbo').
        
        Returns
        -------
        pd.DataFrame
            DataFrame with the calculated NPA charges, renamed by the base atoms.
        """
        if sub_atoms ==[]:
            sub_atoms = None
        # Get the transformed coordinates DataFrame for the given atoms.
        df_trans = self.get_coordination_transformation_df(atoms,sub_atoms)
        # If sub_atoms is provided, select only those rows.
        if sub_atoms is not None:
            df_subset = df_trans.loc[sub_atoms, ['x', 'y', 'z']].astype(float)
        else:
            df_subset = df_trans[['x', 'y', 'z']].astype(float)
        coordinates_array = np.array(df_subset)
        
        # Get the charges from the charge dictionary.
        charges = np.array(self.charge_dict[type])

        # Calculate the NPA charges (assuming calc_npa_charges is defined elsewhere)
        npa_df = calc_npa_charges(coordinates_array, charges)
        npa_df = npa_df.rename(index={0: f'NPA_{atoms[0]}-{atoms[1]}-{atoms[2]}'})
        return npa_df

    def get_npa_df(self, base_atoms_indices: List[int], sub_atoms: Optional[List[int]] = None, type='nbo') -> pd.DataFrame:
        """
        Returns a DataFrame with the NPA charges calculated based on the specified base atoms 
        and optionally only the sub atoms.
        
        Parameters
        ----------
        base_atoms_indices : List[int] or List[List[int]]
            The indices of the base atoms (or groups of atoms) to use for the NPA calculation.
        sub_atoms : Optional[List[int]] or Optional[List[List[int]]]
            The indices of the sub atoms to use for the NPA calculation. If provided, only these rows 
            will be taken from the transformed coordinates.
        type : str, optional
            Type of charge dictionary to use (default 'nbo').
        
        Returns
        -------
        pd.DataFrame
            A DataFrame with the calculated NPA charges.
        """
        # If the second element is a list, then base_atoms_indices is a list of groups.
        if isinstance(base_atoms_indices[1], list):
            # If sub_atoms is provided as a list of lists, then zip them.
            if sub_atoms and isinstance(sub_atoms[0], list):
                npa_list = [
                    self.get_npa_df_single(atoms, sub_atoms=sub_group, type=type)
                    for atoms, sub_group in zip(base_atoms_indices, sub_atoms)
                ]
            else:
                # Otherwise, pass the same sub_atoms list to each group.
                npa_list = [
                    self.get_npa_df_single(atoms, sub_atoms=sub_atoms, type=type)
                    for atoms in base_atoms_indices
                ]
            npa_df = pd.concat(npa_list, axis=0)
        else:
            npa_df = self.get_npa_df_single(base_atoms_indices, sub_atoms=sub_atoms, type=type)
        return npa_df
    
   
    # def get_dipole_gaussian_df_single(self, atoms, visualize_bool: bool = False) -> pd.DataFrame:
    #     """
    #     atoms can be any of:
    #     - [o, y, p]                         # single-origin, single-plane atom
    #     - [o, y, [p1, p2, ...]]             # single-origin, multi-plane set
    #     - [[o1, o2, ...], y, p]             # multi-origin set, single-plane atom
    #     - [[o1, o2, ...], y, [p1, p2, ...]] # multi-origin set, multi-plane set

    #     All indices may be 1-based or 0-based (auto-detected per array).
    #     Origin point used for viz is the centroid of the origin set.
    #     """

    #     # --- helpers (local) ---
    #     def _as_int_list(x):
    #         if isinstance(x, (list, tuple, set)):
    #             return [int(i) for i in x]
    #         return [int(x)]

    #     def _parse_origin_y_plane(group):
    #         # group: [origin, y, plane]
    #         if len(group) < 3:
    #             raise ValueError("atoms must be [origin, y, plane] where origin/plane can be scalars or lists")
    #         origin_set = _as_int_list(group[0])
    #         y_idx      = int(group[1])
    #         plane_set  = _as_int_list(group[2])
    #         return origin_set, y_idx, plane_set

    #     def _to0(idx_list):
    #         arr = np.asarray(idx_list, dtype=int)
    #         # If any are 0, assume already 0-based
    #         return arr if (arr == 0).any() else (arr - 1)

    #     # ---- parse input ----
    #     origin_set, y_idx, plane_set = _parse_origin_y_plane(atoms)

    #     # ---- compute dipole ----
    #     # Normalize the atoms spec we pass to the calculator: [origin_set, y, plane_set]
    #     atoms_norm = [origin_set, y_idx, plane_set]
    #     dipole_df = calc_dipole_gaussian(
    #         self.coordinates_array,
    #         np.array(self.gauss_dipole_df),
    #         atoms_norm
    #     )

    #     # ---- build a descriptive label ----
    #     def _fmt_set(xs):
    #         return f"{{{','.join(map(str, xs))}}}" if len(xs) > 1 else f"{xs[0]}"
    #     label = f"dipole_{_fmt_set(origin_set)}-{y_idx}-{_fmt_set(plane_set)}"
    #     try:
    #         dipole_df = dipole_df.rename(index={0: label})
    #     except Exception:
    #         dipole_df = dipole_df.rename(index={0: "dipole_gaussian_base"})

    #     # ---- optional visualization ----
    #     if visualize_bool:
    #         try:
    #             coords = np.asarray(self.coordinates_array, dtype=float)
    #             origin0 = _to0(origin_set)
    #             origin_point = coords[origin0].mean(axis=0)

    #             # If your transformer expects a flat list [origin..., y, plane...], provide that
    #             flat_spec = [*origin_set, y_idx, *plane_set]
    #             try:
    #                 xyz_df = self.get_coordination_transformation_df(flat_spec)
    #             except Exception:
    #                 xyz_df = getattr(self, "xyz_df", None)

    #             visualize.show_single_molecule(
    #                 molecule_name=self.molecule_name,
    #                 xyz_df=xyz_df,
    #                 dipole_df=dipole_df,
    #                 origin=origin_point
    #             )
    #         except Exception as e:
    #             print(f"[visualize] Skipping visualization due to: {e}")

    #     return dipole_df

    def get_dipole_gaussian_df_single(self, atoms, visualize_bool: bool = False) -> pd.DataFrame:
        """
        atoms can be:
        - [o, y, plane]
        - [o1, o2, ..., y, plane]
        - [[o1, o2, ...], y, plane]
        """
        # --- helpers (local) ---
        def _parse_base(group):
            # Returns origin_set (list[int]), y (int), plane (int)
            if isinstance(group[0], (list, tuple)):           # [[o...], y, plane]
                origin_set = list(group[0]); y = int(group[1]); plane = int(group[2])
            elif len(group) >= 4:                              # [o..., y, plane]
                origin_set = list(group[:-2]); y = int(group[-2]); plane = int(group[-1])
            else:                                              # [o, y, plane]
                origin_set = [int(group[0])]; y = int(group[1]); plane = int(group[2])
            return origin_set, y, plane

        def _to0(idx_list):
            arr = np.asarray(idx_list, dtype=int)
            return arr if (arr == 0).any() else (arr - 1)

        origin_set, y_idx, plane_idx = _parse_base(atoms)

        # Compute dipole (uses new calc_dipole_gaussian WITHOUT origin)
        dipole_df = calc_dipole_gaussian(self.coordinates_array,
                                        np.array(self.gauss_dipole_df),
                                        atoms)

        # Nice index label (support origin-set)
        if len(origin_set) == 1:
            label = f"dipole_{origin_set[0]}-{y_idx}-{plane_idx}"
        else:
            label = f"dipole_{{{','.join(map(str, origin_set))}}}-{y_idx}-{plane_idx}"

        try:
            dipole_df = dipole_df.rename(index={0: label})
        except Exception:
            dipole_df = dipole_df.rename(index={0: "dipole_gaussian_base"})

        # Optional visualization
        if visualize_bool:
            try:
                # 1) Get transformed coords (LOCAL frame)
                xyz_df = self.get_coordination_transformation_df([*origin_set, y_idx, plane_idx])

                # 2) Compute origin IN THE SAME FRAME AS xyz_df
                if isinstance(xyz_df, (pd.DataFrame,)) and {"x", "y", "z"}.issubset(xyz_df.columns):
                    origin_idx = _to0(origin_set)
                    origin_point = xyz_df.loc[origin_idx, ["x", "y", "z"]].to_numpy().mean(axis=0)
                else:
                    # If the helper returns something else but is already centered → use (0,0,0)
                    origin_point = np.zeros(3)

                print(f"[DEBUG] Visualization origin point (LOCAL): {origin_point}")

                # 3) (Nice-to-have) attach coords for auto-scaling downstream
                try:
                    dipole_df.mol_coords = xyz_df[["x", "y", "z"]].to_numpy()
                except Exception:
                    pass

                # 4) Visualize: dipole_df is already in the same LOCAL frame → no basis needed
                visualize.show_single_molecule(
                    molecule_name=self.molecule_name,
                    xyz_df=xyz_df,
                    dipole_df=dipole_df,
                    origin=origin_point
                )

            except Exception as e:
                # Fallback disabled to avoid frame mismatches (local vs world)
                print(f"[visualize] Skipping visualization (frame mismatch likely): {e}")
                
        return dipole_df


    def get_dipole_gaussian_df(self, base_atoms_indices, visualize_bool: bool = False) -> pd.DataFrame:
        """
        Accepts either a single group:       [[1,2,3], 5, 6]  or  [1,2,3, 5, 6]  or  [1, 5, 6]
        Or multiple groups:                  [ [[1,2,3],5,6], [ [4,7],10,11 ], [2,8,9,12] ]
        """
        def _is_group_like(x):
            return isinstance(x, (list, tuple)) and len(x) >= 3

        # Multiple groups if the top-level contains only group-like items
        if isinstance(base_atoms_indices, (list, tuple)) and all(_is_group_like(g) for g in base_atoms_indices):
            dipole_list = [
                self.get_dipole_gaussian_df_single(g, visualize_bool=visualize_bool)
                for g in base_atoms_indices
            ]
            return pd.concat(dipole_list, axis=0)

        # Otherwise treat as a single group (including [[1,2,3],5,6])
        return self.get_dipole_gaussian_df_single(base_atoms_indices, visualize_bool=visualize_bool)



    def get_bond_angle(self, atom_indices: List[int]) -> pd.DataFrame:
        """
        Returns a DataFrame with the bond angles calculated based on the specified atom indices.

        Args:
            atom_indices (List[int]): The indices of the atoms to use for the bond angle calculation.

        Returns:
            pd.DataFrame: A DataFrame with the bond angles.
        """
        return get_angle_df(self.coordinates_array, atom_indices)
    
    def get_bond_length_single(self, atom_pair):
        bond_length = calc_single_bond_length(self.coordinates_array, atom_pair)
        bond_length_df = pd.DataFrame([bond_length], index=[f'bond_length_{atom_pair[0]}-{atom_pair[1]}'])
        return bond_length_df

    def get_bond_length(self, atom_pairs):
        """
            Returns a DataFrame with the bond lengths calculated based on the specified atom pairs.

            Args:
                atom_pairs (Union[List[Tuple[int, int]], Tuple[int, int]]): The pairs of atoms to use for the bond length calculation.

            Returns:
                pd.DataFrame: A DataFrame with the bond lengths.
            """
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            bond_length_list = [self.get_bond_length_single(pair) for pair in atom_pairs]
            bond_df = pd.concat(bond_length_list, axis=0)
        else:
            # If atom_pairs is a single pair, just process that pair
            bond_df = self.get_bond_length_single(atom_pairs)
        return bond_df


    def get_vibration_max_frequency(self, vibration_atom_num: int) -> float:
        """
        Returns the maximum frequency magnitude for the vibration data of the specified atom number.

        Args:
            vibration_atom_num (int): The atom number for which to retrieve the vibration data.

        Returns:
            float: The maximum frequency magnitude for the vibration data of the specified atom number.
        """
        vibration_array: Union[List[float], None] = self.vibration_dict.get('vibration_atom_{}'.format(str(vibration_atom_num)))
        return calc_max_frequency_magnitude(vibration_array, self.info_df.T)

    def get_stretch_vibration_single(self, atom_pair: List[int],threshold=1600,upper_threshold=3000)-> pd.DataFrame:
        
       
        if check_pair_in_bonds(atom_pair, self.bonds_df) == True:
            try:
                extended_vib_df = calc_vibration_dot_product_from_pairs(
                    self.coordinates_array, self.vibration_dict, atom_pair, self.info_df,threshold=threshold,upper_threshold=upper_threshold,vibration_mode_dict=self.vibration_mode_dict
                )
            except TypeError:
                print(f'Strech Vibration Error: no vibration array for the molecule {self.molecule_name} for {atom_pair} - check atom numbering in molecule')
                return None
            # print(extended_vib_df)
            vibration_df, idx = calc_max_frequency_gen_vibration(extended_vib_df) ## remove for later versions
            return vibration_df.rename(index={idx: f'Stretch_{atom_pair[0]}_{atom_pair[1]}'})
            # return extended_vib_df
        else:
            print(f'Strech Vibration Error: the following bonds do not exist-check atom numbering in molecule: \n {self.molecule_name} for {atom_pair} \n')
            
            df=pd.DataFrame([[np.nan,np.nan]],columns=[['Frequency','Amplitude']])
            df.rename(index={0: f'Stretch_{atom_pair[0]}_{atom_pair[1]}'},inplace=True)
            
            return df
    
    ### check if all frequencies are the same for different vibrations and remove, check what atoms are moving the most
    ## return what atoms are moving the most 
    ### split all the vibrations to symmetric and asymmetric
    ## if you found symmetric look for the asymmetric

    def get_stretch_vibration(self, atom_pairs: List[int],threshold=1600,upper_threshold=3500)-> pd.DataFrame:
        """
        Parameters
        ----------
        atom_pairs : molecule_1.get_stretch_vibration([[1,6],[3,4]])
            atom pairs must have a corresponding vibration file and appear in bonds_df.

        Returns
        -------
        dataframe

        [                   0
         IR            0.7310
         Frequency  1655.5756,
                             0
         IR            0.35906
         Frequency  1689.59450]

        """
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            vibration_list = [self.get_stretch_vibration_single(pair,threshold,upper_threshold) for pair in atom_pairs]
            # Filter out None results
            vibration_list = [vib for vib in vibration_list if vib is not None]
            vibration_df = pd.concat(vibration_list, axis=0)
            
        else:
            # If atom_pairs is a single pair, just process that pair
            vibration_df=self.get_stretch_vibration_single(atom_pairs,threshold,upper_threshold)
        
        # go over df to check if there are equal frequencies and remove them, saving the
        return vibration_df
        
    
    def get_ring_vibrations(
        self,
        ring_atom_indices: List[List[int]],
        *,
        # NEW: thresholds (all optional, applied post-filter)
        freq_min = None,
        freq_max = None,
        inten_min = None,
        inten_max = None,
        return_nan_on_empty: bool = True,   # if True -> return a NaN row when filters empty
        verbose: bool = True,
        **filter_kwargs,                    # still passed to get_filtered_ring_df
    ) -> pd.DataFrame:
        """
        Parameters
        ----------
        ring_atom_indices : list
            Either a single list like [1] or [1,4], or a list of such lists, e.g. [[1],[2,5],...].
            Example: molecule_1.get_ring_vibrations([6])

        freq_min/freq_max : float or None
            Optional min/max threshold for vibrational frequency (cm^-1).

        inten_min/inten_max : float or None
            Optional min/max threshold for vibrational intensity (a.u. or km/mol).

        Returns
        -------
        pd.DataFrame or None
            If multiple inputs: stacked rows per atom-set. Index is labeled (e.g., 'cross_[1, 4]').
            Columns: cross, cross_angle, para, para_angle (from calc_min_max_ring_vibration).
            If filters (including thresholds) produce no data and `return_nan_on_empty=True`,
            returns a single NaN row; otherwise returns None.
        """
        import numpy as np
        import pandas as pd

        def _pick_col(df: pd.DataFrame, candidates: list[str]) -> str :
            """Return the first column name present in df from candidates (or None)."""
            for c in candidates:
                if c in df.columns:
                    return c
            return None

        def _apply_thresholds(df: pd.DataFrame) -> pd.DataFrame:
            """Apply freq/intensity thresholds if provided; return filtered df."""
            if df is None or getattr(df, "empty", False):
                return df

            # Try common column names
            freq_col = _pick_col(df, ["frequency", "freq", "wavenumber", "wn"])
            inten_col = _pick_col(df, ["intensity", "inten", "IR_intensity", "IR", "I"])

            out = df
            if freq_col is not None:
                if freq_min is not None:
                    out = out[out[freq_col] >= float(freq_min)]
                if freq_max is not None:
                    out = out[out[freq_col] <= float(freq_max)]

            if inten_col is not None:
                if inten_min is not None:
                    out = out[out[inten_col] >= float(inten_min)]
                if inten_max is not None:
                    out = out[out[inten_col] <= float(inten_max)]

            return out

        def _nan_row(atom_set):
            nan_row = pd.DataFrame(
                [[np.nan, np.nan, np.nan, np.nan]],
                columns=["cross", "cross_angle", "para", "para_angle"]
            )
            nan_row.rename(index={0: f"empty_{atom_set}"}, inplace=True)
            return nan_row

        def _process_one(atom_set):
            # Resolve benzene partner indices
            try:
                z, x, c, v, b, n = get_benzene_ring_indices(self.bonds_df, atom_set)
                ring_atom_indices_group = [[z, x], [c, v], [b, n]]
            except Exception as e:
                if verbose:
                    print(f"[ERROR] get_benzene_ring_indices failed for {atom_set}: {e}")
                return None

            # Filter vibrations (pre-existing filters via filter_kwargs)
            try:
                filtered_df = get_filtered_ring_df(
                    self.info_df, self.coordinates_array, self.vibration_dict, ring_atom_indices_group,
                    **filter_kwargs
                )
            except FileNotFoundError:
                if verbose:
                    print(f"[ERROR] No vibration - check atom numbering in molecule {getattr(self, 'molecule_name', 'Unknown')}")
                return None
            except Exception as e:
                if verbose:
                    print(f"[ERROR] get_filtered_ring_df failed for {atom_set}: {e}")
                return None

            # NEW: apply min/max thresholds (freq/intensity)
            filtered_df = _apply_thresholds(filtered_df)

            # Handle empty after thresholds
            if filtered_df is None or getattr(filtered_df, "empty", False):
                if verbose:
                    print("No Vibration data found within the specified thresholds. Adjust your thresholds.")
                return _nan_row(atom_set) if return_nan_on_empty else None

            # Safe compute min/max summarization
            try:
                df = calc_min_max_ring_vibration(filtered_df, ring_atom_indices)
            except ValueError as e:
                if verbose:
                    print(f"[WARN] Summary failed for {atom_set}: {e}")
                return _nan_row(atom_set) if return_nan_on_empty else None
            except Exception as e:
                if verbose:
                    print(f"[ERROR] calc_min_max_ring_vibration failed for {atom_set}: {e}")
                return None

            return df

        try:
            # Normalize input shape to a list of lists
            if not isinstance(ring_atom_indices, list):
                ring_atom_indices = [ring_atom_indices]
            if ring_atom_indices and not isinstance(ring_atom_indices[0], list):
                ring_atom_indices = [ring_atom_indices]

            out = []
            for atoms in ring_atom_indices:
                res = _process_one(atoms)
                if res is not None:
                    out.append(res)

            if out:
                return pd.concat(out, axis=0)
            else:
                if verbose:
                    print(f"[INFO] No ring vibration data produced for {getattr(self, 'molecule_name', 'Unknown')}.")
                return None

        except Exception as e:
            print(f"[ERROR] Unexpected error in get_ring_vibrations for molecule {getattr(self, 'molecule_name', 'Unknown')}: {e}")
            return None

 
    
    # def get_bend_vibration_single(self, atom_pair: List[int], threshold: float = 1300)-> pd.DataFrame:
    #     # Create the adjacency dictionary for the pair of atoms
    #     atom_pair = [a - 1 for a in atom_pair]
    #     adjacency_dict = create_adjacency_dict_for_pair(self.bonds_df, atom_pair)
  
    #     center_atom_exists = find_center_atom(atom_pair[0], atom_pair[1], adjacency_dict)
    #     if not center_atom_exists:
    #         raise ValueError(f'Bend Vibration - Atoms do not share a center in molecule {self.molecule_name} - for atoms {atom_pair} check atom numbering in molecule')
    #     else:
    #         # Create the extended DataFrame for the vibration modes
    #         extended_df = extended_df_for_stretch(self.vibration_dict, self.info_df, atom_pair, threshold)
    #         # Calculate the cross product and magnitude for each row in the extended DataFrame
    #         cross_list = []
    #         cross_mag_list = []
    #         for i in range(extended_df.shape[0]):
    #             cross_list.append(np.cross(extended_df[[0, 1, 2]].iloc[i], extended_df[[3, 4, 5]].iloc[i]))
    #             cross_mag_list.append(np.linalg.norm(cross_list[-1]))
    #         # Add the cross magnitude column to the extended DataFrame
    #         extended_df['Cross_mag'] = cross_mag_list
    #         extended_df.reset_index(drop=True, inplace=True)
    #         # Find the row with the maximum cross magnitude and extract the bending frequency and cross magnitude

    #         index_max = extended_df['Cross_mag'].idxmax()
    #         max_frequency_vibration = (pd.DataFrame(extended_df.iloc[index_max]).T)[['Frequency', 'Cross_mag']]
    #         max_frequency_vibration = max_frequency_vibration.rename(index={index_max: f'Bending_{atom_pair[0]}-{atom_pair[1]}'})
    #         return max_frequency_vibration

    def get_bend_vibration_single(self, atom_pair: List[int], threshold: float = 1300) -> pd.DataFrame:
        print(f"[get_bend_vibration_single] molecule={self.molecule_name}, atom_pair={atom_pair}, threshold={threshold}")

        # bonds_df is 1-based
        adjacency_dict = create_adjacency_dict_for_pair(self.bonds_df, atom_pair)
        center_atom_exists = find_center_atom(atom_pair[0], atom_pair[1], adjacency_dict)

        if not center_atom_exists:
            raise ValueError(
                f'Bend Vibration - Atoms do not share a center in molecule {self.molecule_name} '
                f'- for atoms {atom_pair} check atom numbering in molecule'
            )

        # extended_df_for_stretch expects 0-based
        # atom_pair_0 = [a - 1 for a in atom_pair]
        print(f"[get_bend_vibration_single] converted to 0-based: {atom_pair}")
        
        extended_df = extended_df_for_stretch(self.vibration_dict, self.info_df, atom_pair, threshold)
        print(f"[get_bend_vibration_single] extended_df shape={extended_df.shape}")
        print(f"[get_bend_vibration_single] extended_df head:\n{extended_df.head()}")

        if extended_df.empty:
            raise ValueError(f"[get_bend_vibration_single] No modes found above threshold {threshold} for {self.molecule_name}")

        cross_list, cross_mag_list = [], []
        for i in range(extended_df.shape[0]):
            cross = np.cross(extended_df[[0, 1, 2]].iloc[i], extended_df[[3, 4, 5]].iloc[i])
            cross_list.append(cross)
            cross_mag_list.append(np.linalg.norm(cross))

        extended_df['Cross_mag'] = cross_mag_list
        extended_df.reset_index(drop=True, inplace=True)

        index_max = extended_df['Cross_mag'].idxmax()
        print(f"[get_bend_vibration_single] index_max={index_max}, max Cross_mag={cross_mag_list[index_max]:.4f}")
        print(f"[get_bend_vibration_single] winning row:\n{extended_df.iloc[index_max]}")

        max_frequency_vibration = pd.DataFrame(extended_df.iloc[index_max]).T[['Frequency', 'Cross_mag']]
        max_frequency_vibration = max_frequency_vibration.rename(
            index={index_max: f'Bending_{atom_pair[0]}-{atom_pair[1]}'}
        )
        print(f"[get_bend_vibration_single] result:\n{max_frequency_vibration}")
        return max_frequency_vibration

    def get_bend_vibration(self, atom_pairs: List[str], threshold: float = 1300) -> pd.DataFrame:
        """"

        Finds the bending frequency for a pair of atoms in a molecule.

        Args from self:
            atom_pair (list): A list of two atom symbols representing the pair of atoms to find the bending frequency for.
            vibration_dict (dict): A dictionary that maps each vibration mode to a list of its frequencies.
            bonds_df (pd.DataFrame): A DataFrame that contains the bond information for the molecule.
            info_df (pd.DataFrame): A DataFrame that contains the information for each vibration mode in the molecule.
            threshold (float): The frequency threshold for selecting vibration modes.

        Returns:
            pd.DataFrame: If the atoms do not share a center, returns a string indicating this. Otherwise, returns a DataFrame that contains the bending frequency and cross magnitude for the pair of atoms.

        """
        
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            vibration_list = [self.get_bend_vibration_single(pair, threshold) for pair in atom_pairs]
            vibration_df = pd.concat(vibration_list, axis=0)
            return vibration_df
        else:
            # If atom_pairs is a single pair, just process that pair
            return self.get_bend_vibration_single(atom_pairs, threshold)



class Molecules():
    
    def __init__(self,molecules_dir_name, renumber=False, threshold=1.82):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[]
        self.failed_molecules=[]
        self.success_molecules=[]
        for feather_file in os.listdir(): 
            if feather_file.endswith('.feather'):
                try:
                    self.molecules.append(Molecule(feather_file, threshold=threshold))
                    self.success_molecules.append(feather_file)
                except Exception as e:
                    self.failed_molecules.append(feather_file)
                    print(f'Error: {feather_file} could not be processed : {e}')
                   
        print(f'Molecules Loaded: {self.success_molecules}',f'Failed Molecules: {self.failed_molecules}')

        self.molecule_names=[molecule.molecule_name for molecule in self.molecules]
        self.old_molecules=self.molecules
        self.old_molecule_names=self.molecule_names
        os.chdir('../')

    def export_all_xyz(self):
        os.makedirs('xyz_files', exist_ok=True)
        os.chdir('xyz_files')
        for mol in self.molecules:
            xyz_df=mol.xyz_df
            mol_name=mol.molecule_name
            data_to_xyz(xyz_df, f'{mol_name}.xyz')
        os.chdir('../')
    
    def renumber_all_molecules(self,renumbering_list, target_idx=None):
        import copy
        # if len(renumbering_list) != len(self.molecules):
        #     raise ValueError("Renumbering list length must match the number of molecules.")
        # if target_idx is -1 translate to 0
        if target_idx is not None and target_idx == -1:
            target_idx = 0
        
        for i, (renumbering_dict, molecule) in enumerate(zip(renumbering_list, self.molecules)):
            if target_idx == i:
                if i + 1 < len(renumbering_list):
                    merged = copy.deepcopy(renumbering_list[i])
                    merged.update(copy.deepcopy(renumbering_dict))  # current has priority
                    renumbering_list[i + 1] = merged
                print(f'Skipping renumbering for molecule {molecule.molecule_name} at index {i}. dict : {renumbering_dict}')
                continue
            # try:
            molecule.renumber_atoms(renumbering_dict, bonds_are_0_based=True)
            print(f'Molecule {molecule.molecule_name} renumbered successfully.')
            # except Exception as e:
            #     print(f'Error renumbering molecule {molecule.molecule_name}: {e}')
            #     log_exception("renumber_all_molecules")
        print('All molecules renumbered.')
    
    def filter_molecules(self, indices):
        self.molecules = [self.molecules[i] for i in indices]
        self.molecule_names = [self.molecule_names[i] for i in indices]

    def get_sterimol_dict(self,atom_indices, radii='CPK',sub_structure=True, drop_atoms=None,visualize_bool=None,mode='all'):
        """
        Returns a dictionary with the Sterimol parameters calculated based on the specified base atoms.

        Args:
            base_atoms (Union[None, List[int, int], List[List[int, int]]]): The indices of the base atoms to use for the Sterimol calculation.
            radii (str, optional): The radii to use for the Sterimol calculation. Defaults to 'CPK'.
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the Sterimol parameters.
        
        input example: molecules.get_sterimol_dict([[1,6],[3,4]])
        output example: 
        Results for LS1716_optimized:
            B1    B5     L  loc_B1  loc_B5
        1-2  1.91  8.12  8.89    5.35    6.51
        3-4  1.70  5.25  9.44    1.70   -1.65


        Results for LS1717_optimized:
            B1    B5     L  loc_B1  loc_B5
        1-2  1.9  7.98  9.02    3.64    6.68
        3-4  1.7  6.45  9.45    1.70   -1.97
        """
        sterimol_dict={}
        
        for molecule in self.molecules:
            try:
                
                sterimol_dict[molecule.molecule_name]=molecule.get_sterimol(atom_indices, radii, sub_structure=sub_structure, drop_atoms=drop_atoms, visualize_bool=visualize_bool, mode=mode)
            except Exception as e:
                print(f'Error: {molecule.molecule_name} sterimol could not be processed: {e}')
                log_exception("get_sterimol_dict")
                pass
        return dict_to_horizontal_df(sterimol_dict)
    
    def get_npa_dict(self,atom_indices,sub_atoms=None):
        """
        Returns a dictionary with the Natural Population Analysis (NPA) charges calculated for the specified base atoms and sub atoms.

        Args:
            base_atoms (List[int]): The indices of the base atoms to use for the NPA calculation.
            sub_atoms (Union[List[int], None], optional): The indices of the sub atoms to use for the NPA calculation. Defaults to None.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the NPA charges.
        
        input example: molecules.get_npa_dict([1, 2, 3], [5, 6, 7])
        output example: 
        Results for LS1716_optimized:
                    dip_x     dip_y     dip_z  total_dipole
            NPA_1-2-3  0.092108  0.181346 -0.300763      0.363082
            NPA_5-6-7  0.191986  0.297839  0.079456      0.363153

        Results for LS1717_optimized:
                    dip_x     dip_y     dip_z  total_dipole
            NPA_1-2-3  0.126370  0.271384  0.354595      0.464065
            NPA_5-6-7  0.225257  0.399616 -0.069960      0.464035
        """
        
        npa_dict={}
        for molecule in self.molecules:
            try:
                npa_dict[molecule.molecule_name]=molecule.get_npa_df(atom_indices,sub_atoms)
            except Exception as e:
                print(f'Error: {molecule.molecule_name} npa could not be processed: {e}')
                log_exception("get_npa_dict")
                pass

        return dict_to_horizontal_df(npa_dict)

    def get_ring_vibration_dict(self, ring_atom_indices, freq_min=1550, freq_max=1700):
        """
        Parameters
        ----------
        ring_atom_indices : list
            Can be a single list of indices [a, b] or a list of lists like [[a, b], [c, d]].
            Example:
                molecule_1.get_ring_vibrations([[8,11],[9,12]])

        For example – for a ring of atoms 1–6 where atom 4 is connected to the main group
        and 1 is para to it (ortho would be 3 & 5, meta would be 2 & 6),
        you’d enter [1,4].

        Returns
        -------
        pd.DataFrame
            Each molecule’s ring vibration data combined horizontally into one table.
        """
        import traceback
        import pandas as pd

        ring_dict = {}

        # normalize to list-of-lists
        if not isinstance(ring_atom_indices, list):
            ring_atom_indices = [ring_atom_indices]
        if ring_atom_indices and not isinstance(ring_atom_indices[0], list):
            ring_atom_indices = [ring_atom_indices]

        for molecule in self.molecules:
            try:
                combined_df_list = []

                # iterate over multiple ring atom sets if given
                for atom_set in ring_atom_indices:
                    res = molecule.get_ring_vibrations(atom_set, freq_min=freq_min, freq_max=freq_max)
                    if res is not None and not getattr(res, "empty", False):
                        # name each subset by atom indices
                        label = "-".join(map(str, atom_set))
                        res.columns = [f"{col}_{label}" for col in res.columns]
                        combined_df_list.append(res)

                # concatenate all ring results for this molecule
                if combined_df_list:
                    ring_dict[molecule.molecule_name] = pd.concat(combined_df_list, axis=1)
                else:
                    ring_dict[molecule.molecule_name] = pd.DataFrame()

            except Exception as e:
                print(f"Error: {molecule.molecule_name} ring vibration could not be processed: {e}")
                traceback.print_exc()
                continue

        return dict_to_horizontal_df(ring_dict)
    

    def get_dipole_dict(self, atom_indices, visualize_bool: bool = False):
        """
        Returns a wide (horizontal) DataFrame aggregating dipole results for each molecule.

        Parameters
        ----------
        atom_indices : sequence
            Base-atom spec (R-style):
            - [o, y, plane]
            - [o1, o2, ..., y, plane]
            - [[o1, o2, ...], y, plane]
            Or a list of such groups: [ [...], [...], ... ]
        visualize_bool : bool
            Whether to visualize per-molecule dipole.

        Returns
        -------
        pd.DataFrame
            A horizontally concatenated DataFrame; each molecule's dipole DataFrame is a block.
        """
        dipole_dict = {}
        for molecule in self.molecules:
            try:
                df = molecule.get_dipole_gaussian_df(atom_indices, visualize_bool=visualize_bool)
                # If you want to match older prints that used 'total_dipole', uncomment next line:
                # df = df.rename(columns={'total': 'total_dipole'})
                dipole_dict[molecule.molecule_name] = df
            except Exception as e:
                print(f"Error: {molecule.molecule_name} dipole could not be processed: {e}")
                log_exception("get_dipole_dict")

        return dict_to_horizontal_df(dipole_dict)

    
    def get_bond_angle_dict(self,atom_indices):
        """
        Returns a dictionary with the bond angles calculated for the specified atom indices.

        Args:
            atom_indices (List[List[int, int]]) or : The indices of the atoms to use for the bond angle calculation.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the bond angle data.
        
        input example: molecules.get_bond_angle_dict([1, 2, 3])
        output example: 
        Results for LS1716_optimized:
                        Angle
            Angle_1-2-3  109.5

        Results for LS1717_optimized:
                        Angle
            Angle_1-2-3  108.7
        """

        bond_angle_dict={}
        for molecule in self.molecules:
            try:
                bond_angle_dict[molecule.molecule_name]=molecule.get_bond_angle(atom_indices)
            except Exception as e:
                print(f'Error: {molecule.molecule_name} Angle could not be processed: {e}')
                log_exception("get_bond_angle_dict")
                pass
        return dict_to_horizontal_df(bond_angle_dict)
    
    def get_bond_length_dict(self,atom_pairs):
        """
        Returns a dictionary with the bond lengths calculated for the specified atom pairs.

        Args:
            atom_pairs  (List[List[int, int]]) : The pairs of atoms to use for the bond length calculation.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the bond length data.
        
        input example: molecules.get_bond_length_dict([[1, 2], [3, 4]])
        output example: 
        Results for LS1716_optimized:
                        0
        bond_length_1-2  1.529754
        bond_length_3-4  1.466212


        Results for LS1717_optimized:
                                0
        bond_length_1-2  1.511003
        bond_length_3-4  1.466089
        """

        bond_length_dict={}
        for molecule in self.molecules:
            try:
                bond_length_dict[molecule.molecule_name]=molecule.get_bond_length(atom_pairs)
            except Exception as e:
                print(f'Error: {molecule.molecule_name} Bond Length could not be processed: {e}')
                
                pass
        return dict_to_horizontal_df(bond_length_dict)
    
    def get_stretch_vibration_dict(self,atom_pairs,threshold=1400,upper_threshold=3500):
        """
        Returns a dictionary with the stretch vibrations calculated for the specified atom pairs.

        Args:
            atom_pairs (List[Tuple[int, int]]): The pairs of atoms to use for the stretch vibration calculation.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the stretch vibration data.
        
        input example: molecules.get_stretch_vibration_dict([[1, 2], [3, 4]])
        output example: 
        Results for LS1716_optimized:
                        Frequency  Amplitude
            Stretch_1_2  3174.3565   0.330304
            Stretch_3_4  3242.4530   0.170556

        Results for LS1717_optimized:
                        Frequency  Amplitude
            Stretch_1_2  3242.4465   0.252313
            Stretch_3_4  3175.4029   0.443073
    """
        stretch_vibration_dict={}
        print(f'Calculating stretch vibration for atoms {atom_pairs} with threshold {threshold} \n Remember : ALWAYS LOOK AT THE RESULTING VIBRATION')
        for molecule in self.molecules:
            try:
                stretch_vibration_dict[molecule.molecule_name]=molecule.get_stretch_vibration(atom_pairs,threshold,upper_threshold)
            except Exception as e:
                print(f'Error: could not calculate strech vibration for {molecule.molecule_name}: {e} ')
                log_exception("get_stretch_vibration_dict")
                pass

        return dict_to_horizontal_df(stretch_vibration_dict)
    
    def get_charge_df_dict(self,atom_indices):
        """
        Returns a dictionary with the Natural Bond Orbital (NBO) charges for the specified atoms.

        Args:
            atoms_indices (List[int]): The indices of the atoms to include in the NBO charge calculation.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the NBO charges.
        
        input example: molecules.get_charge_df_dict([3, 5, 7, 9])
        output example: 
        Results for LS1716_optimized:
                    atom_3   atom_5   atom_7   atom_9
            charge -0.12768 -0.39006  0.14877 -0.00656

        Results for LS1717_optimized:
                    atom_3   atom_5   atom_7   atom_9
            charge -0.12255 -0.38581  0.14691 -0.00681
        """

        nbo_dict={}
        for molecule in self.molecules:
            try:
                nbo_dict[molecule.molecule_name]=molecule.get_charge_df(atom_indices,type='all')
            except Exception as e:
                print(f'Error: could not calculate nbo value for {molecule.molecule_name}: {e} ')
                log_exception("get_charge_df_dict")
                pass
        return charge_dict_to_horizontal_df(nbo_dict)
    
    def get_charge_diff_df_dict(self,atom_indices,type='all'):
        """
        Returns a dictionary with the differences in Natural Bond Orbital (NBO) charges for the specified pairs of atoms.

        Args:
            diff_indices (List[List[int]]): The indices of the atom pairs to calculate the NBO charge differences for.
        
        Returns:
            Dict[str, pd.DataFrame]: A dictionary where each key is a molecule name and each value is a DataFrame with the NBO charge differences.
        
        input example: molecules.get_charge_diff_df_dict([[1, 2], [3, 4]])
        output example: 
        Results for LS1716_optimized:
                        diff_1-2  diff_3-4
            charge -0.12768 -0.39006  0.14877 -0.00656


            Results for LS1717_optimized:
                        diff_1-2  diff_3-4
            charge -0.12255 -0.38581  0.14691 -0.00681
        """
        charge_diff_dict={}
        for molecule in self.molecules:
            try:
                charge_diff_dict[molecule.molecule_name]=molecule.get_charge_diff_df(atom_indices,type)
            except Exception as e:
                print(f'Error: could not calculate nbo difference for {molecule.molecule_name}: {e} ')
                log_exception("get_charge_diff_df_dict")
                pass
        return charge_dict_to_horizontal_df(charge_diff_dict)
    
    def get_bend_vibration_dict(self, atom_pairs, threshold=1300):
        bending_dict = {}

        for molecule in self.molecules:
            molecule_results = []
            
            for pair in atom_pairs:
                try:
                    # Call a version that handles a SINGLE pair
                    res = molecule.get_bend_vibration_single(pair, threshold)
                    molecule_results.append(res)
                except Exception as e:
                    # Log that this specific pair failed, but keep going
                    print(f"Skipping pair {pair} for {molecule.molecule_name}: {e}")
                    continue
            
            # Only add to the dictionary if at least one pair succeeded
            if molecule_results:
                bending_dict[molecule.molecule_name] = pd.concat(molecule_results)

        
        return dict_to_horizontal_df(bending_dict)
    
    def visualize_molecules(self,indices=None):
        if indices is not None:
            for idx in indices:
                self.molecules[idx].visualize_molecule()
        else:
            for molecule in self.molecules:
                molecule.visualize_molecule()
        
    def visualize_smallest_molecule(self):
        """
        Visualizes the smallest molecule based on the number of atoms.

        Args:
            None
        
        Returns:
            None
        
        input example: molecules.visualize_smallest_molecule()
        output example: 
        # This will open a visualization window or generate a visualization file for the smallest molecule.

        """
        idx=0
        smallest= len(self.molecules[0].xyz_df)
        for id, molecule in enumerate(self.molecules[1:]):
            if len(molecule.xyz_df)<smallest:
                smallest=len(molecule.xyz_df)
                idx=id
        html=self.molecules[idx].visualize_molecule()
        return html
    
    def visualize_smallest_molecule_morfeus(self, indices=None):
        idx=0
        mol=self.molecules[0]
        smallest= len(mol.xyz_df)
        for id, molecule in enumerate(self.molecules[1:]):
            if len(molecule.xyz_df)<smallest:
                smallest=len(molecule.xyz_df)
                idx=id
        mol=self.molecules[idx]
        mol.get_sterimol(indices,visualize_bool=True)
        
    def get_molecules_features_set(self, entry_widgets, parameters=None, answers_list=None, save_as=False, csv_file_name='features_output'):
        """
        Gathers user input from entry_widgets, applies parameters,
        extracts features, and optionally saves results to a file.
        Includes NaN and duplicate-value diagnostics.
        """
        import pandas as pd
        import datetime
        import numpy as np
        import re
        import os

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if save_as:
            csv_file_name = f"{csv_file_name}_{timestamp}"

        # --- 1. GUI answers ---
        answers = {}
        if parameters is None:
            parameters = {'Radii': 'CPK', 'Isotropic': True}

        for param_name, entry in entry_widgets.items():
            key = param_name.split()[0].lower().replace('-', '_')
            try:
                answers[key] = entry.get()
            except AttributeError:
                answers[key] = entry

        if answers_list is not None:
            answers_list = answers_list  # placeholder

        radii = parameters.get('Radii', 'CPK')
        iso = parameters.get('Isotropic', True)

        for k, v in answers.items():
            if v != '':
                if isinstance(v, str):
                    answers[k] = convert_to_list_or_nested_list(v)
                elif isinstance(v, int):
                    answers[k] = [v]
                else:
                    answers[k] = v
            else:
                answers[k] = []

        res_df = pd.DataFrame()

        def safe_concat(res_df, new_df):
            if res_df.empty:
                return new_df
            else:
                return pd.concat([res_df, new_df], axis=1)

        # --- 2. Feature extraction steps ---
        feature_steps = [
            ('ring', lambda a: self.get_ring_vibration_dict(a)),
            ('stretching', lambda a: self.get_stretch_vibration_dict(a, answers.get('stretch', [None])[0], answers.get('upper_stretch', [None])[0])),
            ('bending', lambda a: self.get_bend_vibration_dict(a, answers.get('bend', [None])[0])),
            ('npa', lambda a: self.get_npa_dict(a, sub_atoms=answers.get('sub_atoms', []))),
            ('dipole', lambda a: self.get_dipole_dict(a)),
            ('charges', lambda a: self.get_charge_df_dict(a)),
            ('charge_diff', lambda a: self.get_charge_diff_df_dict(a)),
            ('sterimol', lambda a: self.get_sterimol_dict(a, radii=radii, drop_atoms=answers.get('drop_atoms', []))),
            ('bond_angle', lambda a: self.get_bond_angle_dict(a)),
            ('bond_length', lambda a: self.get_bond_length_dict(a)),
        ]

        # --- 3. Apply all extraction steps ---
        for key, handler in feature_steps:
            if answers.get(key):
                try:
                    new_df = handler(answers[key])
                    res_df = safe_concat(res_df, new_df)
                except Exception as e:
                    print(f"Error processing {key} for {getattr(self.molecules[0], 'molecule_name', 'unknown')}: {e}")
                    log_exception(f"get_molecules_comp_set_app – {key}")
                    continue

        # --- 4. Polarizability & energy block ---
        if iso:
            try:
                rows = []
                for molecule in self.molecules:
                    info = pd.DataFrame(index=[molecule.molecule_name])
                    try:
                        polar = molecule.polarizability_df.copy().iloc[[0]]
                        polar.index = [molecule.molecule_name]
                        info = polar
                    except Exception:
                        pass
                    try:
                        info["energy"] = molecule.energy_value.values
                    except Exception:
                        pass
                    rows.append(info)

                if rows:
                    polar_energy_concat = pd.concat(rows, axis=0)
                    polar_energy_concat = polar_energy_concat.dropna(axis=1, how='all')
                    polar_energy_concat = polar_energy_concat.reset_index().set_index('index')
                    res_df = safe_concat(res_df, polar_energy_concat)
            except Exception as e:
                print(f"Error processing polarizability/Energy: {e}")
                log_exception("get_molecules_comp_set_app – polarizability/energy")

        # --- 5. Correlation analysis (optional visualization) ---
        # interactive_corr_heatmap_with_highlights(res_df)
        correlation_table = show_highly_correlated_pairs(res_df, corr_thresh=0.8)

        # --- 6. Post-processing diagnostics ---
        print("\n===== Feature Integrity Diagnostics =====")

        # (a) Check columns with NaN values
        nan_cols = res_df.columns[res_df.isna().any()].tolist()
        if nan_cols:
            print(f"[!] Columns with missing (NaN) values ({len(nan_cols)}):")
            for c in nan_cols:
                missing_pct = res_df[c].isna().mean() * 100
                print(f"    - {c} ({missing_pct:.1f}% missing)")
        else:
            print("[✓] No columns with NaN values.")

        # (b) Check for identical or nearly identical columns
        numeric_df = res_df.select_dtypes(include=[np.number])
        duplicate_groups = []
        corr = numeric_df.corr().abs()

        corr_arr = corr.to_numpy(copy=True)
        corr_arr[np.tril_indices_from(corr_arr)] = np.nan
        corr = pd.DataFrame(corr_arr, index=corr.index, columns=corr.columns)

        near_dupes = corr.stack()[corr.stack() > 0.999].index.tolist()
                
        if near_dupes:
            print(f"[!] Columns with nearly identical values (>0.999 correlation):")
            for (a, b) in near_dupes:
                print(f"    - {a}  ≈  {b}")
        else:
            print("[✓] No nearly identical numeric columns detected.")

        print("========================================\n")

        # --- 7. Save results if requested ---
        res_df = res_df.sort_index(ascending=False)
        if save_as:
            try:
                res_df = res_df.sort_index(
                    key=lambda idx: idx.map(
                        lambda x: int(re.search(r"\d+", x).group()) if re.search(r"\d+", x) else x
                    )
                )
            except Exception:
                pass

            res_df.to_csv(f"{csv_file_name}.csv", index=True)
            correlation_table.to_csv(f"{csv_file_name}_correlation_table.csv", index=True)
            print(f"Features saved to {csv_file_name}.csv and correlation table to {csv_file_name}_correlation_table.csv in {os.getcwd()}")

        return res_df


    def get_renumbering_dict(self):
        
        self.export_all_xyz()
        os.chdir('xyz_files')
        self.renumbering_list, target_idx = batch_renumbering(os.getcwd())
      
        # self.renumbering_list = [
        # {
        #     mol_name: {
        #         **inner, 
        #         **{v: k for k, v in inner.items() if v != k}
        #     }
        #     for mol_name, inner in d.items()
        # }
        # for d in renumbering_list
    # ]
    
        self.renumber_all_molecules(self.renumbering_list, target_idx)
        os.chdir('../')

        
    def extract_all_dfs(self):
        os.chdir(self.molecules_path)
        for molecule in self.molecules:
            try:
                os.mkdir(molecule.molecule_name)
            except FileExistsError:
                pass
            os.chdir(molecule.molecule_name)
            molecule.write_xyz_file()
            molecule.write_csv_files()
            os.chdir('../')
        
    
    def extract_all_xyz(self):
        os.chdir(self.molecules_path)
        try:
            os.mkdir('xyz_files')
        except FileExistsError:
            pass
        os.chdir('xyz_files')
        for molecule in self.molecules:
            data_to_xyz(molecule.xyz_df,(molecule.molecule_name+'.xyz'))
        os.chdir('../')



    # def mols_to_rdkit_mols(self):
    #     mols=[mol for mol in self.molecules]
    #     xyz_list=[xyz.xyz_df for xyz in mols]
    #     molss=xyz_list_to_mols(xyz_list, mol_names=self..molecule_names)
    #     img = visualize_mols(molss)
    #     # display(img) 
    #     fps = mols_to_fingerprint_df(molss, fp_type="morgan", nBits=1024)

    #     fp_df, pca_df, var_ratio = mols_to_fingerprint_df(
    #     mols,
    #     fp_type="morgan",
    #     nBits=2048,
    #     radius=2,
    #     do_pca=True,
    #     n_components=20,
    #     scale_before_pca=True
    # )

