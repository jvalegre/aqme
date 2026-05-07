import pandas as pd
import numpy as np
from enum import Enum
from typing import List, Dict, Any, Union


class Names(Enum):
    DIPOLE_COLUMNS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
    STANDARD_ORIENTATION_COLUMNS = ['atom', 'x', 'y', 'z']
    DF_LIST = ['standard_orientation_df', 'dipole_df', 'pol_df', 'info_df','energy_df']


def df_list_to_dict(df_list: List[pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """Maps list of DataFrames to names specified in Names.DF_LIST."""
    return {name: df for name, df in zip(Names.DF_LIST.value, df_list)}


def split_to_vib_dict(vectors: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Splits vibration vector columns into separate NumPy arrays per atom."""
    num_columns = vectors.shape[1]
    vib_dict = {}
    for i in range(num_columns // 3):
        key = f'vibration_atom_{i + 1}'
        vib_dict[key] = vectors.iloc[:, i*3:(i+1)*3].dropna().to_numpy(dtype=float)
    return vib_dict



def feather_file_handler(feather_file: str) -> List[Any]:
    """
    Reads a feather file with molecular data, processes, and returns:
    [df_dict, vib_dict, charge_dict, ev_df]
    Uses columns by name instead of index.
    """
    # Define the expected columns (for error-proof selection)
    expected_cols = [
        'atom', 'x', 'y', 'z', 'dip_x', 'dip_y', 'dip_z', 'total_dipole', 'aniso', 'iso',
        'nbo_hirshfeld', 'hirshfeld_charge', 'cm5_charge', 'Frequency', 'IR'
    ]
    try:
        data = pd.read_feather(feather_file)
        # Optionally, fix column names in case they are not matching exactly
        data.columns = [c.strip() for c in data.columns]

        # XYZ
        xyz_cols = ['atom', 'x', 'y', 'z']
        xyz = data[xyz_cols].dropna().reset_index(drop=True)
        xyz[['x', 'y', 'z']] = xyz[['x', 'y', 'z']].astype(float)

        # Energy value (optional—get from a suitable column or from metadata)
        try:
            ev = pd.DataFrame({'energy': [float(data.get('energy', [None])[0])]})
        except Exception as e:
            # print(f"Error extracting energy value: {e}")
            ev = pd.DataFrame({'energy': [None]})

        # Charges
        nbo_charge_df = data[['nbo_charge']].dropna().rename(columns={'nbo_charge': 'charge'}).astype(float).reset_index(drop=True)
        hirshfeld_charge_df = data[['hirshfeld_charge']].dropna().rename(columns={'hirshfeld_charge': 'charge'}).astype(float).reset_index(drop=True)
        cm5_charge_df = data[['cm5_charge']].dropna().rename(columns={'cm5_charge': 'charge'}).astype(float).reset_index(drop=True)

        # Align charges with xyz
        for cdf in [nbo_charge_df, hirshfeld_charge_df, cm5_charge_df]:
            if len(cdf) != len(xyz):
                cdf = cdf[cdf.index.isin(xyz.index)].reset_index(drop=True)

        charge_dict = {'nbo': nbo_charge_df, 'hirshfeld': hirshfeld_charge_df, 'cm5': cm5_charge_df}

        # Dipole
        dipole_cols = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
        dipole_df = data[dipole_cols].dropna().astype(float).reset_index(drop=True)

        # Polarizability
        pol_cols = ['aniso', 'iso']
        pol_df = data[pol_cols].dropna().astype(float).reset_index(drop=True)

        # Info
        info_cols = ['Frequency', 'IR']
        info_df = data[info_cols].dropna().astype(float).reset_index(drop=True)

        # Vibration vectors: everything else after IR, or define your own selection
        vib_start = data.columns.get_loc('IR') + 1
        vectors = data.iloc[:, vib_start:].dropna(axis=1, how='all')
        vib_dict = split_to_vib_dict(vectors)

        # Package DataFrames into a dict (use your own helper)
        df_list = [xyz, dipole_df, pol_df, info_df, ev]
        df_list = [df for df in df_list if not df.empty and df.dropna(how='all').shape[0] > 0]
        df_dict = df_list_to_dict(df_list)

        return [df_dict, vib_dict, charge_dict]
    except Exception as e:
        data = pd.read_feather(feather_file)
        data.columns = range(len(data.columns))

        # Extract core DataFrames
        xyz = data.iloc[:, 0:4].dropna()
        try:
            xyz.columns = Names.STANDARD_ORIENTATION_COLUMNS.value
            xyz = xyz.reset_index(drop=True)
            xyz[['x', 'y', 'z']] = xyz[['x', 'y', 'z']].astype(float)
        except Exception:
            # fallback: sometimes header is off by a few rows
            xyz = xyz.iloc[2:].reset_index(drop=True)
            xyz.columns = Names.STANDARD_ORIENTATION_COLUMNS.value
            xyz[['x', 'y', 'z']] = xyz[['x', 'y', 'z']].astype(float)
        xyz = xyz.dropna()
        # Extract energy value as single-row DataFrame
        try:
            ev = pd.DataFrame({'energy': [float(data.iloc[1, 1])]})
        except Exception as e:
            # print(f"Error extracting energy value: {e}")
            ev = pd.DataFrame({'energy': [None]})

        # Extract charge DataFrames
        nbo_charge_df = data.iloc[:, 10:11].dropna().rename(columns={10: 'charge'}).astype(float).reset_index(drop=True)
        hirshfeld_charge_df = data.iloc[:, 11:12].dropna().rename(columns={11: 'charge'}).astype(float).reset_index(drop=True)
        cm5_charge_df = data.iloc[:, 12:13].dropna().rename(columns={12: 'charge'}).astype(float).reset_index(drop=True)

        # Ensure charge DataFrames align with xyz atoms
        for cdf in [nbo_charge_df, hirshfeld_charge_df, cm5_charge_df]:
            if len(cdf) != len(xyz):
                cdf = cdf[cdf.index.isin(xyz.index)].reset_index(drop=True)

        charge_dict = {'nbo': nbo_charge_df, 'hirshfeld': hirshfeld_charge_df, 'cm5': cm5_charge_df}

        # Dipole, polarizability, info, vibration vectors
        dipole_df = data.iloc[:, 4:8].dropna().rename(
            columns={4: 'dip_x', 5: 'dip_y', 6: 'dip_z', 7: 'total_dipole'}
        ).astype(float).reset_index(drop=True)

        pol_df = data.iloc[:, 8:10].dropna().rename(
            columns={8: 'aniso', 9: 'iso'}
        ).astype(float).reset_index(drop=True)

        info_df = data.iloc[:, 13:15].dropna().rename(
            columns={13: 'Frequency', 14: 'IR'}
        ).astype(float).reset_index(drop=True)

        # Vibration vectors as dict of numpy arrays
        vectors = data.iloc[:, 15:].dropna()
        vib_dict = split_to_vib_dict(vectors)

        # Main DataFrames as dict
        df_list = [xyz, dipole_df, pol_df, info_df]
        df_list = [df for df in df_list if not df.empty and df.dropna(how='all').shape[0] > 0]
        df_dict = df_list_to_dict(df_list)

        return [df_dict, vib_dict, charge_dict, ev]


def save_to_feather(
    df_dict: Dict[str, pd.DataFrame],
    vib_dict: Dict[str, Any],
    charge_dict: Dict[str, pd.DataFrame],
    ev_df: pd.DataFrame,
    out_file: str
) -> None:
    """
    Reconstructs a single DataFrame from processed molecular pieces
    and saves it as a feather file in the same format as feather_file_handler expects.

    Parameters
    ----------
    df_dict : dict
        Dictionary of DataFrames (xyz, dipole, pol, info, etc.).
    vib_dict : dict
        Dictionary of vibration vectors (mode -> numpy array / DataFrame).
    charge_dict : dict
        Dictionary of charge DataFrames: {'nbo', 'hirshfeld', 'cm5'}.
    ev_df : pd.DataFrame
        Single-row DataFrame with energy value.
    out_file : str
        Path to save the feather file.
    """

    # --- Core DataFrames reconstruction ---
    xyz = df_dict.get("xyz", pd.DataFrame(columns=["atom", "x", "y", "z"]))
    dipole_df = df_dict.get("dipole", pd.DataFrame(columns=["dip_x", "dip_y", "dip_z", "total_dipole"]))
    pol_df = df_dict.get("polarizability", pd.DataFrame(columns=["aniso", "iso"]))
    info_df = df_dict.get("info", pd.DataFrame(columns=["Frequency", "IR"]))

    # Charges
    nbo_df = charge_dict.get("nbo", pd.DataFrame(columns=["charge"]))
    hirsh_df = charge_dict.get("hirshfeld", pd.DataFrame(columns=["charge"]))
    cm5_df = charge_dict.get("cm5", pd.DataFrame(columns=["charge"]))

    # Rename charge columns to match original names
    nbo_df = nbo_df.rename(columns={"charge": "nbo_charge"})
    hirsh_df = hirsh_df.rename(columns={"charge": "hirshfeld_charge"})
    cm5_df = cm5_df.rename(columns={"charge": "cm5_charge"})

    # Energy
    ev_val = None
    print(ev_df)
    if ev_df is not None and not ev_df.empty:
        try:
            ev_val = float(ev_df.iloc[0, 0])
            # check if
        except Exception as e:
            # print(f"Error extracting energy value: {e}")
            ev_val = None
    ev_col = pd.Series([ev_val] * len(xyz), name="energy")

    # --- Vibrations reconstruction ---
    vib_frames = []
    for mode, arr in vib_dict.items():
        # arr could be numpy array or DataFrame
        if isinstance(arr, pd.DataFrame):
            frame = arr
        else:
            frame = pd.DataFrame(arr)
        frame.columns = [f"{mode}_{i}" for i in range(frame.shape[1])]
        vib_frames.append(frame)

    vib_df = pd.concat(vib_frames, axis=1) if vib_frames else pd.DataFrame()

    # --- Concatenate everything into one big DataFrame ---
    main_df = pd.concat(
        [xyz, dipole_df, pol_df, nbo_df, hirsh_df, cm5_df, info_df, vib_df],
        axis=1
    )

    # Add energy column (single value broadcasted)
    main_df.insert(0, "energy", ev_col)

    # Save as feather
    main_df.reset_index(drop=True).to_feather(out_file)

    print(f"Feather file saved: {out_file}")
