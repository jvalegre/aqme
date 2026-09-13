#!/usr/bin/env python3
"""
extract_chem_json.py

Standalone command-line extractor for Gaussian and ORCA quantum-chemistry
output files. This is the notebook logic from
code/cclib_json_feather_gaussian_orca_hybrid_cleaned.ipynb, pulled out into a
plain script so it can be run from the command line without opening Jupyter.

It uses cclib as the primary parser, with custom regex fallbacks for values
cclib does not expose the same way (NBO/Hirshfeld/CM5 charges,
polarizability, HOMO-LUMO gap, Gaussian vibration blocks). For each input
file it writes:

  <output>/json_files/<name>.json     - structured JSON payload
  <output>/feather_files/<name>.feather - wide table, same layout as the
                                          original pipeline (requires pyarrow)

Install requirements:
    pip install cclib numpy pandas pyarrow --break-system-packages

Usage
-----
Single file (program is auto-detected from the file):

    python extract_chem_json.py path/to/BA044.out

    python extract_chem_json.py path/to/BA044.log -o converted_outputs

Whole folder, only ORCA outputs:

    python extract_chem_json.py path/to/folder --program orca

Whole folder, only Gaussian logs, including subfolders:

    python extract_chem_json.py path/to/folder --program gaussian --recursive

Whole folder, both Gaussian and ORCA (default):

    python extract_chem_json.py path/to/folder --program both

JSON only (skip Feather / skip needing pyarrow):

    python extract_chem_json.py path/to/folder --no-feather

Old parent/com layout (parent_dir/molecule/com/*.out|*.log):

    python extract_chem_json.py path/to/parent_dir --parent-mode
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Pattern, Tuple, Union

import numpy as np
import pandas as pd

try:
    import cclib
    from cclib.io import ccread
except ImportError:
    print("ERROR: cclib is not installed. Run:", file=sys.stderr)
    print("  pip install cclib numpy pandas pyarrow --break-system-packages", file=sys.stderr)
    sys.exit(1)

import json


# ---------------------------------------------------------------------------
# 1. Constants and shared helpers
# ---------------------------------------------------------------------------

FLOAT = r"[-+]?[0-9]+\.[0-9]+"
FLOAT_POL = r"-?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?"
FLOATS_ONLY = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

DIPOLE_START = "Dipole moment"
DIPOLE_END = "Quadrupole moment"
MAIN_CUTOFF = r"\-{69}"
STANDARD_ORIENTATION_START = "Standard orientation:"
POL_START = "iso"
POL_END = "xx"
CHARGE_START = "Summary of"
CHARGE_END = "====="
FREQUENCY_START = "Harmonic frequencies"
FREQUENCY_END = "Thermochemistry"
HIRSH_CHARGE_START = "Hirshfeld"
HIRSH_CHARGE_END = "Hirshfeld charges with hydrogens "

DIPOLE_COLUMNS = ["dip_x", "dip_y", "dip_z", "total_dipole"]
STANDARD_ORIENTATION_COLUMNS = ["atom", "x", "y", "z"]

ATOMIC_SYMBOLS = [
    None,
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]


def atomno_to_symbol(atomno: Union[int, float, str]) -> str:
    try:
        z = int(float(atomno))
        if 0 < z < len(ATOMIC_SYMBOLS) and ATOMIC_SYMBOLS[z] is not None:
            return ATOMIC_SYMBOLS[z]
        return str(z)
    except Exception:
        return str(atomno)


def read_text_file(path: Union[str, Path]) -> str:
    return Path(path).expanduser().resolve().read_text(encoding="utf-8", errors="replace")


def extract_lines_from_text(text: str, re_expression: Union[str, Pattern[str]]) -> List[str]:
    pattern = re_expression if isinstance(re_expression, re.Pattern) else re.compile(re_expression)
    return pattern.findall(text)


def find_all_matches(text: str, key_phrase: Union[str, Pattern[str]]):
    pattern = key_phrase if isinstance(key_phrase, re.Pattern) else re.compile(key_phrase)
    return list(pattern.finditer(text))


def fortran_float(s: str) -> float:
    return float(str(s).replace("D", "E").replace("d", "e"))


def df_with_column(values: Optional[Iterable[Any]], column: str, n_rows: Optional[int] = None) -> pd.DataFrame:
    if values is None:
        if n_rows is None:
            return pd.DataFrame(columns=[column])
        return pd.DataFrame({column: [np.nan] * n_rows})
    values = list(values)
    if not values and n_rows is not None:
        values = [np.nan] * n_rows
    return pd.DataFrame({column: values})


def dataframe_to_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    if df is None or df.empty:
        return []
    cleaned = df.replace({np.nan: None})
    return cleaned.to_dict(orient="records")


def make_json_safe(obj: Any) -> Any:
    if obj is None:
        return None
    if isinstance(obj, np.ndarray):
        return make_json_safe(obj.tolist())
    if isinstance(obj, np.generic):
        return make_json_safe(obj.item())
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, (str, int, bool)):
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, bytes):
        return obj.decode("utf-8", errors="replace")
    if isinstance(obj, dict):
        return {str(make_json_safe(k)): make_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [make_json_safe(v) for v in obj]
    if isinstance(obj, pd.DataFrame):
        return dataframe_to_records(obj)
    if hasattr(obj, "tolist"):
        try:
            return make_json_safe(obj.tolist())
        except Exception:
            pass
    return str(obj)


# ---------------------------------------------------------------------------
# 2. cclib parsing helpers
# ---------------------------------------------------------------------------

def read_cclib_output(output_file: Union[str, Path]) -> Tuple[Optional[Any], Optional[str]]:
    output_file = Path(output_file).expanduser().resolve()
    try:
        data = ccread(str(output_file))
        if data is None:
            return None, "cclib returned None"
        return data, None
    except Exception as exc:
        return None, str(exc)


def get_cclib_attributes(data: Optional[Any]) -> Dict[str, Any]:
    if data is None:
        return {}
    if hasattr(data, "getattributes"):
        try:
            return data.getattributes(tolists=True)
        except TypeError:
            return data.getattributes()
        except Exception:
            pass

    attrs = {}
    for name in dir(data):
        if name.startswith("_"):
            continue
        try:
            value = getattr(data, name)
        except Exception:
            continue
        if not callable(value):
            attrs[name] = value
    return attrs


def detect_program(output_file: Union[str, Path], text: str, cclib_attrs: Dict[str, Any]) -> str:
    metadata = cclib_attrs.get("metadata", {}) if isinstance(cclib_attrs, dict) else {}
    package = str(metadata.get("package", "")).lower() if isinstance(metadata, dict) else ""

    if "orca" in package or "o   r   c   a" in text.lower() or "orca terminated" in text.lower():
        return "orca"
    if "gaussian" in package or "entering gaussian system" in text.lower() or "standard orientation:" in text.lower():
        return "gaussian"

    suffix = Path(output_file).suffix.lower()
    if suffix in {".orcaout"}:
        return "orca"
    return "unknown"


def atoms_from_cclib(data: Optional[Any]) -> List[Dict[str, Any]]:
    if data is None:
        return []

    atomnos = getattr(data, "atomnos", None)
    atomcoords = getattr(data, "atomcoords", None)
    if atomnos is None or atomcoords is None:
        return []

    atomnos = np.asarray(atomnos)
    coords = np.asarray(atomcoords, dtype=float)
    if coords.ndim == 3:
        coords = coords[-1]
    if coords.ndim != 2 or coords.shape[1] != 3:
        return []

    atoms = []
    for idx, (z, xyz) in enumerate(zip(atomnos, coords), start=1):
        atoms.append({
            "atom_index": idx,
            "atomic_number": int(z),
            "atom": atomno_to_symbol(z),
            "x": float(xyz[0]),
            "y": float(xyz[1]),
            "z": float(xyz[2]),
        })
    return atoms


def standard_orientation_df_from_atoms(atoms: List[Dict[str, Any]], include_current_header_rows: bool = True) -> pd.DataFrame:
    if not atoms:
        return pd.DataFrame(columns=STANDARD_ORIENTATION_COLUMNS)

    df = pd.DataFrame([{
        "atom": atom["atom"],
        "x": atom["x"],
        "y": atom["y"],
        "z": atom["z"],
    } for atom in atoms], columns=STANDARD_ORIENTATION_COLUMNS)

    if not include_current_header_rows:
        return df

    n_atoms = df.shape[0]
    headers = pd.DataFrame(
        [[np.nan] * df.shape[1], [n_atoms] + [np.nan] * (df.shape[1] - 1)],
        columns=df.columns,
    )
    return pd.concat([headers, df], ignore_index=True)


def homo_lumo_gap_from_cclib(data: Optional[Any]) -> Optional[float]:
    if data is None:
        return None
    homos = getattr(data, "homos", None)
    moenergies = getattr(data, "moenergies", None)
    if homos is None or moenergies is None:
        return None

    try:
        homos = np.asarray(homos, dtype=int)
        moenergies = [np.asarray(spin, dtype=float) for spin in moenergies]
        spin_index = 0
        homo_index = int(homos[spin_index])
        if homo_index + 1 >= len(moenergies[spin_index]):
            return None
        return float(moenergies[spin_index][homo_index + 1] - moenergies[spin_index][homo_index])
    except Exception:
        return None


def cclib_scf_energies_ev(data: Optional[Any]) -> List[float]:
    if data is None or not hasattr(data, "scfenergies"):
        return []
    try:
        return [float(x) for x in np.asarray(data.scfenergies, dtype=float).ravel()]
    except Exception:
        return []


def dipole_df_from_cclib(data: Optional[Any]) -> pd.DataFrame:
    if data is None or not hasattr(data, "moments"):
        return pd.DataFrame(columns=DIPOLE_COLUMNS)
    try:
        moments = getattr(data, "moments")
        if moments is None or len(moments) < 2:
            return pd.DataFrame(columns=DIPOLE_COLUMNS)
        dip = np.asarray(moments[1], dtype=float).ravel()
        if dip.size < 3:
            return pd.DataFrame(columns=DIPOLE_COLUMNS)
        values = [float(dip[0]), float(dip[1]), float(dip[2]), float(np.linalg.norm(dip[:3]))]
        return pd.DataFrame([values], columns=DIPOLE_COLUMNS)
    except Exception:
        return pd.DataFrame(columns=DIPOLE_COLUMNS)


def info_df_from_cclib(data: Optional[Any]) -> pd.DataFrame:
    if data is None or not hasattr(data, "vibfreqs"):
        return pd.DataFrame(columns=["Frequency", "IR"])

    try:
        freqs = np.asarray(getattr(data, "vibfreqs"), dtype=float).ravel()
    except Exception:
        return pd.DataFrame(columns=["Frequency", "IR"])

    if hasattr(data, "vibirs"):
        try:
            irs = np.asarray(getattr(data, "vibirs"), dtype=float).ravel()
        except Exception:
            irs = np.full(freqs.shape, np.nan)
    else:
        irs = np.full(freqs.shape, np.nan)

    n = min(len(freqs), len(irs))
    return pd.DataFrame({"Frequency": freqs[:n], "IR": irs[:n]}).astype(float)


def vibs_df_from_cclib(data: Optional[Any]) -> pd.DataFrame:
    if data is None or not hasattr(data, "vibdisps"):
        return pd.DataFrame()

    try:
        vibdisps = np.asarray(getattr(data, "vibdisps"), dtype=float)
    except Exception:
        return pd.DataFrame()

    if vibdisps.ndim != 3 or vibdisps.shape[2] != 3:
        return pd.DataFrame()

    frames = []
    for mode_index in range(vibdisps.shape[0]):
        frames.append(pd.DataFrame(
            vibdisps[mode_index],
            columns=[
                f"mode_{mode_index + 1}_x",
                f"mode_{mode_index + 1}_y",
                f"mode_{mode_index + 1}_z",
            ],
        ))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=1).astype(float)


def charges_from_cclib(data: Optional[Any], n_atoms: Optional[int] = None) -> Dict[str, pd.DataFrame]:
    empty = {
        "nbo": df_with_column(None, "nbo_charge", n_atoms),
        "hirshfeld": df_with_column(None, "hirshfeld_charge", n_atoms),
        "cm5": df_with_column(None, "cm5_charge", n_atoms),
    }

    if data is None or not hasattr(data, "atomcharges"):
        return empty

    atomcharges = getattr(data, "atomcharges")
    if not isinstance(atomcharges, dict):
        return empty

    def pick_charge(keys: List[str]) -> Optional[List[float]]:
        lower_map = {str(k).lower(): v for k, v in atomcharges.items()}
        for key in keys:
            if key.lower() in lower_map:
                try:
                    return [float(x) for x in np.asarray(lower_map[key.lower()], dtype=float).ravel()]
                except Exception:
                    return None
        return None

    return {
        "nbo": df_with_column(pick_charge(["nbo", "natural", "mulliken"]), "nbo_charge", n_atoms),
        "hirshfeld": df_with_column(pick_charge(["hirshfeld"]), "hirshfeld_charge", n_atoms),
        "cm5": df_with_column(pick_charge(["cm5"]), "cm5_charge", n_atoms),
    }


def vibrations_json_from_cclib(data: Optional[Any]) -> List[Dict[str, Any]]:
    if data is None:
        return []

    info_df = info_df_from_cclib(data)
    vibs_df = vibs_df_from_cclib(data)

    if info_df.empty:
        return []

    modes = []
    n_modes = info_df.shape[0]
    for mode_index in range(n_modes):
        mode = {
            "mode_index": mode_index + 1,
            "frequency_cm-1": float(info_df.iloc[mode_index]["Frequency"]) if pd.notna(info_df.iloc[mode_index]["Frequency"]) else None,
            "ir_intensity": float(info_df.iloc[mode_index]["IR"]) if pd.notna(info_df.iloc[mode_index]["IR"]) else None,
            "displacements": [],
        }

        base = mode_index * 3
        if not vibs_df.empty and base + 2 < vibs_df.shape[1]:
            cols = vibs_df.iloc[:, base:base + 3]
            for atom_index, row in enumerate(cols.to_numpy(dtype=float), start=1):
                mode["displacements"].append({
                    "atom_index": atom_index,
                    "dx": float(row[0]),
                    "dy": float(row[1]),
                    "dz": float(row[2]),
                })

        modes.append(mode)

    return modes


# ---------------------------------------------------------------------------
# 3. Gaussian fallback parser helpers
# ---------------------------------------------------------------------------

def process_gaussian_charge_text(log_text: str) -> pd.DataFrame:
    starts = find_all_matches(log_text, CHARGE_START)
    ends = find_all_matches(log_text, CHARGE_END)
    if not starts or not ends:
        raise ValueError("Could not locate charge start/end markers in log text.")

    start_idx = starts[-1].end()
    end_idx = ends[-1].start()
    section = log_text[start_idx:end_idx]

    floats = extract_lines_from_text(section, FLOATS_ONLY)
    if not floats:
        raise ValueError("No floating-point numbers found in charge section.")

    arr = np.array(floats, dtype=float)
    charge_values = arr[1::6]
    return pd.DataFrame({"nbo_charge": charge_values})


def process_gaussian_dipole_text(log_text: str) -> pd.DataFrame:
    starts = find_all_matches(log_text, DIPOLE_START)
    ends = find_all_matches(log_text, DIPOLE_END)
    if not starts or not ends:
        raise ValueError("Could not locate dipole start/end markers in log text.")

    start_idx = starts[-1].end()
    end_idx = ends[-1].start()
    block = log_text[start_idx:end_idx]

    floats = extract_lines_from_text(block, FLOATS_ONLY)
    if len(floats) < 4:
        raise ValueError(f"Expected 4 dipole components, found {len(floats)}.")

    values = list(map(float, floats[:4]))
    return pd.DataFrame([values], columns=DIPOLE_COLUMNS)


def gauss_first_split(log_text: str) -> List[str]:
    after_std = re.split(STANDARD_ORIENTATION_START, log_text)[-1]
    return re.split(MAIN_CUTOFF, after_std)


def process_gaussian_standard_orientation_text(log_text: str) -> pd.DataFrame:
    float_lines = extract_lines_from_text(log_text, FLOATS_ONLY)
    arr = np.array(float_lines, dtype=float)

    if arr.size % 6 != 0:
        raise ValueError("Floating-point count is not a multiple of 6.")

    mat = arr.reshape(-1, 6)
    coords = np.delete(mat, (0, 2), axis=1)
    df = pd.DataFrame(coords, columns=STANDARD_ORIENTATION_COLUMNS)

    df["atom"] = df["atom"].astype(int).map(atomno_to_symbol)
    for col in ["x", "y", "z"]:
        df[col] = df[col].astype(float)

    n_atoms = df.shape[0]
    headers = pd.DataFrame(
        [[np.nan] * df.shape[1], [n_atoms] + [np.nan] * (df.shape[1] - 1)],
        columns=df.columns,
    )
    return pd.concat([headers, df], ignore_index=True)


def search_phrase_in_text_pol(text: str, key_phrase: str) -> Optional[int]:
    pattern = re.compile(rf"{re.escape(str(key_phrase))}\s")
    match = pattern.search(text)
    return match.start() if match else None


def process_gaussian_pol_text(log_text_or_lines: Union[str, List[str]]) -> pd.DataFrame:
    log_text = "\n".join(log_text_or_lines) if isinstance(log_text_or_lines, list) else log_text_or_lines

    pol_start = search_phrase_in_text_pol(log_text, POL_START)
    pol_end = search_phrase_in_text_pol(log_text, POL_END)

    if pol_start is None or pol_end is None or pol_end <= pol_start:
        return pd.DataFrame([100.0, 100.0], index=["iso", "aniso"]).T

    chunk = log_text[pol_start:pol_end]
    pol = extract_lines_from_text(chunk, FLOAT_POL)
    if len(pol) < 4:
        return pd.DataFrame([100.0, 100.0], index=["iso", "aniso"]).T

    iso = fortran_float(pol[0])
    aniso = fortran_float(pol[3])
    return pd.DataFrame({"iso": [iso], "aniso": [aniso]}, dtype=float)


def process_hirshfeld_charges(log_text: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    charges_start = find_all_matches(log_text, HIRSH_CHARGE_START)
    charges_end = find_all_matches(log_text, HIRSH_CHARGE_END)

    if not charges_start or not charges_end:
        raise ValueError("Could not locate Hirshfeld start/end markers.")

    start_idx = charges_start[-2].end()
    end_idx = charges_end[-1].start()
    section = log_text[start_idx:end_idx]

    selected = extract_lines_from_text(section, FLOAT)
    if not selected:
        raise ValueError("No floats found in Hirshfeld section.")

    arr = np.array(selected, dtype=float)
    rem = arr.size % 6
    if rem != 0:
        arr = arr[:-rem]

    table = arr.reshape(-1, 6)
    hirsh = table[:, 0]
    cm5 = table[:, 5]

    return pd.DataFrame(hirsh, columns=["hirshfeld_charge"]), pd.DataFrame(cm5, columns=["cm5_charge"])


def remove_floats_until_first_int(input_list: List[str]) -> List[str]:
    output_list = []
    encountered_integer = False
    for item in input_list:
        try:
            int(item)
            encountered_integer = True
            output_list.append(item)
        except ValueError:
            try:
                float(item)
                if encountered_integer:
                    output_list.append(item)
            except ValueError:
                output_list.append(item)
    return output_list


def process_gaussian_vibs_string(log_text: str) -> List[str]:
    pattern = re.compile(rf"{FREQUENCY_START}([\s\S]*?){FREQUENCY_END}")
    match = pattern.search(log_text)
    if not match:
        raise ValueError("Could not locate Harmonic frequencies...Thermochemistry section.")

    vibration_section = match.group(1).strip()
    frequencies_blocks = re.split(r"Frequencies --", vibration_section)
    return [block.split("-------------------")[0].strip() for block in frequencies_blocks[1:]]


def process_gaussian_info(frequency_blocks: List[str]) -> pd.DataFrame:
    ir, frequency = [], []
    for data in frequency_blocks:
        match = extract_lines_from_text(data, FLOATS_ONLY)
        frequency.append(match[0:3])
        ir.append(match[6:9])

    info_df = pd.DataFrame()
    info_df["Frequency"] = [item for sub in frequency for item in sub]
    info_df["IR"] = [item for sub in ir for item in sub]
    return info_df.astype(float)


def vib_array_list_to_df(array_list: List[np.ndarray]) -> pd.DataFrame:
    frames = []
    for i, array in enumerate(array_list, 1):
        new_array = np.delete(array, [0, 1], axis=1).reshape(-1, 3)
        frames.append(pd.DataFrame(new_array, columns=[f"mode_{i}_x", f"mode_{i}_y", f"mode_{i}_z"]))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=1).astype(float)


def process_gaussian_frequency_string(final_blocks: List[str]) -> pd.DataFrame:
    vibs_list = []

    for i, data in enumerate(final_blocks):
        match = re.findall(FLOATS_ONLY, data)
        match = match[12:]
        match = np.array(match, dtype=str)

        if i == 1:
            cleaned = remove_floats_until_first_int(list(match))
            match = np.array(cleaned, dtype=str)

        rem = len(match) % 11
        if rem != 0:
            match = match[:-rem]

        if len(match) == 0:
            raise ValueError(f"Block {i} has no vibrational data after cleaning.")

        vibs_list.append(match.reshape(-1, 11))

    vibs = np.vstack(vibs_list)
    atom_idx_numeric = np.array([int(float(x)) for x in vibs[:, 0]])
    final_atom = int(atom_idx_numeric.max())

    ordered_vibs = []
    for atom_idx in range(1, final_atom + 1):
        ordered_vibs.append(vibs[atom_idx_numeric == atom_idx])

    return vib_array_list_to_df(ordered_vibs)


def process_gaussian_energy_text(log_text: str) -> pd.DataFrame:
    if "SCF Done" not in log_text:
        raise ValueError("No SCF Done found in log text.")

    cut = re.split("SCF Done", log_text)[1]
    energy = extract_lines_from_text(cut, FLOAT)
    if not energy:
        raise ValueError("No energy float found after SCF Done.")

    return pd.DataFrame([[float(energy[0])]], columns=["energy"]).astype(float)


def extract_homo_lumo_gap_from_text(log_text: str) -> Optional[float]:
    lines = log_text.splitlines()
    last_occ_line_idx = None

    for i, line in enumerate(lines):
        if "Alpha  occ. eigenvalues" in line:
            last_occ_line_idx = i

    if last_occ_line_idx is None:
        return None

    homo_values = re.findall(r"-?\d+\.\d+", lines[last_occ_line_idx])
    if not homo_values:
        return None
    homo = float(homo_values[-1])

    lumo = None
    for j in range(last_occ_line_idx + 1, len(lines)):
        if "Alpha virt. eigenvalues" in lines[j]:
            lumo_values = re.findall(r"-?\d+\.\d+", lines[j])
            if not lumo_values:
                return None
            lumo = float(lumo_values[0])
            break

    if lumo is None:
        return None

    return (lumo - homo) * 27.2114


def gaussian_standard_orientation_from_text(log_text: str) -> pd.DataFrame:
    parts = gauss_first_split(log_text)
    if len(parts) < 3:
        raise ValueError("gauss_first_split did not yield expected parts.")
    return process_gaussian_standard_orientation_text(parts[2])


# ---------------------------------------------------------------------------
# 4. ORCA raw-text helpers
# ---------------------------------------------------------------------------

def extract_orca_text_summary(orca_text: str) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}

    version_match = re.search(r"Program Version\s+([0-9][^\s]*)", orca_text)
    if version_match:
        summary["orca_version"] = version_match.group(1)

    energy_matches = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
        orca_text,
    )
    if energy_matches:
        summary["final_single_point_energy_hartree"] = float(energy_matches[-1])

    scf_matches = re.findall(
        r"Total Energy\s*:\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s*Eh",
        orca_text,
    )
    if scf_matches:
        summary["last_total_energy_hartree"] = float(scf_matches[-1])

    summary["terminated_normally"] = "ORCA TERMINATED NORMALLY" in orca_text
    summary["has_error_termination"] = "ORCA finished by error termination" in orca_text

    return summary


# ---------------------------------------------------------------------------
# 5. Feather-compatible table builder
# ---------------------------------------------------------------------------

def prefer_nonempty(primary: pd.DataFrame, fallback: pd.DataFrame) -> pd.DataFrame:
    if primary is not None and not primary.empty:
        return primary
    if fallback is not None:
        return fallback
    return pd.DataFrame()


def build_feather_compatible_df(
    output_file: Union[str, Path],
    data: Optional[Any],
    text: str,
    program: str,
    warnings: List[str],
) -> pd.DataFrame:
    atoms = atoms_from_cclib(data)
    n_atoms = len(atoms) if atoms else None

    standard_orientation_df = standard_orientation_df_from_atoms(atoms, include_current_header_rows=True)
    dipole_df = dipole_df_from_cclib(data)
    pol_df = pd.DataFrame(columns=["iso", "aniso"])
    charge_dfs = charges_from_cclib(data, n_atoms=n_atoms)
    nbo_charge_df = charge_dfs["nbo"]
    hirsh_charge_df = charge_dfs["hirshfeld"]
    cm5_charge_df = charge_dfs["cm5"]
    info_df = info_df_from_cclib(data)
    vibs_df = vibs_df_from_cclib(data)

    if program == "gaussian":
        try:
            custom_std = gaussian_standard_orientation_from_text(text)
            standard_orientation_df = prefer_nonempty(custom_std, standard_orientation_df)
            n_atoms = custom_std.dropna(subset=["x", "y", "z"]).shape[0]
        except Exception as exc:
            warnings.append(f"Gaussian custom standard orientation failed: {exc}")

        try:
            dipole_df = prefer_nonempty(process_gaussian_dipole_text(text).astype(float), dipole_df)
        except Exception as exc:
            warnings.append(f"Gaussian custom dipole failed: {exc}")

        try:
            pol_df = process_gaussian_pol_text(text).astype(float)
        except Exception as exc:
            warnings.append(f"Gaussian custom polarizability failed: {exc}")

        try:
            nbo_charge_df = process_gaussian_charge_text(text).astype(float)
        except Exception as exc:
            warnings.append(f"Gaussian custom NBO charge failed: {exc}")

        try:
            hirsh_charge_df, cm5_charge_df = process_hirshfeld_charges(text)
            hirsh_charge_df = hirsh_charge_df.astype(float)
            cm5_charge_df = cm5_charge_df.astype(float)
        except Exception as exc:
            warnings.append(f"Gaussian custom Hirshfeld/CM5 failed: {exc}")
            if n_atoms is not None:
                hirsh_charge_df = df_with_column(None, "hirshfeld_charge", n_atoms)
                cm5_charge_df = df_with_column(None, "cm5_charge", n_atoms)

        try:
            frequency_blocks = process_gaussian_vibs_string(text)
            custom_info_df = process_gaussian_info(frequency_blocks).astype(float)
            custom_vibs_df = process_gaussian_frequency_string(frequency_blocks)
            info_df = prefer_nonempty(info_df, custom_info_df)
            vibs_df = prefer_nonempty(vibs_df, custom_vibs_df)
        except Exception as exc:
            warnings.append(f"Gaussian custom vibrations failed: {exc}")

    gap_ev = None
    if program == "gaussian":
        gap_ev = extract_homo_lumo_gap_from_text(text)
    if gap_ev is None:
        gap_ev = homo_lumo_gap_from_cclib(data)

    ev_df = pd.DataFrame({"energy": [gap_ev]}) if gap_ev is not None else pd.DataFrame(columns=["energy"])

    pieces = [
        standard_orientation_df,
        dipole_df,
        pol_df,
        ev_df,
        nbo_charge_df,
        hirsh_charge_df,
        cm5_charge_df,
        info_df,
        vibs_df,
    ]

    pieces = [piece for piece in pieces if piece is not None]
    if not pieces:
        return pd.DataFrame()

    return pd.concat(pieces, axis=1)


def save_current_feather(df: pd.DataFrame, feather_file: Union[str, Path]) -> Path:
    feather_file = Path(feather_file).expanduser().resolve()
    feather_file.parent.mkdir(parents=True, exist_ok=True)

    df = df.copy()
    if "atom" in df.columns:
        df["atom"] = df["atom"].astype(str)
    df.reset_index(drop=True, inplace=True)
    df.to_feather(feather_file)
    return feather_file


def save_xyz(atoms: List[Dict[str, Any]], xyz_file: Union[str, Path], comment: str = "") -> Path:
    xyz_file = Path(xyz_file).expanduser().resolve()
    xyz_file.parent.mkdir(parents=True, exist_ok=True)

    lines = [str(len(atoms)), comment]
    for atom in atoms:
        atom_label = str(atom.get('atom', 'X'))
        lines.append(
            f"{atom_label:>2} "
            f"{float(atom.get('x', 0.0)):.8f} "
            f"{float(atom.get('y', 0.0)):.8f} "
            f"{float(atom.get('z', 0.0)):.8f}"
        )

    xyz_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return xyz_file


# ---------------------------------------------------------------------------
# 6. Structured JSON builder
# ---------------------------------------------------------------------------

def extract_gaussian_custom_json(text: str, warnings: List[str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    try:
        dipole_df = process_gaussian_dipole_text(text)
        row = dipole_df.iloc[0]
        out["dipole"] = {
            "x": float(row["dip_x"]),
            "y": float(row["dip_y"]),
            "z": float(row["dip_z"]),
            "total": float(row["total_dipole"]),
            "source": "gaussian_custom_parser",
        }
    except Exception as exc:
        warnings.append(f"JSON Gaussian dipole fallback failed: {exc}")

    try:
        pol_df = process_gaussian_pol_text(text)
        row = pol_df.iloc[0]
        out["polarizability"] = {
            "iso": float(row["iso"]),
            "aniso": float(row["aniso"]),
            "source": "gaussian_custom_parser",
        }
    except Exception as exc:
        warnings.append(f"JSON Gaussian polarizability fallback failed: {exc}")

    try:
        energy_df = process_gaussian_energy_text(text)
        out["scf_energy_hartree"] = float(energy_df.iloc[0]["energy"])
    except Exception as exc:
        warnings.append(f"JSON Gaussian SCF energy fallback failed: {exc}")

    try:
        gap_ev = extract_homo_lumo_gap_from_text(text)
        if gap_ev is not None:
            out["homo_lumo_gap_ev"] = float(gap_ev)
    except Exception as exc:
        warnings.append(f"JSON Gaussian HOMO-LUMO fallback failed: {exc}")

    try:
        nbo_df = process_gaussian_charge_text(text)
        out["nbo_charges"] = [float(x) for x in nbo_df["nbo_charge"].tolist()]
    except Exception as exc:
        warnings.append(f"JSON Gaussian NBO charge fallback failed: {exc}")

    try:
        hirsh_df, cm5_df = process_hirshfeld_charges(text)
        out["hirshfeld_charges"] = [float(x) for x in hirsh_df["hirshfeld_charge"].tolist()]
        out["cm5_charges"] = [float(x) for x in cm5_df["cm5_charge"].tolist()]
    except Exception as exc:
        warnings.append(f"JSON Gaussian Hirshfeld/CM5 fallback failed: {exc}")

    return out


def extract_feather_df_sections(feather_df: pd.DataFrame) -> Dict[str, Any]:
    if feather_df is None or feather_df.empty:
        return {}

    sections: Dict[str, Any] = {}

    try:
        xyz = feather_df.iloc[:, 0:4].copy()
        xyz.columns = ["atom", "x", "y", "z"]
        xyz = xyz.dropna(subset=["x", "y", "z"], how="any")
        sections["current_feather_atoms"] = dataframe_to_records(xyz)
    except Exception:
        pass

    try:
        dipole = feather_df.iloc[:, 4:8].dropna(how="all").copy()
        if dipole.shape[1] == 4:
            dipole.columns = DIPOLE_COLUMNS
            sections["current_feather_dipole"] = dataframe_to_records(dipole)
    except Exception:
        pass

    try:
        pol = feather_df.iloc[:, 8:10].dropna(how="all").copy()
        if pol.shape[1] == 2:
            pol.columns = ["iso", "aniso"]
            sections["current_feather_polarizability"] = dataframe_to_records(pol)
    except Exception:
        pass

    try:
        gap = feather_df.iloc[:, 10:11].dropna(how="all").copy()
        if gap.shape[1] == 1:
            gap.columns = ["homo_lumo_gap_ev"]
            sections["current_feather_homo_lumo_gap"] = dataframe_to_records(gap)
    except Exception:
        pass

    return sections


def build_structured_json(
    output_file: Union[str, Path],
    data: Optional[Any],
    text: str,
    program: str,
    feather_df: pd.DataFrame,
    warnings: List[str],
    include_full_cclib: bool = True,
) -> Dict[str, Any]:
    output_file = Path(output_file).expanduser().resolve()
    cclib_attrs = get_cclib_attributes(data)

    atoms = atoms_from_cclib(data)
    if not atoms and feather_df is not None and not feather_df.empty:
        try:
            xyz = feather_df.iloc[:, 0:4].copy()
            xyz.columns = ["atom", "x", "y", "z"]
            xyz = xyz.dropna(subset=["x", "y", "z"], how="any")
            atoms = []
            for idx, row in enumerate(xyz.to_dict(orient="records"), start=1):
                atoms.append({
                    "atom_index": idx,
                    "atom": row["atom"],
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "z": float(row["z"]),
                })
        except Exception as exc:
            warnings.append(f"Could not build atoms from feather table: {exc}")

    charge_dfs = charges_from_cclib(data, n_atoms=len(atoms) if atoms else None)
    nbo = charge_dfs["nbo"]["nbo_charge"].dropna().astype(float).tolist() if "nbo_charge" in charge_dfs["nbo"] else []
    hirsh = charge_dfs["hirshfeld"]["hirshfeld_charge"].dropna().astype(float).tolist() if "hirshfeld_charge" in charge_dfs["hirshfeld"] else []
    cm5 = charge_dfs["cm5"]["cm5_charge"].dropna().astype(float).tolist() if "cm5_charge" in charge_dfs["cm5"] else []

    dipole_records = dataframe_to_records(dipole_df_from_cclib(data))
    dipole = dipole_records[0] if dipole_records else None

    pol = None
    gaussian_custom = {}
    orca_summary = {}
    if program == "gaussian":
        gaussian_custom = extract_gaussian_custom_json(text, warnings)
        if "dipole" in gaussian_custom:
            dipole = gaussian_custom["dipole"]
        if "polarizability" in gaussian_custom:
            pol = gaussian_custom["polarizability"]
        if gaussian_custom.get("nbo_charges"):
            nbo = gaussian_custom["nbo_charges"]
        if gaussian_custom.get("hirshfeld_charges"):
            hirsh = gaussian_custom["hirshfeld_charges"]
        if gaussian_custom.get("cm5_charges"):
            cm5 = gaussian_custom["cm5_charges"]

    if program == "orca":
        orca_summary = extract_orca_text_summary(text)

    gap_ev = gaussian_custom.get("homo_lumo_gap_ev")
    if gap_ev is None:
        gap_ev = homo_lumo_gap_from_cclib(data)

    scf_energies_ev = cclib_scf_energies_ev(data)
    final_scf_energy_ev = scf_energies_ev[-1] if scf_energies_ev else None

    energies = {
        "scf_energies_ev_cclib": scf_energies_ev,
        "final_scf_energy_ev_cclib": final_scf_energy_ev,
        "scf_energy_hartree_gaussian_custom": gaussian_custom.get("scf_energy_hartree"),
        "final_single_point_energy_hartree_orca": orca_summary.get("final_single_point_energy_hartree"),
        "last_total_energy_hartree_orca": orca_summary.get("last_total_energy_hartree"),
    }

    vibrations = vibrations_json_from_cclib(data)
    if not vibrations:
        try:
            frequency_blocks = process_gaussian_vibs_string(text)
            info_df = process_gaussian_info(frequency_blocks)
            vibs_df = process_gaussian_frequency_string(frequency_blocks)
            n_modes = min(info_df.shape[0], vibs_df.shape[1] // 3)
            vibrations = []
            for mode_index in range(n_modes):
                base = mode_index * 3
                mode = {
                    "mode_index": mode_index + 1,
                    "frequency_cm-1": float(info_df.iloc[mode_index]["Frequency"]),
                    "ir_intensity": float(info_df.iloc[mode_index]["IR"]),
                    "displacements": [],
                    "source": "gaussian_custom_parser",
                }
                for atom_index, row in enumerate(vibs_df.iloc[:, base:base + 3].to_numpy(dtype=float), start=1):
                    mode["displacements"].append({
                        "atom_index": atom_index,
                        "dx": float(row[0]),
                        "dy": float(row[1]),
                        "dz": float(row[2]),
                    })
                vibrations.append(mode)
        except Exception as exc:
            if program == "gaussian":
                warnings.append(f"JSON Gaussian vibration fallback failed: {exc}")

    payload = {
        "schema_version": "chem-json-v1",
        "source": {
            "file": str(output_file),
            "name": output_file.name,
            "program": program,
            "cclib_version": cclib.__version__,
            "parser_strategy": "cclib_primary_custom_parser_fallback",
        },
        "status": {
            "success": True,
            "warnings": warnings,
        },
        "molecule": {
            "n_atoms": len(atoms),
            "atoms": atoms,
        },
        "properties": {
            "energies": energies,
            "dipole": dipole,
            "polarizability": pol,
            "orbitals": {
                "homo_lumo_gap_ev": gap_ev,
            },
        },
        "charges": {
            "nbo": nbo,
            "hirshfeld": hirsh,
            "cm5": cm5,
            "cclib_atomcharges": make_json_safe(cclib_attrs.get("atomcharges", {})),
        },
        "vibrations": {
            "modes": vibrations,
        },
        "orca_summary": orca_summary if program == "orca" else {},
        "gaussian_custom": gaussian_custom if program == "gaussian" else {},
        "current_feather_compatibility": extract_feather_df_sections(feather_df),
    }

    if include_full_cclib:
        payload["cclib"] = {
            "metadata": make_json_safe(cclib_attrs.get("metadata", {})),
            "attributes": make_json_safe(cclib_attrs),
        }
    else:
        payload["cclib"] = {
            "metadata": make_json_safe(cclib_attrs.get("metadata", {})),
            "selected_attributes": make_json_safe({
                key: cclib_attrs.get(key)
                for key in [
                    "atomnos", "atomcoords", "charge", "mult", "scfenergies",
                    "vibfreqs", "vibirs", "vibdisps", "homos", "moenergies",
                    "atomcharges", "moments",
                ]
                if key in cclib_attrs
            }),
        }

    return make_json_safe(payload)


def save_json(payload: Dict[str, Any], json_file: Union[str, Path]) -> Path:
    json_file = Path(json_file).expanduser().resolve()
    json_file.parent.mkdir(parents=True, exist_ok=True)
    with json_file.open("w", encoding="utf-8") as f:
        json.dump(make_json_safe(payload), f, indent=2, ensure_ascii=False, allow_nan=False)
    return json_file


def _emit(log: Optional[Any], message: str) -> None:
    if log is not None:
        log.write(message)
    else:
        print(message)


# ---------------------------------------------------------------------------
# 7. Main conversion entry points
# ---------------------------------------------------------------------------

def convert_one_output(
    output_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    make_json: bool = True,
    make_feather: bool = True,
    include_full_cclib: bool = True,
    milo_layout: bool = False,
    log: Optional[Any] = None,
) -> Dict[str, Any]:
    output_file = Path(output_file).expanduser().resolve()
    if output_folder is None:
        output_folder = output_file.parent / "converted_outputs"
    else:
        output_folder = Path(output_folder).expanduser().resolve()

    json_dir = output_folder if milo_layout else output_folder / "json_files"
    feather_dir = output_folder / "feather_files"
    xyz_file = output_folder / f"{output_file.stem}.xyz" if milo_layout else None

    text = read_text_file(output_file)
    data, cclib_error = read_cclib_output(output_file)
    cclib_attrs = get_cclib_attributes(data)
    program = detect_program(output_file, text, cclib_attrs)

    warnings: List[str] = []
    if cclib_error:
        warnings.append(f"cclib parse warning: {cclib_error}")

    feather_df = build_feather_compatible_df(output_file, data, text, program, warnings)

    result = {
        "input_file": str(output_file),
        "program": program,
        "json_file": None,
        "xyz_file": None,
        "feather_file": None,
        "warnings": warnings,
    }

    if make_feather:
        feather_file = feather_dir / f"{output_file.stem}.feather"
        save_current_feather(feather_df, feather_file)
        result["feather_file"] = str(feather_file)

    if make_json:
        payload = build_structured_json(
            output_file=output_file,
            data=data,
            text=text,
            program=program,
            feather_df=feather_df,
            warnings=warnings,
            include_full_cclib=include_full_cclib,
        )
        json_file = json_dir / f"{output_file.stem}.json"
        save_json(payload, json_file)
        result["json_file"] = str(json_file)

        if milo_layout:
            atoms = payload.get("molecule", {}).get("atoms", [])
            if atoms:
                saved_xyz = save_xyz(atoms, xyz_file, comment=output_file.name)
                result["xyz_file"] = str(saved_xyz)
            else:
                warnings.append(f"No atomic coordinates found for {output_file.name}; XYZ file was not written.")

    return result


def convert_folder(
    input_folder: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    patterns: Iterable[str] = ("*.log", "*.out", "*.orcaout"),
    recursive: bool = False,
    make_json: bool = True,
    make_feather: bool = True,
    include_full_cclib: bool = True,
    milo_layout: bool = False,
    log: Optional[Any] = None,
) -> List[Dict[str, Any]]:
    input_folder = Path(input_folder).expanduser().resolve()
    if output_folder is None:
        output_folder = input_folder / "converted_outputs"
    else:
        output_folder = Path(output_folder).expanduser().resolve()

    results = []
    seen = set()

    for pattern in patterns:
        iterator = input_folder.rglob(pattern) if recursive else input_folder.glob(pattern)
        for file_path in iterator:
            file_path = file_path.resolve()
            if not file_path.is_file() or file_path in seen:
                continue
            seen.add(file_path)

            try:
                result = convert_one_output(
                    output_file=file_path,
                    output_folder=output_folder,
                    make_json=make_json,
                    make_feather=make_feather,
                    include_full_cclib=include_full_cclib,
                    milo_layout=milo_layout,
                    log=log,
                )
                results.append(result)
                _emit(
                    log,
                    f"Converted {file_path.name}: "
                    f"JSON={bool(result['json_file'])}, "
                    f"XYZ={bool(result.get('xyz_file'))}, "
                    f"Feather={bool(result['feather_file'])}, "
                    f"program={result['program']}"
                )
            except Exception as exc:
                result = {
                    "input_file": str(file_path),
                    "program": "unknown",
                    "json_file": None,
                    "xyz_file": None,
                    "feather_file": None,
                    "warnings": [str(exc)],
                }
                results.append(result)
                _emit(log, f"Failed {file_path.name}: {exc}")

    return results


def convert_parent_mode(
    parent_dir: Union[str, Path],
    com_dirname: str = "com",
    collected_dir: str = "converted_collected",
    patterns: Iterable[str] = ("*.log", "*.out", "*.orcaout"),
    recursive: bool = False,
    overwrite: bool = False,
    include_full_cclib: bool = True,
    log: Optional[Any] = None,
) -> List[Dict[str, Any]]:
    parent_dir = Path(parent_dir).expanduser().resolve()
    collected_root = parent_dir / collected_dir
    collected_root.mkdir(parents=True, exist_ok=True)

    all_results = []

    for entry in sorted(parent_dir.iterdir()):
        if not entry.is_dir():
            continue

        com_path = entry / com_dirname
        if not com_path.is_dir():
            continue

        dest = collected_root / entry.name
        if dest.exists() and overwrite:
            shutil.rmtree(dest)
        dest.mkdir(parents=True, exist_ok=True)

        _emit(log, f"[PARENT MODE] Processing: {com_path}")
        results = convert_folder(
            input_folder=com_path,
            output_folder=dest,
            patterns=patterns,
            recursive=recursive,
            make_json=True,
            make_feather=True,
            include_full_cclib=include_full_cclib,
            log=log,
        )
        all_results.extend(results)

    _emit(log, f"[PARENT MODE] Collected outputs into: {collected_root}")
    return all_results


# ---------------------------------------------------------------------------
# 8. Command-line interface
# ---------------------------------------------------------------------------

PROGRAM_PATTERNS = {
    "orca": ("*.out", "*.orcaout"),
    "gaussian": ("*.log",),
    "both": ("*.log", "*.out", "*.orcaout"),
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract structured JSON (and optional Feather tables) from Gaussian and/or ORCA output files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input_path",
        help="Path to a single output file (.log/.out/.orcaout), or a folder containing many.",
    )
    parser.add_argument(
        "-p", "--program",
        choices=["orca", "gaussian", "both"],
        default="both",
        help="When input_path is a folder, which file types to include (default: both). "
             "Ignored for a single file (program is auto-detected from its content).",
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_folder",
        default=None,
        help="Output folder for json_files/ and feather_files/ (default: <input>/converted_outputs).",
    )
    parser.add_argument(
        "-r", "--recursive",
        action="store_true",
        help="When input_path is a folder, also search subfolders.",
    )
    parser.add_argument(
        "--no-json",
        action="store_true",
        help="Skip writing JSON output.",
    )
    parser.add_argument(
        "--no-feather",
        action="store_true",
        help="Skip writing Feather output (avoids needing pyarrow).",
    )
    parser.add_argument(
        "--minimal-cclib",
        action="store_true",
        help="Store only a selected subset of cclib attributes in the JSON instead of the full attribute dump (smaller files).",
    )
    parser.add_argument(
        "--parent-mode",
        action="store_true",
        help="Old parent/com layout: input_path/<molecule>/com/*.out|*.log -> input_path/converted_collected/<molecule>/...",
    )
    parser.add_argument(
        "--com-dirname",
        default="com",
        help="Subfolder name to look for in --parent-mode (default: com).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    input_path = Path(args.input_path).expanduser().resolve()
    if not input_path.exists():
        print(f"ERROR: input path does not exist: {input_path}", file=sys.stderr)
        return 1

    make_json = not args.no_json
    make_feather = not args.no_feather
    include_full_cclib = not args.minimal_cclib

    if args.parent_mode:
        if not input_path.is_dir():
            print("ERROR: --parent-mode requires input_path to be a directory.", file=sys.stderr)
            return 1
        patterns = PROGRAM_PATTERNS[args.program]
        convert_parent_mode(
            parent_dir=input_path,
            com_dirname=args.com_dirname,
            patterns=patterns,
            recursive=args.recursive,
            include_full_cclib=include_full_cclib,
        )
        return 0

    if input_path.is_file():
        result = convert_one_output(
            output_file=input_path,
            output_folder=args.output_folder,
            make_json=make_json,
            make_feather=make_feather,
            include_full_cclib=include_full_cclib,
        )
        print(json.dumps(result, indent=2))
        if result["warnings"]:
            print(f"\n{len(result['warnings'])} warning(s):", file=sys.stderr)
            for w in result["warnings"]:
                print(f"  - {w}", file=sys.stderr)
        return 0

    if input_path.is_dir():
        patterns = PROGRAM_PATTERNS[args.program]
        results = convert_folder(
            input_folder=input_path,
            output_folder=args.output_folder,
            patterns=patterns,
            recursive=args.recursive,
            make_json=make_json,
            make_feather=make_feather,
            include_full_cclib=include_full_cclib,
        )
        n_ok = sum(1 for r in results if r.get("json_file") or r.get("feather_file"))
        n_fail = len(results) - n_ok
        print(f"\nDone. {n_ok} converted, {n_fail} failed, out of {len(results)} file(s).")
        return 0 if n_fail == 0 else 1

    print(f"ERROR: input path is neither a file nor a directory: {input_path}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
