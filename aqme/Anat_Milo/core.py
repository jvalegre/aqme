import os
import time
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import yaml

from .engines import (
    compute_charge_difference,
    compute_cube_sterimol,
    compute_rigid_sterimol,
    compute_rotated_dipole,
    compute_standard_geometry,
)
from .extract_chem_json import PROGRAM_PATTERNS, convert_folder
from .loader import PairedMoleculeLoader


def _emit(log: Optional[Any], message: str) -> None:
    if log is not None:
        log.write(message)
    else:
        print(message)


def prepare_milo_folder(
    outputs_dir: str,
    recursive: bool = False,
    program: str = "both",
    log: Optional[Any] = None,
) -> Path:
    """
    Convert Gaussian/ORCA output files into paired JSON + XYZ files inside a MILO folder.
    """
    source_dir = Path(outputs_dir).expanduser().resolve()
    if not source_dir.is_dir():
        message = f"Output folder target path is unreachable: {source_dir}"
        _emit(log, f"x  {message}")
        raise NotADirectoryError(message)

    milo_dir = source_dir.parent / "MILO"
    milo_dir.mkdir(parents=True, exist_ok=True)

    patterns = PROGRAM_PATTERNS.get(program, PROGRAM_PATTERNS["both"])
    output_files = []
    for pattern in patterns:
        iterator = source_dir.rglob(pattern) if recursive else source_dir.glob(pattern)
        for file_path in iterator:
            if file_path.is_file():
                output_files.append(file_path)

    if not output_files:
        message = (
            f"No Gaussian/ORCA output files were found in {source_dir}. "
            f"Expected patterns: {', '.join(patterns)}"
        )
        _emit(log, f"x  {message}")
        raise FileNotFoundError(message)

    _emit(
        log,
        f"\nStarting MILO with {len(output_files)} job(s) "
        "(Gaussian and ORCA outputs will be converted to JSON and XYZ)\n",
    )
    _emit(log, f"o  Creating MILO folder at {milo_dir}")
    results = convert_folder(
        input_folder=source_dir,
        output_folder=milo_dir,
        patterns=patterns,
        recursive=recursive,
        make_json=True,
        make_feather=False,
        include_full_cclib=True,
        milo_layout=True,
        log=log,
    )
    n_ok = sum(1 for r in results if r.get("json_file") or r.get("xyz_file"))
    n_fail = len(results) - n_ok
    _emit(log, f"\nDone. {n_ok} converted, {n_fail} failed, out of {len(results)} file(s).\n")

    return milo_dir


def _build_descriptor_matrix(
    yaml_file: str,
    data_dir: str,
    output_csv: Optional[str] = None,
    log: Optional[Any] = None,
) -> pd.DataFrame:
    if not os.path.exists(yaml_file):
        message = f"YAML Configuration file missing: {yaml_file}"
        _emit(log, f"x  {message}")
        raise FileNotFoundError(message)
    if not os.path.isdir(data_dir):
        message = f"Data folder target path is unreachable: {data_dir}"
        _emit(log, f"x  {message}")
        raise NotADirectoryError(message)

    with open(yaml_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    is_one_based = config.get("atom_indexing", "one_based") == "one_based"
    desc_config = config.get("descriptors", {})

    compiled_dataset = []
    all_files = os.listdir(data_dir)
    molecule_bases = {os.path.splitext(f)[0] for f in all_files if f.endswith(".xyz")}

    _emit(
        log,
        f"\nStarting MILO descriptor calculation with {len(molecule_bases)} structure(s)\n",
    )

    if not molecule_bases:
        message = f"No .xyz files were found in the MILO folder: {data_dir}"
        _emit(log, f"x  {message}")
        raise FileNotFoundError(message)

    for base in sorted(molecule_bases):
        xyz_path = os.path.join(data_dir, f"{base}.xyz")
        json_path = os.path.join(data_dir, f"{base}.json")
        cube_path = os.path.join(data_dir, f"{base}.cube")

        if not os.path.exists(json_path):
            continue

        mol = PairedMoleculeLoader(base, xyz_path, json_path, cube_path=cube_path)

        name_parts = base.split("_")
        code_name = "_".join(name_parts[:2]) if len(name_parts) >= 2 else base
        row_metrics = {"code_name": code_name}

        row_metrics.update(compute_standard_geometry(mol, desc_config, is_one_based))
        row_metrics.update(compute_rigid_sterimol(mol, desc_config, is_one_based))
        row_metrics.update(compute_cube_sterimol(mol, desc_config, is_one_based))
        row_metrics.update(compute_rotated_dipole(mol, desc_config, is_one_based))
        row_metrics.update(compute_charge_difference(mol, desc_config, is_one_based))

        if "global_properties" in desc_config and desc_config["global_properties"].get("enabled", False):
            for key in desc_config["global_properties"].get("keys", []):
                val = mol.get_prop(key)
                if val is not None:
                    row_metrics[key.lower().replace(" ", "_")] = float(val) if isinstance(val, (int, float)) else val

        if "atomic_properties" in desc_config and desc_config["atomic_properties"].get("enabled", False):
            for lbl, defs in desc_config["atomic_properties"].get("definitions", {}).items():
                atom_idx = defs.get("atom") - 1 if is_one_based else defs.get("atom")
                array_data = mol.get_prop(defs.get("key"))
                if array_data and 0 <= atom_idx < len(array_data):
                    row_metrics[lbl] = array_data[atom_idx]

        compiled_dataset.append(row_metrics)

    if not compiled_dataset:
        message = f"No matching .xyz/.json file pairs were found in: {data_dir}"
        _emit(log, f"x  {message}")
        raise FileNotFoundError(message)

    df_final = pd.DataFrame(compiled_dataset)

    if output_csv:
        output_csv_path = Path(output_csv).expanduser().resolve()
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        df_final.to_csv(output_csv_path, index=False)
        _emit(log, f"[Anat_Milo] Matrix processing complete. Outputs stored inside: {output_csv_path}")
    else:
        _emit(log, "[Anat_Milo] Matrix processing complete.")

    return df_final


def run_anat_milo_workflow(
    yaml_file: Optional[str] = None,
    outputs_dir: Optional[str] = None,
    data_dir: Optional[str] = None,
    output_csv: Optional[str] = None,
    recursive: bool = False,
    program: str = "both",
    log: Optional[Any] = None,
):
    """
    Run the MILO workflow.

    Stage 1:
        Convert Gaussian/ORCA outputs into paired JSON + XYZ files in a MILO folder.
    Stage 2:
        If a YAML file is supplied, compute descriptors from the generated JSON/XYZ pairs
        and export the CSV.
    """
    milo_dir: Optional[Path] = None

    if outputs_dir:
        milo_dir = prepare_milo_folder(outputs_dir, recursive=recursive, program=program, log=log)

    if data_dir and milo_dir is None:
        milo_dir = Path(data_dir).expanduser().resolve()
        if not milo_dir.is_dir():
            message = f"Data folder target path is unreachable: {milo_dir}"
            _emit(log, f"x  {message}")
            raise NotADirectoryError(message)

    if milo_dir is None:
        message = "Provide either outputs_dir or data_dir so MILO has something to process."
        _emit(log, f"x  {message}")
        raise ValueError(message)

    if output_csv is None and yaml_file:
        output_csv = str(milo_dir / "aqme_milo_features.csv")

    if not yaml_file:
        _emit(log, f"[Anat_Milo] MILO JSON and XYZ ready in: {milo_dir}")
        return milo_dir

    return _build_descriptor_matrix(
        yaml_file=yaml_file,
        data_dir=str(milo_dir),
        output_csv=output_csv,
        log=log,
    )
