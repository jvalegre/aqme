import os
# import re
import yaml
import pandas as pd
from .loader import PairedMoleculeLoader
from .engines import (
    compute_standard_geometry,
    compute_rigid_sterimol,
    compute_cube_sterimol,
    compute_rotated_dipole,
    compute_charge_difference
)

def run_anat_milo_workflow(yaml_file: str, data_dir: str, output_csv: str = "aqme_milo_features.csv"):
    """
    Main execution orchestrator linking the config matrix with file storage matrices.
    
    Parameters
    ----------
    yaml_file : str
        Path to the external YAML configuration file.
    data_dir : str
        Directory path containing paired .xyz, .json, and optional .cube files.
    output_csv : str, optional
        Output CSV filename for the generated descriptor matrix.
    """
    if not os.path.exists(yaml_file):
        raise FileNotFoundError(f"YAML Configuration file missing: {yaml_file}")
    if not os.path.isdir(data_dir):
        raise NotADirectoryError(f"Data folder target path is unreachable: {data_dir}")

    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f) or {}

    is_one_based = config.get("atom_indexing", "one_based") == "one_based"
    desc_config = config.get("descriptors", {})
    
    compiled_dataset = []
    
    # Isolate shared unique molecular filenames inside the data directory
    all_files = os.listdir(data_dir)
    molecule_bases = {os.path.splitext(f)[0] for f in all_files if f.endswith('.xyz')}
    
    for base in sorted(molecule_bases):
        xyz_path = os.path.join(data_dir, f"{base}.xyz")
        json_path = os.path.join(data_dir, f"{base}.json")
        cube_path = os.path.join(data_dir, f"{base}.cube")
        
        if not os.path.exists(json_path):
            continue # Skipped dynamically if computational JSON file pair is missing
            
        mol = PairedMoleculeLoader(base, xyz_path, json_path, cube_path=cube_path)
        
        # Parse the code_name column content until the second underscore (excluding it)
        # Example: 'mol_9_rdkit' split by '_' picks ['mol', '9'] -> 'mol_9'
        name_parts = base.split('_')
        code_name = "_".join(name_parts[:2]) if len(name_parts) >= 2 else base
        
        row_metrics = {"code_name": code_name}
        
        # Sequentially map operational parameters driven by config switches
        row_metrics.update(compute_standard_geometry(mol, desc_config, is_one_based))
        row_metrics.update(compute_rigid_sterimol(mol, desc_config, is_one_based))
        row_metrics.update(compute_cube_sterimol(mol, desc_config, is_one_based))
        row_metrics.update(compute_rotated_dipole(mol, desc_config, is_one_based))
        row_metrics.update(compute_charge_difference(mol, desc_config, is_one_based))
        
        # Pull direct raw scalar parameters out of the computational JSON roots
        if "global_properties" in desc_config and desc_config["global_properties"].get("enabled", False):
            for key in desc_config["global_properties"].get("keys", []):
                val = mol.get_prop(key)
                if val is not None:
                    row_metrics[key.lower().replace(" ", "_")] = float(val) if isinstance(val, (int, float)) else val
                    
        # Pull atom index arrays data fields
        if "atomic_properties" in desc_config and desc_config["atomic_properties"].get("enabled", False):
            for lbl, defs in desc_config["atomic_properties"].get("definitions", {}).items():
                atom_idx = defs.get("atom") - 1 if is_one_based else defs.get("atom")
                array_data = mol.get_prop(defs.get("key"))
                if array_data and 0 <= atom_idx < len(array_data):
                    row_metrics[lbl] = array_data[atom_idx]

        compiled_dataset.append(row_metrics)
        
    # Compile results and export directly to a tabular CSV format
    if compiled_dataset:
        df_final = pd.DataFrame(compiled_dataset)
        df_final.to_csv(output_csv, index=False)
        print(f"[Anat_Milo] Matrix processing complete. Outputs stored inside: {output_csv}")
    else:
        print("[Anat_Milo] Pipeline aborted. No matching structural file pairs detected.")