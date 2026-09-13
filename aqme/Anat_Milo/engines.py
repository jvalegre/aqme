import numpy as np
import pandas as pd
# import igraph as ig
# from typing import Sequence, Tuple, List

# --- RADIUS REFERENCE DATABASES ---
BONDI_RADII = {'H': 1.10, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.83, 'I': 1.98}
CPK_RADII = {'C': 1.50, 'C3': 1.60, 'C6/N6': 1.70, 'H': 1.00, 'N': 1.50, 'N4': 1.45, 'O': 1.35, 'O2': 1.35, 'P': 1.40, 'S': 1.70, 'S1': 1.00, 'F': 1.35, 'Cl': 1.80, 'Br': 1.95, 'I': 2.15}

# --- NATIVE UTILITY MATH MATRICES ---
def adjust_indices(element):
    """Convert user-facing 1-based atom indices into 0-based indices."""
    if isinstance(element, (list, tuple, np.ndarray)):
        return [adjust_indices(sub) for sub in element]
    if isinstance(element, (int, np.integer)):
        return int(element) - 1
    if isinstance(element, float) and element.is_integer():
        return int(element) - 1
    raise ValueError(f"Unsupported index element: {element}")

def calc_angle(p1: np.ndarray, p2: np.ndarray, degrees: bool = False) -> float:
    """Calculate the angle between two vectors."""
    denom = np.linalg.norm(p1) * np.linalg.norm(p2)
    if denom == 0: return 0.0
    cosine = np.clip(np.dot(p1, p2) / denom, -1.0, 1.0)
    theta = np.arccos(cosine)
    return float(np.degrees(theta) if degrees else theta)

def generate_circle(center_x, center_y, radius, n_points=100):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))

# --- GEOMETRICAL FEATURE EXTRACTION ENGINE ---
def compute_standard_geometry(mol, config_desc, is_one_based) -> dict:
    """
    Batch extract standard geometric parameters (bond lengths, bond angles, 
    and dihedral angles) based on targets enabled in the configuration file.
    """
    results = {}
    coords = mol.coordinates
    
    if "bond_length" in config_desc:
        for lbl, atoms in config_desc["bond_length"].get("definitions", {}).items():
            i1 = atoms[0] - 1 if is_one_based else atoms[0]
            i2 = atoms[1] - 1 if is_one_based else atoms[1]
            results[f"{lbl}_length"] = float(np.linalg.norm(coords[i1] - coords[i2]))
            
    if "bond_angle" in config_desc:
        for lbl, atoms in config_desc["bond_angle"].get("definitions", {}).items():
            i1, i2, i3 = [a - 1 if is_one_based else a for a in atoms[:3]]
            v1, v2 = coords[i1] - coords[i2], coords[i3] - coords[i2]
            results[f"{lbl}_angle"] = calc_angle(v1, v2, degrees=True)
            
    if "dihedral_angle" in config_desc:
        for lbl, atoms in config_desc["dihedral_angle"].get("definitions", {}).items():
            p0, p1, p2, p3 = [coords[a - 1 if is_one_based else a] for a in atoms[:4]]
            b0, b1, b2 = -1.0 * (p1 - p0), p2 - p1, p3 - p2
            b1_n = b1 / np.linalg.norm(b1)
            v = b0 - np.dot(b0, b1_n) * b1_n
            w = b2 - np.dot(b2, b1_n) * b1_n
            results[f"{lbl}_dihedral"] = float(np.degrees(np.arctan2(np.dot(np.cross(b1_n, v), w), np.dot(v, w))))
            
    return results

# --- COORDINATE AND ROTATION CALCULATORS ---
def calc_new_base_atoms(coordinates_array, atom_indices):
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices[0] is list, compute the new origin as the middle of the first atoms.
    """
    if isinstance(atom_indices[0], list):
        new_origin = np.mean(coordinates_array[atom_indices[0]], axis=0)
    else:
        new_origin = coordinates_array[atom_indices[0]]
    new_y = (coordinates_array[atom_indices[1]] - new_origin) / np.linalg.norm(coordinates_array[atom_indices[1]] - new_origin)
    coplane = (coordinates_array[atom_indices[2]] - new_origin) / np.linalg.norm(coordinates_array[atom_indices[2]] - new_origin + 1e-8)
    return new_origin, new_y, coplane

def calc_basis_vector(origin, y, coplane):
    """
    Calculate the new basis vector.
    
    Parameters
    ----------
    origin : array-like
        The origin of the new basis.
    y : array-like
        The new basis's y direction.
    coplane : array-like
        A vector coplanar with the new y direction.
    
    Returns
    -------
    new_basis : np.array
        The computed new basis matrix.
    """
    coef_mat = np.vstack([coplane, y, np.cross(coplane, y)])
    angle_new_y_coplane = calc_angle(coplane, y)
    result_vector = [np.cos(angle_new_y_coplane - (np.pi/2)), 0, 0]
    new_x, _, _, _ = np.linalg.lstsq(coef_mat, result_vector, rcond=None)
    return np.vstack([new_x, y, np.cross(new_x, y)])

def preform_coordination_transformation(xyz_df, indices):
    """
    Perform a coordination transformation on the xyz DataFrame.
    
    Parameters
    ----------
    xyz_df : pd.DataFrame
        DataFrame containing columns 'x', 'y', 'z'.
    indices : array-like, optional
        Atom indices to use for the new basis. If None, default indices [1,2,3] are used.
    
    Returns
    -------
    xyz_copy : pd.DataFrame
        DataFrame with transformed coordinates.
    """
    xyz_copy = xyz_df.copy()
    coords = np.array(xyz_copy[['x', 'y', 'z']].values)
    origin, y_vec, coplane = calc_new_base_atoms(coords, indices)
    basis = calc_basis_vector(origin, y_vec, coplane)
    transformed = np.array([np.dot(basis, row - origin) for row in coords])
    xyz_copy[['x', 'y', 'z']] = transformed
    return xyz_copy

# --- RIGID VDW STERIMOL EXTRACTION ENGINE ---
def b1s_for_loop_function(degree_list, plane):
    """
    a function that gets a plane transform it and calculate the b1s for each degree.
    checks if the plane is in the x or z axis and calculates the b1s accordingly.
    Parameters:
    ----------
    extended_df : pd.DataFrame
    b1s : list
    b1s_loc : list
    degree_list : list
    plane : np.array
    """
    results = []
    for degree in degree_list:
        theta = np.deg2rad(degree)
        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        transformed_plane = plane @ rot_matrix.T
        
        max_x, min_x = np.max(transformed_plane[:, 0]), np.min(transformed_plane[:, 0])
        max_y, min_y = np.max(transformed_plane[:, 1]), np.min(transformed_plane[:, 1])
        avs = np.abs([max_x, min_x, max_y, min_y])
        
        results.append({
            'degree': degree, 'B1': np.min(avs), 'plane': transformed_plane,
            'b1_normal': np.array([np.cos(theta), np.sin(theta)]) if np.argmin(avs) in (0, 1) else np.array([-np.sin(theta), np.cos(theta)])
        })
    return pd.DataFrame(results)

def direction_atoms_for_sterimol(bonds_df, base_atoms) -> list:
    """
    Return [origin, direction, third] for Sterimol coordinate setup.
    - bonds_df: DataFrame with two columns (0,1), undirected bonds between integer atom indices.
    - base_atoms: [origin, direction] or [origin, direction, third].
       If third == origin, choose a neighbor of 'direction' different from origin.
    Raises ValueError with a clear message if it can't determine 'third'.
    """
    if isinstance(base_atoms[0], list): base_atoms = base_atoms[0]
    origin, direction = int(base_atoms[0]), int(base_atoms[1])
    col0, col1 = bonds_df[0].astype(int), bonds_df[1].astype(int)
    nbr_dir = set(col1[col0 == direction]).union(set(col0[col1 == direction]))
    cand = [a for a in nbr_dir if a != origin]
    third = int(cand[0]) if cand else direction + 1
    return [origin, direction, third]

def compute_rigid_sterimol(mol, config_desc, is_one_based) -> dict:
    """
    Calculate B1 (minimum width), B5 (maximum width), L (length), and loc_B5 steric parameters 
    using the native coarse and fine rotational scanning loops on circle-expanded atom vectors.
    """
    if "sterimol" not in config_desc: return {}
    cfg = config_desc["sterimol"]
    radii_map = CPK_RADII if cfg.get("radii", "CPK") == "CPK" else BONDI_RADII
    results = {}

    for lbl, defs in cfg.get("definitions", {}).items():
        atoms = defs.get("atoms", [])
        adj_atoms = [a - 1 if is_one_based else a for a in atoms[:2]]
        
        bonds_direction = adjust_indices(direction_atoms_for_sterimol(mol.bonds_df, atoms))
        transformed_xyz = preform_coordination_transformation(mol.xyz_df, bonds_direction)
        
        transformed_xyz['radius'] = transformed_xyz['atom'].map(radii_map).fillna(1.70)
        transformed_xyz['B5'] = transformed_xyz['radius'] + np.linalg.norm(transformed_xyz[['x', 'z']].values, axis=1)
        transformed_xyz['L'] = transformed_xyz['y'] + transformed_xyz['radius']
        
        circles = [generate_circle(row['x'], row['z'], row['radius']) for _, row in transformed_xyz.iterrows()]
        plane_xz = np.vstack(circles)
        
        # Coarse Sweep
        df_coarse = b1s_for_loop_function(list(range(18, 108, 5)), plane_xz)
        best_ang = int(df_coarse.loc[df_coarse["B1"].idxmin(), "degree"])
        
        # Fine Sweep Window
        df_fine = b1s_for_loop_function(list(range(best_ang - 5, best_ang + 6)), plane_xz)
        best_row = df_fine.loc[df_fine["B1"].idxmin()]
        
        results[f"{lbl}_B1"] = float(best_row["B1"])
        results[f"{lbl}_B5"] = float(np.max(transformed_xyz["B5"]))
        results[f"{lbl}_L"] = float(np.max(transformed_xyz["L"]))
        results[f"{lbl}_loc_B5"] = float(transformed_xyz['y'].iloc[transformed_xyz['B5'].idxmax()])
        
    return results

# --- ELECTRON DENSITY GRID CUBE STERIMOL ---
def compute_cube_sterimol(mol, config_desc, is_one_based) -> dict:
    """
    Process Gaussian electron density grid data (.cube files) to extract B1, B5, and L parameters 
    by mapping coordinates matching a defined electronic isovalue (default 0.003 a.u.) over rotated planes.
    """
    if "cube_sterimol" not in config_desc or not mol.cube_path: return {}
    cfg = config_desc["cube_sterimol"]
    isovalue = cfg.get("isovalue", 0.003)
    results = {}

    # Standalone parse for electron density blocks
    with open(mol.cube_path, 'r') as f:
        cube_lines = f.readlines()
    
    n_atoms = int(cube_lines[2].split()[0])
    structure_lines = [cube_lines[2:][i] for i in range(4 + n_atoms)]
    
    # Isolate boundary origins and grid steps
    x_origin, y_origin, z_origin = [float(x) for x in structure_lines[0].split()[1:4]]
    x_step, x_size = int(structure_lines[1].split()[0]), float(structure_lines[1].split()[1])
    y_step, y_size = int(structure_lines[2].split()[0]), float(structure_lines[2].split()[2])
    z_step, z_size = int(structure_lines[3].split()[0]), float(structure_lines[3].split()[3])

    # Reconstruct raw points mapping to spatial grids
    density_values = []
    for line in cube_lines[(6 + n_atoms):]:
        density_values.extend([float(val) for val in line.split()])
    
    density_array = np.array(density_values).reshape(x_step, y_step, z_step)
    matching_points = np.argwhere((density_array >= isovalue) & (density_array <= isovalue * 1.1))
    
    if len(matching_points) == 0: return {}

    # Translate voxel locations to Cartesian Angstrom space
    grid_coords = []
    for pt in matching_points:
        cx = (x_origin + pt[0] * x_size) * 0.529177249
        cy = (y_origin + pt[1] * y_size) * 0.529177249
        cz = (z_origin + pt[2] * z_size) * 0.529177249
        grid_coords.append([cx, cy, cz])
    grid_coords = np.array(grid_coords)

    for lbl, defs in cfg.get("definitions", {}).items():
        atoms = defs.get("atoms", [])
        if len(atoms) < 2: continue
        o, d = atoms[0] - 1 if is_one_based else atoms[0], atoms[1] - 1 if is_one_based else atoms[1]
        
        axis_y = mol.coordinates[d] - mol.coordinates[o]
        axis_y /= np.linalg.norm(axis_y)
        
        shifted_grid = grid_coords - mol.coordinates[o]
        y_proj = np.dot(shifted_grid, axis_y)
        
        perp_components = shifted_grid - np.outer(y_proj, axis_y)
        perp_widths = np.linalg.norm(perp_components, axis=1)
        
        results[f"{lbl}_cube_L"] = float(np.max(y_proj))
        results[f"{lbl}_cube_B5"] = float(np.max(perp_widths))
        results[f"{lbl}_cube_B1"] = float(np.min(perp_widths))
        
    return results

# --- ELECTRONIC & ROTATED VECTOR PROPAGATOR ---
def compute_rotated_dipole(mol, config_desc, is_one_based) -> dict:
    """
    Transform Gaussian dipole(s) into a molecular frame defined by base_atoms_indices.
    
    Semantics (R-equivalent):
      - len==3: [origin_atom, y_atom, plane_atom]
      - len>=4: [origin_set..., y_atom, plane_atom]  (origin = centroid(origin_set))
      - [[origin_set], y_atom, plane_atom]
      - [[origin_set], [direction_set], [plane_set]]
    """
    if "rotated_dipole" not in config_desc: return {}
    global_dipole = mol.get_prop("Dipole moment")
    if not global_dipole or len(global_dipole) < 3: return {}
    
    results = {}
    dip_vec = np.array(global_dipole[:3], dtype=float)
    
    for lbl, defs in config_desc["rotated_dipole"].get("definitions", {}).items():
        atoms = defs.get("atoms", [])
        if not isinstance(atoms, list) or len(atoms) < 3: continue
        
        orig_set = adjust_indices(atoms[0] if isinstance(atoms[0], list) else [atoms[0]])
        y_atom = adjust_indices(atoms[1])
        p_atom = adjust_indices(atoms[2])
        
        origin_pt = mol.coordinates[orig_set].mean(axis=0)
        basis = calc_basis_vector(origin_pt, mol.coordinates[y_atom] - origin_pt, mol.coordinates[p_atom] - origin_pt)
        rotated = np.dot(basis, dip_vec)
        
        results[f"{lbl}_dip_x"] = float(rotated[0])
        results[f"{lbl}_dip_y"] = float(rotated[1])
        results[f"{lbl}_dip_z"] = float(rotated[2])
    return results

def compute_charge_difference(mol, config_desc, is_one_based) -> dict:
    """
    Returns a DataFrame with the charge differences for specified atom pairs for the requested type(s).
    The difference computes atom_i - atom_j from the parsed partial charge dictionary arrays.
    """
    if "charge_difference" not in config_desc: return {}
    charges = mol.get_prop("Partial charge")
    if not charges: return {}
    
    results = {}
    for lbl, defs in config_desc["charge_difference"].get("definitions", {}).items():
        atoms = defs.get("atoms", [])
        i1 = atoms[0] - 1 if is_one_based else atoms[0]
        i2 = atoms[1] - 1 if is_one_based else atoms[1]
        
        if max(i1, i2) < len(charges):
            results[f"{lbl}_charge_delta"] = float(charges[i1] - charges[i2])
    return results