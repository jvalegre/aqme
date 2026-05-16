from typing import List, Optional
import numpy as np
import numpy.typing as npt
import pandas as pd


try:
    from ...utils.help_functions import *
    from ...utils.visualize import *
except ImportError:
    from utils.help_functions import * 
    from utils.visualize import *


def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices[0] is list, compute the new origin as the middle of the first atoms.
    """
   
    if isinstance(atom_indices[0], list):
        new_origin=np.mean(coordinates_array[atom_indices[0]], axis=0)
    else:
        new_origin=coordinates_array[atom_indices[0]]
    new_y=(coordinates_array[atom_indices[1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[1]]-new_origin))
    coplane=((coordinates_array[atom_indices[2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[2]]-new_origin)+0.00000001))

    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

from numpy.typing import ArrayLike

def calc_basis_vector(origin, y: ArrayLike, coplane: ArrayLike):
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
    coef_mat = np_cross_and_vstack(coplane, y)
    angle_new_y_coplane = calc_angle(coplane, y)
    cop_ang_x = angle_new_y_coplane - (np.pi/2)
    result_vector = [np.cos(cop_ang_x), 0, 0]
    new_x, _, _, _ = np.linalg.lstsq(coef_mat, result_vector, rcond=None)
    new_basis = np_cross_and_vstack(new_x, y)
 
    return new_basis



def calc_npa_charges(coordinates_array: npt.ArrayLike,charge_array: npt.ArrayLike):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charge_array: np.array
        array of npa charges
    base_atoms_indices:list
        3/4 atom indices for coordinates transformation
        
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.       
    Returns:
    -------
    dipole_df=calc_npa_charges(coordinates_array,charges,base_atoms_indices,sub_atoms)
    Output:
    dipole_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """

    
    dipole_xyz = np.vstack([(row[0] * row[1])for row in
                            list(zip(coordinates_array, charge_array))])
    dipole_vector=np.sum(dipole_xyz,axis=0)
    array_dipole=np.hstack([dipole_vector,np.linalg.norm(dipole_vector)])

    dipole_df=pd.DataFrame(array_dipole,index=XYZConstants.DIPOLE_COLUMNS.value).T
 
    return dipole_df

# def calc_dipole_gaussian(coordinates_array, gauss_dipole_array, base_atoms_indices, origin):
#     """
#     A function that receives coordinates and gaussian dipole, transforms the coordinates
#     by the new base atoms and calculates the dipole in each axis.
    
#     Parameters
#     ----------
#     coordinates_array: np.array
#         Contains x, y, z atom coordinates
#     gauss_dipole_array: np.array
#         Array of Gaussian dipole values
#     base_atoms_indices: list
#         3/4 atom indices for coordinates transformation or empty list to skip transformation
#     origin: list or None
#         Atom indices used to calculate the origin
        
#     Returns
#     -------
#     dipole_df: pd.DataFrame
#         DataFrame containing dipole values in x, y, z directions and total magnitude
#     """
    
#     if not base_atoms_indices:
#         # If indices are empty, return the gauss dipole array as is in a dataframe
#         if isinstance(gauss_dipole_array, pd.Series) or (isinstance(gauss_dipole_array, np.ndarray) and gauss_dipole_array.ndim == 1):
#             gauss_dipole_array = np.expand_dims(gauss_dipole_array, axis=0)
#         return pd.DataFrame(gauss_dipole_array, columns=['dipole_x', 'dipole_y', 'dipole_z', 'total'])

#     indices=adjust_indices(base_atoms_indices)
#     if origin is None:
#         new_origin = coordinates_array[indices[0]]
#     else:
#         orig_indices = adjust_indices(origin)
#         new_origin = np.mean(coordinates_array[orig_indices], axis=0)

#     # Recenter the coordinates array using the new origin
#     recentered_coords = coordinates_array - new_origin
#     basis_vector=calc_basis_vector(*calc_new_base_atoms(recentered_coords, indices))

#     gauss_dipole_array = [np.concatenate((np.matmul(basis_vector, gauss_dipole_array[0, 0:3]), [gauss_dipole_array[0, 3]]))]

#     if isinstance(gauss_dipole_array, pd.Series) or (isinstance(gauss_dipole_array, np.ndarray) and gauss_dipole_array.ndim == 1):
#         gauss_dipole_array = np.expand_dims(gauss_dipole_array, axis=0)

#     dipole_df=pd.DataFrame(gauss_dipole_array,columns=['dipole_x','dipole_y','dipole_z','total'])
    
#     return dipole_df



def calc_dipole_gaussian(
    coordinates_array: np.ndarray,
    gauss_dipole_array,
    base_atoms_indices: Sequence,
) -> pd.DataFrame:
    """
    Transform Gaussian dipole(s) into a molecular frame defined by base_atoms_indices.
    
    Semantics (R-equivalent):
      - len==3: [origin_atom, y_atom, plane_atom]
      - len>=4: [origin_set..., y_atom, plane_atom]  (origin = centroid(origin_set)
      - or: [[origin_set], y_atom, plane_atom]

    Parameters
    ----------
    coordinates_array : (N,3) array
        Atom coordinates.
    gauss_dipole_array : array-like
        Dipole as [dx, dy, dz, total] or [[...], ...]. If 'total' missing, it is computed.
    base_atoms_indices : sequence
        Indices defining origin/y/plane per semantics above. Supports 1- or 0-based.

    Returns
    -------
    pd.DataFrame with columns ['dipole_x','dipole_y','dipole_z','total']
    """

    # -------- helpers --------
    def _normalize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        n = np.linalg.norm(v)
        return v * 0.0 if n < eps else (v / n)

    def _to0(idx_list: Sequence[int]) -> np.ndarray:
        """Convert maybe-1-based indices to 0-based (heuristic: if any 0 present, assume already 0-based)."""
        idx = np.asarray(idx_list, dtype=int)
        return idx if (idx == 0).any() else (idx - 1)

    def _parse_base(group: Sequence) -> Tuple[np.ndarray, int, int]:
        """
        Returns (origin_set_0based, y_idx_0based, plane_idx_0based).
        Accepts [[o...], y, plane] OR [o, y, plane] OR [o..., y, plane].
        """
        if not isinstance(group, (list, tuple)) or len(group) < 3:
            raise ValueError("base_atoms_indices must be a sequence with at least 3 entries.")

        if isinstance(group[0], (list, tuple)):  # [[o...], y, plane]
            origin_set = _to0(group[0])
            y_idx = int(group[1]); plane_idx = int(group[2])
        elif len(group) >= 4:                    # [o..., y, plane]
            origin_set = _to0(group[:-2])
            y_idx = int(group[-2]); plane_idx = int(group[-1])
        else:                                    # [o, y, plane]
            origin_set = _to0([group[0]])
            y_idx = int(group[1]); plane_idx = int(group[2])

        y_idx0, plane_idx0 = _to0([y_idx])[0], _to0([plane_idx])[0]
        return origin_set, y_idx0, plane_idx0

    def _build_basis(coords: np.ndarray, origin_pt: np.ndarray, y_idx: int, plane_idx: int) -> np.ndarray:
        """Rows are basis axes (X, Y, Z) expressed in the original coordinates."""
        new_y = _normalize(coords[y_idx] - origin_pt)
        coplane = _normalize(coords[plane_idx] - origin_pt)

        # X = projection of 'coplane' onto plane ⟂ Y (with robust fallback)
        x_raw = coplane - np.dot(coplane, new_y) * new_y
        if np.linalg.norm(x_raw) < 1e-12:
            fallback = np.array([1.0, 0.0, 0.0]) if abs(new_y[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            x_raw = np.cross(fallback, new_y)
        new_x = _normalize(x_raw)

        new_z = _normalize(np.cross(new_x, new_y))
        basis = np.vstack([new_x, new_y, new_z])

        # Ensure right-handed basis
        if np.linalg.det(basis) < 0:
            new_x = -new_x
            basis = np.vstack([new_x, new_y, new_z])

        return basis

    # -------- inputs & shapes --------
    coords = np.asarray(coordinates_array, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coordinates_array must be shape (N, 3).")

    g = np.asarray(gauss_dipole_array, dtype=float)
    if g.ndim == 1:
        g = g.reshape(1, -1)
    if g.shape[1] < 3:
        raise ValueError("gauss_dipole_array must have at least 3 components per row (dx, dy, dz).")

    # No transform requested → return as-is (pad/trim to 4 cols)
    if not base_atoms_indices:
        if g.shape[1] == 3:
            totals = np.linalg.norm(g[:, :3], axis=1, keepdims=True)
            g = np.hstack([g[:, :3], totals])
        elif g.shape[1] > 4:
            g = g[:, :4]
        return pd.DataFrame(g, columns=['dipole_x','dipole_y','dipole_z','total'])

    # -------- build basis from base_atoms_indices only --------
    origin_set, y_idx, plane_idx = _parse_base(base_atoms_indices)
    origin_pt = coords[origin_set].mean(axis=0)         # centroid of origin set
    basis = _build_basis(coords, origin_pt, y_idx, plane_idx)  # rows: X,Y,Z

    # -------- rotate dipole(s) --------
    d_xyz = g[:, :3]                 # (M,3)
    d_rot = (basis @ d_xyz.T).T      # (M,3)

    if g.shape[1] >= 4:
        total = g[:, 3:4]
    else:
        total = np.linalg.norm(d_xyz, axis=1, keepdims=True)  # invariant anyway

    out = np.hstack([d_rot, total])
    return pd.DataFrame(out, columns=['dipole_x','dipole_y','dipole_z','total'])
