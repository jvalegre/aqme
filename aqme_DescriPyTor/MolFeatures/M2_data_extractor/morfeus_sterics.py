from __future__ import annotations
from typing import Iterable, Sequence, Optional, Dict, Any, Union
import numpy as np


def wrapper(mols):

    results_list=[]
    for mol in mols.molecules:
        res=morfeus_sterics(
            elements=mol.xyz_df['atom'].tolist(),
            coordinates=mol.xyz_df[['x','y','z']].values,
            metal_index=12,
            cone_center_index=12,
            excluded_atoms=None,
            include_hs=True,
            radius_first=3.5,
            radius_second=5.5,
            radii_type="crc",
            radii_scale=1.17,
            density=0.001,
            z_axis_atoms=[12,5],
            xz_plane_atoms=[11,4],
            cone_method="libconeangle",
        )
        results_list.append(res)
    return results_list

def morfeus_sterics(
    *,
    # geometry input (provide either `geometry_file` OR `elements`+`coordinates`)
    geometry_file: Optional[str] = None,
    elements: Optional[Iterable[Union[str, int]]] = None,
    coordinates: Optional[Iterable[Iterable[float]]] = None,

    # indices are 1-INDEXED per morfeus API
    metal_index: int = 1,                  # center for buried volume
    cone_center_index: Optional[int] = None,  # center for cone angle (e.g., donor atom). Defaults to `metal_index`.

    # buried volume options
    excluded_atoms: Optional[Sequence[int]] = None,  # atoms (1-indexed) to exclude besides metal
    include_hs: bool = False,
    radius_first: float = 3.5,             # “first-sphere” (Cavallo/SambVca default)
    radius_second: float = 5.5,            # “second-sphere” (common choice; make it 5.0–5.5 as needed)
    radii_type: str = "bondi",             # 'alvarez' | 'bondi' | 'crc' | 'truhlar'
    radii_scale: float = 1.17,
    density: float = 0.001,                # sphere point density (Å^3 per point)

    # orientation (only needed if you want quadrants/octants & steric maps)
    z_axis_atoms: Optional[Sequence[int]] = None,    # 1-indexed
    xz_plane_atoms: Optional[Sequence[int]] = None,  # 1-indexed

    # cone angle calc
    cone_method: str = "libconeangle",     # 'libconeangle' (fast) or 'internal'
) -> Dict[str, Any]:
    """
    Returns a dict with:
      - cone_angle_deg, cone_tangent_atoms (1-indexed)
      - vbur_first_{fraction,percent,volume_A3,free_volume_A3}
      - vbur_second_{fraction,percent,volume_A3,free_volume_A3}
      - distal_volume_first_A3 (if computed), molecular_volume_first_A3
      - (optional) quadrants_first, octants_first if orientation provided
    """
    try:
        from morfeus import read_xyz, BuriedVolume, ConeAngle
    except Exception as e:
        raise ImportError("morfeus must be installed and importable") from e

    # Load geometry
    if geometry_file is not None:
        elements, coordinates = read_xyz(geometry_file)
    if elements is None or coordinates is None:
        raise ValueError("Provide `geometry_file` OR both `elements` and `coordinates`.")

    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("`coordinates` must be (N, 3).")

    # ---- Cone angle (center is usually the donor atom, *not* the metal) ----
    if cone_center_index is None:
        cone_center_index = metal_index
    cone = ConeAngle(
        elements, coordinates, cone_center_index,
        radii_type="crc",   # CRC is ConeAngle default in morfeus
        method=cone_method
    )
    cone_angle_deg = float(cone.cone_angle)
    cone_tangent_atoms = list(cone.tangent_atoms)  # already 1-indexed

    # ---- Buried volume: first sphere ----
    bv1 = BuriedVolume(
        elements, coordinates, metal_index,
        excluded_atoms=excluded_atoms, include_hs=include_hs,
        radius=radius_first, radii_type=radii_type, radii_scale=radii_scale,
        density=density, z_axis_atoms=z_axis_atoms, xz_plane_atoms=xz_plane_atoms
    )
    # Optional quadrant/octant + distal volume
    if z_axis_atoms is not None:
        bv1.octant_analysis()
    try:
        bv1.compute_distal_volume(method="buried_volume", octants=(z_axis_atoms is not None))
    except Exception:
        # distal volume is optional; ignore if not available
        pass

    # ---- Buried volume: second sphere (larger radius) ----
    bv2 = BuriedVolume(
        elements, coordinates, metal_index,
        excluded_atoms=excluded_atoms, include_hs=include_hs,
        radius=radius_second, radii_type=radii_type, radii_scale=radii_scale,
        density=density, z_axis_atoms=z_axis_atoms, xz_plane_atoms=xz_plane_atoms
    )

    result: Dict[str, Any] = {
        # Cone angle
        "cone_angle_deg": cone_angle_deg,
        "cone_tangent_atoms": cone_tangent_atoms,  # 1-indexed

        # First-sphere buried volume (3.5 Å by default)
        "vbur_first_fraction": float(bv1.fraction_buried_volume),
        "vbur_first_percent": float(bv1.fraction_buried_volume * 100.0),
        "vbur_first_volume_A3": float(bv1.buried_volume),
        "free_volume_first_A3": float(bv1.free_volume),
        "distal_volume_first_A3": getattr(bv1, "distal_volume", None),
        "molecular_volume_first_A3": getattr(bv1, "molecular_volume", None),

        # Second-sphere buried volume (e.g., 5.5 Å)
        "vbur_second_fraction": float(bv2.fraction_buried_volume),
        "vbur_second_percent": float(bv2.fraction_buried_volume * 100.0),
        "vbur_second_volume_A3": float(bv2.buried_volume),
        "free_volume_second_A3": float(bv2.free_volume),
    }

    if hasattr(bv1, "quadrants"):
        result["quadrants_first"] = bv1.quadrants
    if hasattr(bv1, "octants"):
        result["octants_first"] = bv1.octants

    return result


import re
import numpy as np
import pandas as pd
from typing import List, Dict, Any

def morfeus_results_to_df(mols, sort_indices: bool = True) -> pd.DataFrame:
    """
    Build a tidy DataFrame from morfeus_sterics results.
    Index = molecule name (if available), else mol_#.
    """
    rows = {}
    oct_labels = ["+x+y+z","+x+y−z","+x−y+z","+x−y−z","−x+y+z","−x+y−z","−x−y+z","−x−y−z"]
    results_list = wrapper(mols)
    for i, (mol, res) in enumerate(zip(mols.molecules, results_list), start=1):
        name = (
            getattr(mol, "name", None)
            or getattr(mol, "molecule_name", None)
            or getattr(mol, "id", None)
            or f"mol_{i}"
        )

        row = {
            "cone_angle_deg":          res.get("cone_angle_deg", np.nan),
            "%Vbur_first":             res.get("vbur_first_percent", np.nan),
            "%Vbur_second":            res.get("vbur_second_percent", np.nan),
            "Δ%Vbur(second-first)":    (
                res.get("vbur_second_percent", np.nan) - res.get("vbur_first_percent", np.nan)
            ),
            "Vbur1_volume_A3":         res.get("vbur_first_volume_A3", np.nan),
            "Free1_volume_A3":         res.get("free_volume_first_A3", np.nan),
            "Vbur2_volume_A3":         res.get("vbur_second_volume_A3", np.nan),
            "Free2_volume_A3":         res.get("free_volume_second_A3", np.nan),
            "donor_index_used":        res.get("donor_index_used", np.nan),
            "donor_guess_used":        bool(res.get("donor_guess_used", False)),
        }

        # Quadrants (Q1..Q4) if present
        q = res.get("quadrants_first", None)
        if q is not None:
            if isinstance(q, dict):
                for k in ["Q1","Q2","Q3","Q4"]:
                    row[f"{k}_percent"] = float(q.get(k, np.nan))
            else:
                q = list(q)
                for idx, k in enumerate(["Q1","Q2","Q3","Q4"]):
                    row[f"{k}_percent"] = float(q[idx]) if idx < len(q) else np.nan

        # Octants (8 bins) if present
        o = res.get("octants_first", None)
        if o is not None:
            if isinstance(o, dict):
                for lab in oct_labels:
                    row[f"oct_{lab}_percent"] = float(o.get(lab, np.nan))
            else:
                o = list(o)
                for idx, lab in enumerate(oct_labels):
                    row[f"oct_{lab}_percent"] = float(o[idx]) if idx < len(o) else np.nan

        rows[name] = row

    df = pd.DataFrame.from_dict(rows, orient="index")

    # Optional: numeric-aware sort for IDs like "LS2002"
    if sort_indices:
        try:
            df = df.sort_index(
                key=lambda idx: idx.map(lambda x: int(re.search(r"\d+", str(x)).group())
                                        if re.search(r"\d+", str(x)) else float("inf"))
            )
        except Exception:
            pass

    return df
