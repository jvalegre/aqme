# Feature Extraction: `Molecules` and `Molecule`

This document explains the feature‑extraction layer used to turn quantum‑chemistry outputs and geometries into model‑ready descriptors. It focuses on the **`Molecules`** batch orchestration class and uses **`Molecule`** methods where needed.

---

## Overview

* **Input**: a directory of per‑molecule **`.feather`** files (produced by your Gaussian/processing pipeline).
* **Core classes**:

  * **`Molecules`** — loads many molecules, provides batch extractors that return wide, analysis‑ready DataFrames.
  * **`Molecule`** — wraps one molecule’s geometry, connectivity, charges, dipoles, vibrations, etc., and exposes descriptor calculators.
* **Output**: horizontally concatenated DataFrames (per descriptor family) keyed by molecule name, plus optional visualizations and per‑molecule exports.

> **Indexing note:** Most public methods accept human‑friendly **1‑based** atom indices; internally, indices are normalized (via `adjust_indices`) to match underlying arrays.

---

## `Molecules`: batch orchestrator

### Construction

```python
from data_extractor import Molecules

molset = Molecules("/path/to/feather_dir", threshold=1.82)
```

* Scans the directory for `*.feather`, creates a `Molecule` for each.
* Tracks which files loaded successfully (`success_molecules`) and which failed (`failed_molecules`).
* `threshold` controls the **distance cutoff** used to infer connectivity (bonds) for all created molecules.

### Selection & export helpers

```python
molset.filter_molecules([0, 3, 5])      # keep a subset by position
molset.export_all_xyz()                 # write per‑molecule XYZ files to ./xyz_files
molset.extract_all_dfs()                # write each molecule’s internal DataFrames as CSVs
molset.visualize_molecules([0, 1])      # quick 3D previews
molset.visualize_smallest_molecule()    # pick smallest by atom count and visualize
```

### Batch descriptor extractors (wide tables)

Each method iterates over all loaded molecules, calls the respective `Molecule` method, and returns a **horizontally concatenated** DataFrame (one block per molecule). Typical usage:

```python
# Sterimol on atom pairs [[1,6], [3,4]] using CPK radii
steri = molset.get_sterimol_dict([[1,6],[3,4]], radii='CPK')

# NPA on base trio [1,2,3] (optionally also sub‑atoms)
npa = molset.get_npa_dict([1,2,3], sub_atoms=[5,6,7])

# Ring vibrations on para/cross pattern seeds
rings = molset.get_ring_vibration_dict([[8],[14]]])

# Dipoles for base specs (single or multiple groups)
dips = molset.get_dipole_dict([[1,2,3], 5, 6])  # origin set [1,2,3], y=5, plane=6

# Bond metrics
angles = molset.get_bond_angle_dict([1,2,3])
lengths = molset.get_bond_length_dict([[1,2],[3,4]])

# Vibrational descriptors
stretch = molset.get_stretch_vibration_dict([[1,2],[3,4]], threshold=3000)
bend    = molset.get_bend_vibration_dict([[1,2],[3,4]],  threshold=1300)

# Charges & differences (supports multiple charge types)
charges     = molset.get_charge_df_dict([3,5,7,9])          # type='all' inside
charge_diff = molset.get_charge_diff_df_dict([[1,2],[3,4]]) # type='all' by default
```

**Return shape:** A wide table where each molecule contributes a column block. When a molecule’s computation fails, it is skipped (and a message is logged/printed).

### Turn GUI inputs into a feature set: `get_molecules_features_set(...)`

This convenience routine consumes text‑box‑like inputs (strings or widgets), parses them into lists of atom indices, runs the requested extractors, and returns a single **feature matrix** (`res_df`).

```python
res = molset.get_molecules_features_set(entry_widgets={
    'Ring':        '[[8],[14]]',
    'Stretching':  '[[1,2],[3,4]]',
    'Bending':     '[[1,2],[3,4]]',
    'NPA':         '[1,2,3]',
    'Dipole':      '[[1,2,3], 5, 6]',
    'Charges':     '[3,5,7,9]',
    'Charge-Diff': '[[1,2],[3,4]]',
    'Sterimol':    '[[1,6],[3,4]]',
    'Bond-Angle':  '[1,2,3]',
    'Bond-Length': '[[1,2],[3,4]]',
}, parameters={'Radii': 'CPK', 'Isotropic': True}, save_as=True,
   csv_file_name='features_csv')
```

* Parses strings like `'[[1,6],[3,4]]'` into nested lists automatically.
* Runs only the steps for which you provided inputs (ring, stretching, bending, NPA, dipole, charges, charge\_diff, sterimol, bond\_angle, bond\_length).
* When `Isotropic=True`, appends **polarizability** and **energy** per molecule.
* Produces an interactive correlation heatmap and a CSV with highly correlated pairs.
* If `save_as=True`, writes `features_csv.csv` and `features_csv_correlation_table.csv`.

---

## `Molecule`: per‑molecule descriptors

### Initialization

```python
mol = Molecule("/path/to/LS1716_optimized.feather", threshold=1.82)
```

* Loads the `.feather` through your `feather_file_handler` (no need to change CWD manually; class handles it).
* Populates:

  * **Geometry**: `xyz_df` (atom, x, y, z), `coordinates_array`.
  * **Connectivity**: `bonds_df` (inferred via distance cutoff), `atype_list` (atom types).
  * **Quantum outputs**: `gauss_dipole_df`, `polarizability_df`, `info_df`, `energy_value`, `charge_dict`, `vibration_dict` (plus `vibration_mode_dict`).
* The `threshold` controls bond detection; increase slightly for longer bonds, decrease for stricter graphs.

### Core utilities

* **Renumber atoms**: `renumber_atoms({old: new, ...})`
  Rebuilds `xyz_df`, connectivity, types, and reindexes `charge_dict`/`vibration_dict` to match the new order. A helper `reindex_vibration_dict()` is available when you have a `new_index_order` prepared.

* **Export & visualize**:
  `write_xyz_file()` → write per‑molecule XYZ;
  `write_csv_files()` → dump all internal DataFrames;
  `visualize_molecule(vector=None)` → 3D viewer (optionally override dipole vector).

* **Descriptor summary**:
  `get_molecule_descriptors_df()` → concatenates energy, Gaussian dipole, and Sterimol into one small table and sorts by `energy_value`.

### Sterimol

```python
mol.get_sterimol(base_atoms=None, radii='CPK', sub_structure=True, drop_atoms=None, mode='all')
```

* If `base_atoms=None`, an automatic base is detected.
* Accepts either a **single pair** `[i,j]` or a **list of pairs** `[[i,j], [k,l], ...]` and concatenates results.
* `process_sterimol_atom_group(...)` is the worker that computes Sterimol for one pair.

### Charges (NBO/Hirshfeld/CM5…)

```python
mol.get_charge_df([3,5,7,9], type='all')             # dict of DataFrames per charge type
mol.get_charge_df([3,5,7,9], type=['nbo','cm5'])     # select types

# Differences per pair (single or many)
mol.get_charge_diff_df([[1,2],[3,4]], type='all')
```

* `type='all'` returns a dictionary of DataFrames for every available type in `charge_dict` (silently skips missing types).
* Differences compute `atom_i − atom_j` per mode and label the columns as `diff_i-j`.

### Dipoles (Gaussian‑based)

```python
# Atoms can be one of:
#   [o, y, plane]
#   [o1, o2, ..., y, plane]
#   [[o1, o2, ...], y, plane]

df = mol.get_dipole_gaussian_df([ [1,2,3], 5, 6 ], visualize_bool=True)
```

* Flexible “R‑style” base specification: single origin, multiple origins, or an explicit list of origins.
* Returns a one‑row DataFrame per spec; `get_dipole_gaussian_df(...)` handles a **list of specs** and concatenates rows.
* Optional visualization uses the centroid of the origin set as the vector origin.

### Bond metrics

```python
ang = mol.get_bond_angle([1,2,3])
len1 = mol.get_bond_length([1,2])
lenN = mol.get_bond_length([[1,2],[3,4]])
```

### Vibrational descriptors

**Stretch (pair‑aligned; typical high‑frequency cut)**

```python
vst = mol.get_stretch_vibration([1,2], threshold=3000)
# or many: mol.get_stretch_vibration([[1,2],[3,4]], threshold=3000)
```

* Validates that the pair is bonded.
* Computes dot‑product‑based alignment per mode, selects the maximum within the threshold.

**Bend (pair with common center)**

```python
vbd = mol.get_bend_vibration([1,2], threshold=1300)
```

* Builds a local adjacency, finds a shared center, then uses cross‑product magnitudes to identify the strongest bending mode.

**Ring vibrations (para / cross analysis)**

```python
rings = mol.get_ring_vibrations([[1,4],[2,5]], return_nan_on_empty=True)
```

* Resolves a benzene‑like pattern around the provided seed(s).
* Returns `cross`, `cross_angle`, `para`, `para_angle`. If filters eliminate everything and `return_nan_on_empty=True`, a NaN row is returned (instead of `None`).

### Geometry transforms and helpers

* `get_coordinates_mean_point([i,j,k])` → centroid of selected atoms.
* `get_coordination_transformation_df([i,j,k], origin=None)` → axis transform (e.g., for frame‑aligned analyses).
* `swap_atom_pair((i,j))` → swap coordinates of two atoms (useful for quick what‑ifs).

---

## Practical tips

* **Indexing**: You can pass 1‑based indices; the code standardizes them internally. When you rename or reorder atoms, re‑extract descriptors to ensure consistency.
* **Connectivity cutoff**: If bonds are missing/extra, adjust `threshold` in the constructor (default 1.82 Å).
* **Robustness**: Most batch methods skip failing molecules and continue, printing a helpful message.
* **Correlations**: Use the GUI aggregator (`get_molecules_features_set`) to quickly assemble a matrix and auto‑inspect highly correlated features before modeling.

---

## Minimal end‑to‑end example

```python
molset = Molecules("./mols", threshold=1.82)

# Extract a compact, model‑ready feature set
feat = molset.get_molecules_features_set(
    entry_widgets={
        'Sterimol': '[[1,6],[3,4]]',
        'Stretching': '[[1,2],[3,4]]',
        'Bond-Angle': '[1,2,3]',
        'Dipole': '[[1,2,3], 5, 6]'
    },
    parameters={'Radii': 'CPK', 'Isotropic': True},
    save_as=True,
    csv_file_name='molecule_features'
)

print(feat.shape)
```

