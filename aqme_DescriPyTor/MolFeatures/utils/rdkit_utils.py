from __future__ import annotations
import numpy as np
import pandas as pd
from collections import OrderedDict
from typing import List, Dict, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdFingerprintGenerator as rdFG  # NEW

# =================== core helper ===================
def _build_mol_from_xyz_df(xyz_df: pd.DataFrame,
                           mol_id: str,
                           atom_col="atom",
                           x_col="x", y_col="y", z_col="z",
                           charge: int = 0) -> Chem.Mol:
    """Convert a single XYZ dataframe into an RDKit molecule."""
    symbols = xyz_df[atom_col].tolist()
    coords = xyz_df[[x_col, y_col, z_col]].values

    mol = Chem.RWMol()
    for s in symbols:
        mol.AddAtom(Chem.Atom(s))

    conf = Chem.Conformer(len(symbols))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    mol.AddConformer(conf)

    # Bond perception (try rdDetermineBonds first)
    try:
        from rdkit.Chem import rdDetermineBonds
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except Exception:
        pass

    mol = mol.GetMol()
    mol.SetProp("_Name", mol_id)
    Chem.SanitizeMol(mol, catchErrors=True)
    return mol


# =================== main converters ===================
def xyz_list_to_mols(xyz_list: List[pd.DataFrame],
                     mol_names: Optional[List[str]] = None) -> Dict[str, Chem.Mol]:
    """Take a list of XYZ dataframes and return {name: Mol}."""
    mols = OrderedDict()
    for i, xyz_df in enumerate(xyz_list):
        name = mol_names[i] if mol_names else f"Mol_{i+1}"
        mols[name] = _build_mol_from_xyz_df(xyz_df, name)
    return mols


def visualize_mols(mols: Dict[str, Chem.Mol],
                   n_cols: int = 4,
                   size: Tuple[int, int] = (250, 250)):
    """Visualize multiple molecules as a grid (returns PIL Image)."""
    mols_list, legends = [], []
    for name, mol in mols.items():
        mol_copy = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol_copy)
        mols_list.append(mol_copy)
        legends.append(name)
    return MolsToGridImage(mols_list, legends=legends, molsPerRow=n_cols, subImgSize=size)


# Optional: scikit-learn for PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def mols_to_fingerprint_df(
    mols: Dict[str, Chem.Mol],
    fp_type: str = "morgan",        # "morgan" | "rdkit" | "maccs"
    nBits: int = 2048,
    radius: int = 2,
    use_chirality: bool = True,
    use_features: bool = False,
    *,
    # --- PCA options ---
    do_pca: bool = False,
    n_components: Optional[int] = 50,   # None -> keep all components
    scale_before_pca: bool = True,      # standardize bit columns before PCA
    random_state: int = 0
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.Series]]:
    """
    Compute fingerprint bit vectors for all molecules.
    If do_pca=True, also returns (pca_df, explained_variance_ratio).

    Returns
    -------
    fp_df : DataFrame (n_mols x n_bits)
    pca_df : DataFrame | None
    pca_var_ratio : Series | None     (explained variance ratio per component)
    """
    data, names = [], []

    # Prepare FP function once
    fp_type_l = fp_type.lower()
    if fp_type_l == "morgan":
        morgan_gen = rdFG.GetMorganGenerator(
            radius=radius,
            fpSize=nBits,
            useChirality=use_chirality,
            useFeatures=use_features
        )
        fp_getter = lambda m: morgan_gen.GetFingerprint(m)  # ExplicitBitVect
    elif fp_type_l == "rdkit":
        fp_getter = lambda m: Chem.RDKFingerprint(m, fpSize=nBits)
    elif fp_type_l == "maccs":
        from rdkit.Chem import MACCSkeys
        fp_getter = lambda m: MACCSkeys.GenMACCSKeys(m)  # fixed 167 bits
    else:
        raise ValueError(f"Unknown fp_type: {fp_type}")

    # Build bit matrix
    for name, mol in mols.items():
        if mol is None:
            continue
        fp = fp_getter(mol)
        arr = np.zeros((fp.GetNumBits(),), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        data.append(arr)
        names.append(name)

    fp_df = pd.DataFrame(data, index=names)
    fp_df.columns = [f"bit_{i}" for i in range(fp_df.shape[1])]

    # Early exit if no PCA requested
    if not do_pca:
        return fp_df, None, None

    # --- PCA pipeline ---
    X = fp_df.values.astype(float)

    # (Optional) standardize; helps PCA treat bits on comparable scale
    if scale_before_pca:
        scaler = StandardScaler(with_mean=True, with_std=True)
        X = scaler.fit_transform(X)

    # n_components=None => all components; if int, min with rank
    pca = PCA(n_components=n_components, random_state=random_state)
    X_pca = pca.fit_transform(X)

    # Build PCA dataframe with clear column names
    n_comp = X_pca.shape[1]
    pca_cols = [f"PC{i+1}" for i in range(n_comp)]
    pca_df = pd.DataFrame(X_pca, index=fp_df.index, columns=pca_cols)

    var_ratio = pd.Series(pca.explained_variance_ratio_, index=pca_cols, name="explained_variance_ratio")

    return fp_df, pca_df, var_ratio
# =================== usage example ===================
# xyz_list = [df1, df2, df3]  # list of pandas DataFrames (each with symbol, x, y, z)
# mols = xyz_list_to_mols(xyz_list, mol_names=["mol1", "mol2", "mol3"])
# img = visualize_mols(mols)
# display(img)
# fps = mols_to_fingerprint_df(mols, fp_type="morgan", nBits=1024)
# print(fps.head())
