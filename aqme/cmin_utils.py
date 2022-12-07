#####################################################.
#             This file stores functions            #
#             used in in the CMIN module            #
#####################################################.

import sys
from rdkit.Chem import AllChem as Chem
import pandas as pd

hartree_to_kcal = 627.509


def rdkit_sdf_read(file, args):
    """
    Reads sdf files and stops the execution if the file was not accessible.
    rdkit.Chem.Mol objects
    """
    
    inmols = Chem.SDMolSupplier(file, removeHs=False, sanitize=False)

    if inmols is None:
        args.log.write(f"Could not open {file}")
        args.log.finalize()
        sys.exit()
    return inmols


def creation_of_dup_csv_cmin(cmin):
    """
    Generates a pandas.DataFrame object with the appropiate columns for the
    conformational search and the minimization.

    Parameters
    ----------
    cmin : str
        Minimization method. Current valid methods are: ['xtb','ani']

    Returns
    -------
    pandas.DataFrame
    """

    # Boolean aliases from args
    is_xtb = cmin == "xtb"
    is_ani = cmin == "ani"

    # column blocks definitions

    xtb_columns = [
        "xTB-Initial-samples",
        "xTB-energy-window",
        "xTB-initial_energy_threshold",
        "xTB-RMSD-and-energy-duplicates",
        "xTB-Unique-conformers",
    ]
    ANI_columns = [
        "ANI-Initial-samples",
        "ANI-energy-window",
        "ANI-initial_energy_threshold",
        "ANI-RMSD-and-energy-duplicates",
        "ANI-Unique-conformers",
    ]
    end_columns = ["CMIN time (seconds)", "Overall charge"]

    # Check Minimization Method
    if is_ani:
        columns = ANI_columns
    if is_xtb:  # is_ani and is_xtb will not happen, but this is what was written
        columns = xtb_columns

    columns += end_columns
    return pd.DataFrame(columns=columns)
