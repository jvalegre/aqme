import os, glob
from pathlib import Path
from aqme.csearch import csearch
from aqme.qdescp import qdescp

sdf_rdkit_files = glob.glob(f'CSEARCH/rdkit/*.sdf')

qdescp(files=sdf_rdkit_files, boltz=True, program='xtb')