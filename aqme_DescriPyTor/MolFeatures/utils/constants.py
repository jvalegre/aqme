from enum import Enum

class FileExtensions(Enum):
    """
    Hold commonly used file extensions
    """
    SMI='.smi'
    XYZ='.xyz'
    CSV='.csv'
    ZIP='.zip'
    PPT='.ppt'
    CIF='.cif'
    MOL='.mol'
    PDB='.pdb'

class LinuxCommands(Enum):
    COPY='cp'
    OBABEL_PREFIX='/gpfs0/gaus/users/itamarwa/transfer_to_itamar/build/bin/obabel '
    OBABEL_XYZ_SETTINGS_1=' -O '
    OBABEL_XYZ_SETTINGS_2=' --gen3d'
    XTB_INPUT_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb --input ' # xtb_di.inp 
    XTB_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb '
    XTB_SUFIX=' --ohess --dipole --pop'
    CREST_INPUT_PREFIX='/gpfs0/gaus/projects/crest --input '
    CREST_PREFIX='/gpfs0/gaus/projects/crest '	
    GAUSS_SUFIX='/gpfs0/gaus/users/kozuch/home/scripts/gaussian/g16'    
    #GAUSS_SUFIX='/gpfs0/gaus/users/edenspec/g16'

class GaussianDict(Enum):
    gaussian_dict = {
        "basis_sets": {
            "STO-3G": "3-21G",
            "6-31G": "6-31G(d)",
            "6-31G*": "6-31G(d,p)",
            "6-311G*": "6-311G(d,p)",
            "cc-pVDZ": "cc-pVDZ",
            "cc-pVTZ": "cc-pVTZ",
            "aug-cc-pVDZ": "aug-cc-pVDZ",
            "aug-cc-pVTZ": "aug-cc-pVTZ",
            "def2-SVP": "def2-SVP",
            "def2-TZVP": "def2-TZVP",
            "def2-QZVP": "def2-QZVP"
        },
        "functionals": {
            "B3LYP": "B3LYP",
            "PBE0": "PBE0",
            "M062X": "M062X",
            "M06-2X": "M06-2X",
            "MP2": "MP2",
            "CCSD": "CCSD",
            "CCSD(T)": "CCSD(T)"
        },
        "opt": {
            "method": "default",
            "max_cycles": 200,
            "max_step_size": 0.3,
            "rms_force": 0.0001,
            "rms_displacement": 0.001,
            "tight": False,
            "calc_hessian": True,
            "freq": True
        }

    }