import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#lib imports
try:

    from ..utils.constants import LinuxCommands
    from ..utils.file_handlers import *

except ImportError:
    from utils.constants import LinuxCommands
    from utils.file_handlers import *

from typing import List, Optional

def run_gaussian_calculation(xyz_filename: str) -> None:
    """
    Placeholder for running a Gaussian calculation.

    Args:
        xyz_filename (str): The filename of the XYZ file.
    """
    pass

def run_xtb_calcuation(xyz_filename: str, detailed_input: Optional[str] = None) -> None:
    """
    Runs an xTB calculation on a given XYZ file with optional detailed input.

    Args:
        xyz_filename (str): The filename of the XYZ file.
        detailed_input (str, optional): Detailed input string for the xTB calculation. Defaults to None.
    """
    command = LinuxCommands.XTB_INPUT_PREFIX.value + detailed_input + ' ' + xyz_filename + LinuxCommands.XTB_SUFIX.value if detailed_input else LinuxCommands.XTB_PREFIX.value + xyz_filename + LinuxCommands.XTB_SUFIX.value
    os.system(command)

def run_crest_calculation(xyz_filename: str, detailed_input: Optional[str] = None) -> None:
    """
    Runs a CREST calculation on a given XYZ file with optional detailed input.

    Args:
        xyz_filename (str): The filename of the XYZ file.
        detailed_input (str, optional): Detailed input string for the CREST calculation. Defaults to None.
    """
    command = LinuxCommands.CREST_INPUT_PREFIX.value + detailed_input + ' ' + xyz_filename + LinuxCommands.CREST_SUFIX.value if detailed_input else LinuxCommands.CREST_PREFIX.value + xyz_filename + LinuxCommands.CREST_SUFIX.value
    os.system(command)

def submit_to_gaussian_calculation(path=None) -> None:
    """
    Submits files to Gaussian calculation. Assumes .com files are located in the 'com' directory.
    """
    os.chdir(path) if path else os.chdir(os.getcwd())
    os.system(LinuxCommands.COPY.value + ' ' + LinuxCommands.GAUSS_SUFIX.value + ' .')
    for file in os.listdir():
        if file.endswith('.com'):
            output_name = file.split('.')[0] + '.log'
            os.system(f'./g16 {file} {output_name} -q m')
    os.remove('g16')
    os.chdir('..')

def submit_to_xtb_calculation(files_to_send: List[str]) -> str:
    """
    Submits files to xTB calculation.

    Args:
        files_to_send (List[str]): List of file names to be sent for xTB calculation.

    Returns:
        str: Path to the output folder.
    """
    output_folder_filepath = ''
    # Implementation for xTB calculation submission
    return output_folder_filepath

def submit_to_crest_calculation(files_to_send: List[str]) -> str:
    """
    Submits files to CREST calculation.

    Args:
        files_to_send (List[str]): List of file names to be sent for CREST calculation.

    Returns:
        str: Path to the output folder.
    """
    output_folder_filepath = ''
    # Implementation for CREST calculation submission
    return output_folder_filepath

def submit_to_calculation(calculation_method: str = 'gaussian') -> Optional[str]:
    """
    Submits files to a specified calculation method.

    Args:
        calculation_method (str): The calculation method to use. Defaults to 'gaussian'.

    Returns:
        Optional[str]: Path to the output folder, if applicable.
    """
    if calculation_method == 'gaussian':
        return submit_to_gaussian_calculation()
    elif calculation_method == 'xtb':
        return submit_to_xtb_calculation([])
    elif calculation_method == 'crest':
        return submit_to_crest_calculation([])
    return None
