import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from ..utils.constants import *
    from ..utils.file_handlers import *
    from ..M1_pre_calculations.outwards_senders import *
except ImportError:
    from utils.constants import *
    from utils.file_handlers import *
    from M1_pre_calculations.outwards_senders import *

from typing import List, Tuple

class Module1Handler:
    """
    A handler class for managing molecular calculations using Gaussian input files.

    Attributes:
        names_list (List[str]): List of molecule names extracted from the CSV file.
        identifier_list (List[str]): List of molecular identifiers (e.g., SMILES strings).
        working_dir (str): Directory path for working with files.
        smiles_list (List[str]): List of SMILES strings, representing molecules.
        inital_xyz_list (List[str]): List of file paths or data for initial XYZ coordinates.
        calculation_method (str): Method of calculation, e.g., 'gaussian'.
    """

    def __init__(self, id_csv_filepath: str, working_dir: str = None, calculation_method: str = 'gaussian'):
        """
        Initializes the Module1Handler with necessary parameters.

        Args:
            id_csv_filepath (str): File path to the CSV file containing molecular identifiers.
            working_dir (str, optional): Working directory for file operations. Defaults to the current directory.
            calculation_method (str): The method of calculation to be used. Defaults to 'gaussian'.
        """
                  

        
        self.names_list, self.identifier_list = process_identifiers_csv_file(id_csv_filepath)
        self.working_dir =os.path.dirname(id_csv_filepath) if working_dir is None else working_dir
        self.smiles_list = self.identifier_list  # Update this if the identifier isn't SMILES
        self.inital_xyz_list = smiles_to_coordinates(self.identifier_list)
        self.calculation_method = calculation_method

    def make_gaussian_input(self) -> None:
        """
        Generates Gaussian input files from the initial XYZ list.
        """
        xyz_files_list_to_gaussian_files(self.inital_xyz_list)

    def submit_to_calculation(self) -> None:
        """
        Submits the generated Gaussian input for calculation.
        """
        submit_to_calculation(self.calculation_method)

    def execute_full_operation(self) -> None:
        """
        Executes the full operation from generating Gaussian input to submitting for calculation.
        """
        self.make_gaussian_input()
        self.submit_to_calculation()

if __name__ == '__main__':
    csv_filepath = r'smiles_list.csv'  # Replace with get_csv_filepath() if necessary
    calculation_method = 'gaussian'
    module_handler = Module1Handler(csv_filepath, calculation_method)
    module_handler.execute_full_operation()
