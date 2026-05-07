from .constants import FileExtensions, LinuxCommands
from .file_handlers import *

def convert_smile_file_to_xyz_obabel(smile_filepath, xyz_filepath):
    os.system(LinuxCommands.OBABEL_PREFIX.value+smile_filepath+LinuxCommands.OBABEL_XYZ_SETTINGS_1.value+xyz_filepath+LinuxCommands.OBABEL_XYZ_SETTINGS_2.value)

def convert_smarts_to_xyz_obabel(smarts):
    smarts=''
    return smarts

def convert_smile_str_to_xyz_obabel(smile_str, molecule_name, working_dir=None):
    adjust_working_dir(working_dir)
    xyz_filepath=molecule_name+FileExtensions.XYZ.value
    smile_filepath=convert_smile_str_to_file(smile_str)
    convert_smile_file_to_xyz_obabel(smile_filepath, xyz_filepath)
    os.remove(smile_filepath)
    return xyz_filepath

def convert_csv_to_xyz_obabel(csv_filepath):
    xyz_filepath=csv_filepath
    return xyz_filepath