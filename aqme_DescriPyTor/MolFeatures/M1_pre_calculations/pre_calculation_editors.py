import os
import sys
import numpy as np
from ..utils.file_handlers import *
from ..utils.help_functions import *

def prepare_xyz_file_for_xtb(xyz_filepath):
    os.system('echo >> '+xyz_filepath)
    os.system('echo >> '+xyz_filepath)
    os.system('echo "\$write" >> '+xyz_filepath)
    os.system('echo "   output file=properties.out" >> '+xyz_filepath)

def preprocess_data_for_xtb(xyz_filepaths):
    for xyz_filepath in xyz_filepaths:
        prepare_xyz_file_for_xtb(xyz_filepath)
    processed_data=xyz_filepaths
    return processed_data


def prepare_xyz_files_for_gaussian(xyz_filepath, functional='HF', basisset='6-31G(d)', charge='0 1', nbo_answer='n',title='title',task='sp'):
    """
    methods->B3LYP, HF
    basis sets->STO-3G, 3-21G, 6-31G, 6-31G(d), 6-31G(d,p), 6-31+G(d,p), 6-31+G(d,p),
     6-311G(d,p), 6-311+G(d,p), 6-311++G(d,p), 6-311++G(2d,p), 6-311++G(3df,2p), https://gaussian.com/basissets/
    """
    os.chdir(xyz_filepath)
    xyz_list=[filename for filename in os.listdir() if filename.endswith('.xyz')]
    names_list=[filename.split('.')[0] for filename in xyz_list]
    for file, name in zip(xyz_list, names_list):
        with open(f'{name}.com', 'w') as my_file:
            my_file.write(f"%mem=100GB\n%nproc=32\n%chk={name}.chk\n")
            if task=='sp':
                my_file.write(f"#P {functional}/{basisset}\n\n{title}\n\n{charge}\n")
            elif task=='opt':
                my_file.write(f"#P {functional}/{basisset} Opt Freq pop=nbo\n\n{title}\n\n{charge}\n")
            xyz_df=get_df_from_file(file)
            atoms_np_array = np.array(xyz_df)
            for atom_np_array in atoms_np_array:
                try:
                    my_file.write("{:1} {:11.5} {:11.5} {:11.5}\n".format(*atom_np_array))
                except:
                    print('error')
            if nbo_answer == 'y':
                my_file.write("# pop=(full,nbo)\n")
            my_file.write("\n")
    os.mkdir('com')
    for file in os.listdir('.'):
        if file.endswith('.com'):
            os.rename(file, f"com/{file}")    

def preprocess_data_for_gaussian(xyz_filepath):
    ## fine a use for this function
    prepare_xyz_files_for_gaussian(xyz_filepath)
    return 

def prepare_xyz_file_for_crest(xyz_filepath):
    pass

def preprocess_data_for_crest(xyz_filepaths):
    processed_data=''
    return processed_data

def preprocess_data_for_calculation(xyz_filepaths, calculation_method='gaussian'):
    if calculation_method=='gaussian':
        processed_data=preprocess_data_for_gaussian(xyz_filepaths)
    if calculation_method=='xtb':
        processed_data=preprocess_data_for_xtb(xyz_filepaths)
    if calculation_method=='crest':
        processed_data=preprocess_data_for_crest(xyz_filepaths)
    return processed_data

if __name__=='__main__':
    os.chdir(r'C:\Users\edens\Documents\GitHub\Automation_code-main\XYZ files for runs\Easy')
    preprocess_data_for_gaussian(r'C:\Users\edens\Documents\GitHub\Automation_code-main\XYZ files for runs\Easy')
