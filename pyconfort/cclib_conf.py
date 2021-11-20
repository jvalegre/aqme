#####################################################.
#            This file stores all the functions         #
#      used for generating all paramets from cclib       #
#####################################################.

import subprocess, os

from pyconfort.utils import move_file_from_folder
import json
import numpy as np

def calculate_boltz_for_cclib(log_files,name,w_dir_initial,lot,bs):
    """
    Runs goodvibes externally to calculate the --boltz and renames the output 
    file accordingly.

    Parameters
    ----------
    log_files : [type]
        [description]
    name : [type]
        [description]
    w_dir_initial : [type]
        [description]
    lot : [type]
        [description]
    bs : [type]
        [description]
    """

    # GoodVibes must be installed as a module (through pip or conda)
    cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
    for file in log_files:
        cmd_boltz.append(file)
    subprocess.call(cmd_boltz)

    #writing to correct places
    if '/' in str(bs):
        end = str(bs).split('/')[0]
    else:
        end = str(bs)
    destination = f'{w_dir_initial}/QPRED/cclib-json/boltz/{lot}-{end}'

    move_file_from_folder(destination, os.getcwd(),'Goodvibes_{name}.dat')

def log_json(log_file,w_dir_initial,lot,bs):
    if '/' in str(bs):
        end = str(bs).split('/')[0]
    else:
        end = str(bs)
    folder = f'{w_dir_initial}/QPRED/cclib-json/all_confs_cclib/{lot}-{end}'
    try:
        os.makedirs(folder)
    except OSError:
        if os.path.isdir(folder):
            pass

    file_out = f"{folder}/{log_file.split('.log')[0]}.json"
    cmd = ['ccframe',log_file,'-O',file_out]
    cmd_sub = " ".join(cmd)
    os.system(cmd_sub)

def calculate_cclib(log_files,w_dir_initial,lot,bs):
    for file in log_files:
        log_json(file,w_dir_initial,lot,bs)

def calcualte_average_cclib_parameter(json_files,name,w_dir_initial,lot,bs):
    charge,dipole = [],[]
    for file in json_files:
        with open(file,'r') as molfile:
            mol_data = json.load(molfile)

        # get the Boltzmann probabilities of each conformer from a GoodVibes file
        if '/' in str(bs):
            end = str(bs).split('/')[0]
        else:
            end = str(bs)
        boltz_file = f'{w_dir_initial}/QPRED/cclib-json/boltz/{lot}-{end}/Goodvibes_{name}.dat'
        with open(boltz_file,'r') as F:
            boltz_outlines = F.readlines()
        
        token = file.rsplit('.',1)[0]
        for line in boltz_outlines:
            # I remove the NMR from the file names using [0:-4]
            if token in line:
                boltz_factor = float(line.split()[-1])
                break
        #charge
        charge.append(np.array(mol_data['atomcharges']['0']['mulliken']).astype(float)*boltz_factor)
        #dipole
        dipole.append(np.linalg.norm(np.array(mol_data['moments']['0'][1]))*boltz_factor)

    charge = np.sum(charge, axis=0).tolist()
    dipole = np.sum(dipole)

    dict_param = {'name': name,'atomcharges': charge, 'dipole': dipole}

    if '/' in str(bs):
        end = str(bs).split('/')[0]
    else:
        end = str(bs)
    folder = f'{w_dir_initial}/QPRED/cclib-json/average_cclib/{lot}-{end}'

    try:
        os.makedirs(folder)
    except OSError:
        if os.path.isdir(folder):
            pass

    with open(f'{folder}/{name}.json', 'w') as outfile:
        json.dump(dict_param, outfile)
