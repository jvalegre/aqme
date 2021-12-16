#####################################################.
#        This file stores all the functions         #
#    used for genrating all nics input/output       #
#####################################################.
import os
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
from numpy.linalg.linalg import LinAlgError

from pyconfort.utils import periodic_table, move_file_from_folder

def find_coeffplane(ringatoms, CARTESIANS,log):
    rotated = 0
    xyz = np.array([CARTESIANS[i] for i in ringatoms])
    n,_ = xyz.shape
    centroid = xyz.sum(axis=0)/n
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    A = np.array([[np.sum(x**2), np.sum(x*y), np.sum(x)],
                  [np.sum(x*y) , np.sum(y**2),np.sum(y)],
                  [np.sum(x)   , np.sum(y)   ,    n    ]])
    b = np.array([[np.sum(x*z)],
                  [np.sum(y*z)],
                  [np.sum(z)  ]])

    if A[0,2] == 0.0 and A[1,2] == 0.0:
        rotated = 3
        log.write("x  Can't define a ring by points in a line")
    if A[0,2] == 0.0:
        rotated = 1
    if A[1,2] == 0.0:
        rotated = 2

    try: 
        coeff = np.linalg.solve(A, b)
    except LinAlgError: 
        coeff = np.zeros(shape=(3,1))

    return coeff, centroid, rotated

def update_coord(natoms,atomtypes,cartesians,args,log,name,w_dir_initial,type):
    #find the ring atoms in the File
    with open(w_dir_initial+'/'+args.nics_atoms_file,'r') as F:
        filelines = F.readlines()
    ringatoms = []
    for line in filelines:
        split_line = line.strip().split(',')
        if split_line[0] == name.strip().split()[0]:
            for i in range(1,len(split_line)):
                ringatoms.append(int(split_line[i])-1)
            break

    coeffplane, centroid, rotated = find_coeffplane(ringatoms,cartesians,log)

    xcoeff= coeffplane.tolist()[0][0]
    ycoeff= coeffplane.tolist()[1][0]
    cval= coeffplane.tolist()[2][0]

    rawvector=np.array([xcoeff,ycoeff,-1]) #Need to make into unit vector
    vector_norm = rawvector/np.linalg.norm(rawvector)

    if vector_norm[2] < 0: 
        vector_norm = -vector_norm

    if rotated == 1:
        log.write("************ coordinated system was rotated! ***********")
        if vector_norm[2] < 0: 
            vector_norm = -vector_norm
        log.write(f"Unit vector: {vector_norm}")
    if rotated == 2:
        log.write("************ coordinated system was rotated! ***********")
        if vector_norm[2] < 0: 
            vector_norm = -vector_norm
        log.write(f"Unit vector: {vector_norm}")
    if rotated == 3:
        log.write("didn't I tell you this was a bad idea?")
    
    if type=='write':
        spacing = float(args.nics_range)/float(args.nics_number)
        for w in range(-args.nics_number,args.nics_number+1):
            scalefactor = w*spacing
            vector_final = vector_norm*scalefactor + centroid
            x,y,z = vector_final
            natoms += 1
            atomtypes.append("Bq")
            cartesians.append([x,y,z])

        return natoms,atomtypes,cartesians

    if type=='read':
        return vector_norm,centroid

def getSHIELDING(outlines,args,log):
  NATOMS = 0
  stand_or = 0
  ATOMTYPES, CARTESIANS,NMR  = [],[],[]
  #find NATOMS
  for i in range(0,len(outlines)):
      if outlines[i].find("Input orientation") > -1:
          stand_or = i
      if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
          NATOMS = i-stand_or-6
          break

  ATOMTYPES, CARTESIANS,stand_or = get_coords(outlines, stand_or, NATOMS, periodic_table, ATOMTYPES, CARTESIANS)

  NMR = []
  for i in range(0,len(outlines)):
      if outlines[i].find(" SCF GIAO Magnetic shielding tensor (ppm):") > -1:
          j = 0
          while j < NATOMS*5:
              item = {}
              item['atom_index'] = int(outlines[i+1+j].split()[0])
              item['elementID'] = outlines[i+1+j].split()[1]
              item['isotropic'] = float(outlines[i+1+j].split()[4])
              item['anisotropy'] = float(outlines[i+1+j].split()[7])
              item['xx'] = float(outlines[i+2+j].split()[1])
              item['yx'] = float(outlines[i+2+j].split()[3])
              item['zx'] = float(outlines[i+2+j].split()[5])
              item['xy'] = float(outlines[i+3+j].split()[1])
              item['yy'] = float(outlines[i+3+j].split()[3])
              item['zy'] = float(outlines[i+3+j].split()[5])
              item['xz'] = float(outlines[i+4+j].split()[1])
              item['yz'] = float(outlines[i+4+j].split()[3])
              item['zz'] = float(outlines[i+4+j].split()[5])
              item['eigenvalues'] = [float(outlines[i+5+j].split()[1]), float(outlines[i+5+j].split()[2]), float(outlines[i+5+j].split()[3])]
              NMR.append(item)
              j += 5

  return NMR, CARTESIANS,NATOMS,ATOMTYPES


def calculate_nics_parameters(qm_files,args,log,w_dir_initial,name_mol,lot,bs):

    directory = Path(w_dir_initial)

    total_data = pd.DataFrame()
    dist_data = pd.DataFrame()

    for counter,log in enumerate(qm_files):
        name = log.split('.log')[0]

        with open(log,'r') as F: 
            outlines =  F.readlines()

        NMR, cartesians,NATOMS,ATOMTYPES = getSHIELDING(outlines,args,log)
        vector_norm,centroid = update_coord(NATOMS,ATOMTYPES,cartesians,args,log,name_mol,w_dir_initial,'read')

        total_data.at[counter,'Name'] = name_mol
        total_data.at[counter,'log'] = name
        dist_data.at[counter,'Name'] = name_mol
        dist_data.at[counter,'log'] = name

        for tensor in NMR:
            if tensor['elementID'] == "Bq":
                at_idx = tensor['atom_index']
                
                nics_iso = tensor['isotropic']*-1
                w = np.array([tensor['xx'],tensor['yy'],tensor['zz']])
                nics_oop = -np.dot(w,vector_norm)
                nics_ip = (nics_iso * 3 - nics_oop)/2

                total_data.at[counter,f'iso-Bq-{at_idx}'] = nics_iso
                total_data.at[counter,f'oop-Bq-{at_idx}'] = nics_oop
                total_data.at[counter,f'ip-Bq-{at_idx}' ] = nics_ip

                v = np.array(cartesians[(at_idx-1)])
                dist = np.dot((v - centroid),vector_norm)
                dist_data.at[counter,str(at_idx)] = dist


    #creating folder for all molecules to write geom parameter
    
    if str(bs).find('/') > -1:
        folder = directory / f'/QPRED/nics/all_confs_nics/{lot}-{bs.split("/")[0]}'
    else:
        folder = directory / f'/QPRED/nics/all_confs_nics/{lot}-{bs}'

    folder.mkdir(parents=True,exist_ok=True)

    total_data.to_csv(folder/f'{name_mol}-all-nics-data.csv',index=False)
    dist_data.to_csv(folder/f'{name_mol}-all-nics-dist-data.csv',index=False)

def calculate_boltz_for_nics(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):

    # GoodVibes must be installed as a module (through pip or conda)
    cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
    for file in val:
        cmd_boltz.append(file)
    subprocess.call(cmd_boltz)

    #writing to coorect places
    if str(bs).find('/') > -1:
        destination = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs).split('/')[0]
    else:
        destination = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs)
    move_file_from_folder(destination, os.getcwd(),'Goodvibes_'+name+'.dat')

def calculate_avg_nics(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):
    if str(bs).find('/') > -1:
        nics_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+name+'-all-nics-data.csv'
        dist_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+name+'-all-nics-dist-data.csv'
        file = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs).split('/')[0]+'/Goodvibes_'+name+'.dat'
    else:
        nics_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs)+'/'+name+'-all-nics-data.csv'
        dist_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs)+'/'+name+'-all-nics-dist-data.csv'
        file = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs)+'/Goodvibes_'+name+'.dat'

    df_nics =  pd.read_csv(nics_file)
    df_dist =  pd.read_csv(dist_file)
    outlines = open(file,"r").readlines()

    #reading the data from boltz file
    for i in range(len(outlines)):
        # I remove the NMR from the file names using [0:-4]
        if outlines[i].find('   ***************************************************************************************************************************************\n') > -1 and outlines[i-1].find('   Structure') > -1:
            start_line = i+1
        elif outlines[i].find('   ***************************************************************************************************************************************\n') > -1:
            end_line = i

    boltz_values = pd.DataFrame()
    counter = 0
    for i in range(start_line,end_line):
        values = outlines[i].split()
        boltz_values.at[counter,'log'] = values[1]+'_nics'
        boltz_values.at[counter,'boltz'] = values[-1]
        counter+=1

    # multiply actual file and sum and add write to new FIL
    df_nics = df_nics.drop(['Name'], axis=1)
    df_dist = df_dist.drop(['Name'], axis=1)

    df_nics = df_nics.set_index('log')
    df_dist = df_dist.set_index('log')
    boltz_values =  boltz_values.set_index('log')

    avg_data_iso = pd.DataFrame()
    avg_data_oop = pd.DataFrame()
    avg_data_ip = pd.DataFrame()
    for col in df_nics.columns:
        df_nics[col] = df_nics[col].astype(float)*(boltz_values['boltz'].astype(float))
        if col.split('-')[0]=='oop':
            avg_data_oop.at[col,'oop'] = df_nics[col].sum()
            for col2 in df_dist.columns:
                if col.split('-')[2]==col2:
                    avg_data_oop.at[col,'dist'] = df_dist[col2].mean()
        if col.split('-')[0]=='ip':
            avg_data_ip.at[col,'ip'] = df_nics[col].sum()
            for col2 in df_dist.columns:
                if col.split('-')[2]==col2:
                    avg_data_ip.at[col,'dist'] = df_dist[col2].mean()
        if col.split('-')[0]=='iso':
            avg_data_iso.at[col,'iso'] = df_nics[col].sum()
            for col2 in df_dist.columns:
                if col.split('-')[2]==col2:
                    avg_data_iso.at[col,'dist'] = df_dist[col2].mean()

    avg_data_oop['Atom-Number'] = avg_data_oop.index
    avg_data_ip['Atom-Number'] = avg_data_ip.index
    avg_data_iso['Atom-Number'] = avg_data_iso.index
    if str(bs).find('/') > -1:
        folder = w_dir_initial + '/QPRED/nics/average_nics/'+str(lot)+'-'+str(bs).split('/')[0]
    else:
        folder = w_dir_initial + '/QPRED/nics/average_nics/'+str(lot)+'-'+str(bs)
    try:
        os.makedirs(folder)
    except OSError:
        if os.path.isdir(folder):
            pass

    avg_data_oop.to_csv(folder+'/'+name+'-avg-nics-oop.csv',index=False)
    avg_data_ip.to_csv(folder+'/'+name+'-avg-nics-ip.csv',index=False)
    avg_data_iso.to_csv(folder+'/'+name+'-avg-nics-iso.csv',index=False)
