#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	       used for genrating all parameters        #
#####################################################.

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms
import os
import numpy as np
import pandas as pd
from pyconfort.confgen_functions import getDihedralMatches


def get_data(rdkit_mols,min_mols,type,args):
    geom_data = pd.DataFrame()
    for j, mol_j in enumerate(rdkit_mols):
        name = mol_j.GetProp('_Name')
        geom_data.at[j,'Name'] = name
        geom_data.at[j,'Energy-rdkit'] = mol_j.GetProp('Energy')
        if len(args.dihedral) != 0:
            for d,dh in enumerate(args.dihedral):
                dihedral_rdkit = rdMolTransforms.GetDihedralDeg(mol_j.GetConformer(),dh[0],dh[1],dh[2],dh[3])
                geom_data.at[j,args.geom_par_name+'-Dihedral-rdkit-'+str(dh[0])+'-'+str(dh[1])+'-'+str(dh[2])+'-'+str(dh[3])] = dihedral_rdkit
        if len(args.angle) != 0:
            for a,an in enumerate(args.angle):
                angle_rdkit = rdMolTransforms.GetAngleDeg(mol_j.GetConformer(),an[0],an[1],an[2])
                geom_data.at[j,args.geom_par_name+'-Angle-rdkit-'+str(an[0])+'-'+str(an[1])+'-'+str(an[2])] = angle_rdkit
        if len(args.bond) != 0:
            for b,bd in enumerate(args.angle):
                bond_rdkit = rdMolTransforms.GetBondLength(mol_j.GetConformer(),bd[0],bd[1])
                geom_data.at[j,args.geom_par_name+'-Bond-rdkit-'+str(bd[0])+'-'+str(bd[1])] = bond_rdkit
        if min_mols is not None:
            if type =='ani' or type=='xtb':
                for i, mol_i in enumerate(min_mols):
                    if mol_i.GetProp('_Name') == name+' '+type:
                        geom_data.at[j,'Energy-'+type] = mol_i.GetProp('Energy')
                        if len(args.dihedral) != 0:
                            for d,dh in enumerate(args.dihedral):
                                dihedral_min = rdMolTransforms.GetDihedralDeg(mol_i.GetConformer(),dh[0],dh[1],dh[2],dh[3])
                                geom_data.at[j,args.geom_par_name+'-Dihedral-'+type+'-'+str(dh[0])+'-'+str(dh[1])+'-'+str(dh[2])+'-'+str(dh[3])] = dihedral_min
                        if len(args.angle) != 0:
                            for a,an in enumerate(args.angle):
                                angle_min = rdMolTransforms.GetAngleDeg(mol_i.GetConformer(),an[0],an[1],an[2])
                                geom_data.at[j,args.geom_par_name+'-Angle-'+type+'-'+str(an[0])+'-'+str(an[1])+'-'+str(an[2])] = angle_min
                        if len(args.bond) != 0:
                            for b,bd in enumerate(args.angle):
                                bond_min = rdMolTransforms.GetBondLength(mol_i.GetConformer(),bd[0],bd[1])
                                geom_data.at[j,args.geom_par_name+'-Bond-'+type+'-'+str(bd[0])+'-'+str(bd[1])] = bond_min
    return geom_data


def calculate_parameters(sdf_rdkit,sdf_ani,sdf_xtb,args,log,w_dir_initial,name_mol):
    #creating folder for all molecules to write geom parameter
    folder = w_dir_initial + '/geom_parameters'
    try:
        os.makedirs(folder)
        os.chdir(folder)
    except OSError:
        if os.path.isdir(folder):
            os.chdir(folder)
        else:
            raise
    #get mol objects
    rdkit_mols = Chem.SDMolSupplier(sdf_rdkit, removeHs=False)
    if args.rot_dihedral:
        args.dihedral = getDihedralMatches(rdkit_mols[0], args.heavyonly,log)
    if sdf_ani is not None:
        ani_mols = Chem.SDMolSupplier(sdf_ani, removeHs=False)
    if sdf_xtb is not None:
        xtb_mols = Chem.SDMolSupplier(sdf_xtb, removeHs=False)

    #loop over only uniques
    if args.xtb:
        geom_data = get_data(rdkit_mols,xtb_mols,'xtb',args)
        geom_data.to_csv(name_mol+'-all-geom-data-with-xtb.csv',index=False)
    if args.ANI1ccx:
        geom_data = get_data(rdkit_mols,ani_mols,'ani',args)
        geom_data.to_csv(name_mol+'-all-geom-data-with-ani.csv',index=False)
    if not args.xtb and not args.ANI1ccx:
        geom_data = get_data(rdkit_mols,None,'rdkit-only',args)
        geom_data.to_csv(name_mol+'-all-geom-data-with-rdkit-only.csv',index=False)
