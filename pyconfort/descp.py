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
from pyconfort.csearch import getDihedralMatches


def get_data(rdkit_mols,min_mols,dft_mols,lot,bs,name_mol,args,type_csearch,type_min,w_dir_initial):
	geom_data = pd.DataFrame()
	for j, mol_j in enumerate(rdkit_mols):
		name = mol_j.GetProp('_Name')
		name_dft=  '_'.join(mol_j.GetProp('_Name').split())+'_'+type_min
		geom_data.at[j,'Name'] = name
		if len(args.dihedral) != 0:
			for d,dh in enumerate(args.dihedral):
				dihedral_rdkit = rdMolTransforms.GetDihedralDeg(mol_j.GetConformer(),dh[0],dh[1],dh[2],dh[3])
				geom_data.at[j,args.geom_par_name+'-Dihedral-'+type_csearch+'-'+str(dh[0])+'-'+str(dh[1])+'-'+str(dh[2])+'-'+str(dh[3])] = dihedral_rdkit
		if len(args.angle) != 0:
			for a,an in enumerate(args.angle):
				angle_rdkit = rdMolTransforms.GetAngleDeg(mol_j.GetConformer(),an[0],an[1],an[2])
				geom_data.at[j,args.geom_par_name+'-Angle-'+type_csearch+'-'+str(an[0])+'-'+str(an[1])+'-'+str(an[2])] = angle_rdkit
		if len(args.bond) != 0:
			for b,bd in enumerate(args.angle):
				bond_rdkit = rdMolTransforms.GetBondLength(mol_j.GetConformer(),bd[0],bd[1])
				geom_data.at[j,args.geom_par_name+'-Bond-'+type_csearch+'-'+str(bd[0])+'-'+str(bd[1])] = bond_rdkit
		if min_mols is not None:
			if type_min =='ani' or type_min=='xtb':
				for i, mol_i in enumerate(min_mols):
					if mol_i.GetProp('_Name') == name+' '+type_min:
						if len(args.dihedral) != 0:
							for d,dh in enumerate(args.dihedral):
								dihedral_min = rdMolTransforms.GetDihedralDeg(mol_i.GetConformer(),dh[0],dh[1],dh[2],dh[3])
								geom_data.at[j,args.geom_par_name+'-Dihedral-'+type_min+'-'+str(dh[0])+'-'+str(dh[1])+'-'+str(dh[2])+'-'+str(dh[3])] = dihedral_min
						if len(args.angle) != 0:
							for a,an in enumerate(args.angle):
								angle_min = rdMolTransforms.GetAngleDeg(mol_i.GetConformer(),an[0],an[1],an[2])
								geom_data.at[j,args.geom_par_name+'-Angle-'+type_min+'-'+str(an[0])+'-'+str(an[1])+'-'+str(an[2])] = angle_min
						if len(args.bond) != 0:
							for b,bd in enumerate(args.angle):
								bond_min = rdMolTransforms.GetBondLength(mol_i.GetConformer(),bd[0],bd[1])

		if dft_mols is not None:
			if type_min =='ani' or type_min=='xtb':
				for i, mol_i in enumerate(dft_mols):
					if mol_i.GetProp('_Name').split('/')[-1].split('.log')[0] == name_dft:
						if len(args.dihedral) != 0:
							for d,dh in enumerate(args.dihedral):
								dihedral_min = rdMolTransforms.GetDihedralDeg(mol_i.GetConformer(),dh[0],dh[1],dh[2],dh[3])
								geom_data.at[j,args.geom_par_name+'-Dihedral-'+lot+'-'+bs+'-'+str(dh[0])+'-'+str(dh[1])+'-'+str(dh[2])+'-'+str(dh[3])] = dihedral_min
						if len(args.angle) != 0:
							for a,an in enumerate(args.angle):
								angle_min = rdMolTransforms.GetAngleDeg(mol_i.GetConformer(),an[0],an[1],an[2])
								geom_data.at[j,args.geom_par_name+'-Angle-'+lot+'-'+bs+'-'+str(an[0])+'-'+str(an[1])+'-'+str(an[2])] = angle_min
						if len(args.bond) != 0:
							for b,bd in enumerate(args.angle):
								bond_min = rdMolTransforms.GetBondLength(mol_i.GetConformer(),bd[0],bd[1])
								geom_data.at[j,args.geom_par_name+'-Bond-'+lot+'-'+bs+'-'+str(bd[0])+'-'+str(bd[1])] = bond_min
	return geom_data


def calculate_parameters(sdf_rdkit,sdf_ani,sdf_xtb,log_files,args,log,w_dir_initial,name_mol,lot,bs):
	#creating folder for all molecules to write geom parameter
	folder = w_dir_initial + '/QSTAT/geom_parameters'
	try:
		os.makedirs(folder)
		os.chdir(folder)
	except OSError:
		if os.path.isdir(folder):
			os.chdir(folder)
		else:
			raise
	#get mol objects
	dft_mols= []
	rdkit_mols = Chem.SDMolSupplier(sdf_rdkit, removeHs=False)
	if args.rot_dihedral:
		args.dihedral = getDihedralMatches(rdkit_mols[0], args.heavyonly,log)
	if sdf_ani is not None:
		ani_mols = Chem.SDMolSupplier(sdf_ani, removeHs=False)
	if sdf_xtb is not None:
		xtb_mols = Chem.SDMolSupplier(sdf_xtb, removeHs=False)
	ob_compat = True
	try:
		import openbabel as ob
	except (ModuleNotFoundError,AttributeError):
		ob_compat = False
		log.write('\nx  Open Babel is not installed correctly, it is not possible to get molecular descriptors')

	if ob_compat:
		obConversion = ob.OBConversion()
		obConversion.SetInAndOutFormats("log", "mol")
		ob_mol = ob.OBMol()
		for file in log_files:
			if str(bs).find('/') > -1:
				obConversion.ReadFile(ob_mol, args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files/'+file)
				obConversion.WriteFile(ob_mol, args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files/'+file.split('.')[0]+'.mol')
				obConversion.CloseOutFile()
				dft_mols.append(Chem.MolFromMolFile(args.path + str(lot) + '-' + str(bs).split('/')[0] +'/success/output_files/'+file.split('.')[0]+'.mol', removeHs=False))
			else:
				obConversion.ReadFile(ob_mol, args.path + str(lot) + '-' + str(bs) +'/success/output_files/'+file)
				obConversion.WriteFile(ob_mol, args.path + str(lot) + '-' + str(bs) +'/success/output_files/'+file.split('.')[0]+'.mol')
				obConversion.CloseOutFile()
				dft_mols.append(Chem.MolFromMolFile(args.path + str(lot) + '-' + str(bs) +'/success/output_files/'+file.split('.')[0]+'.mol', removeHs=False))

		if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf'):
			geom_data = get_data(rdkit_mols,xtb_mols,dft_mols,lot,bs,name_mol,args,'rdkit','xtb',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-rdkit-xtb.csv',index=False)

		if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf'):
			geom_data = get_data(rdkit_mols,ani_mols,dft_mols,lot,bs,name_mol,args,'rdkit','ani',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-rdkit-ani.csv',index=False)

			##########

		if  os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf'):
			geom_data = get_data(rdkit_mols,xtb_mols,dft_mols,lot,bs,name_mol,args,'summ','xtb',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-summ-xtb.csv',index=False)

		if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf'):
			geom_data = get_data(rdkit_mols,ani_mols,dft_mols,lot,bs,name_mol,args,'summ','ani',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-summ-ani.csv',index=False)

			#############

		if  os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf'):
			geom_data = get_data(rdkit_mols,xtb_mols,dft_mols,lot,bs,name_mol,args,'fullmonte','xtb',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-fullmonte-xtb.csv',index=False)


		if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf'):
			geom_data = get_data(rdkit_mols,ani_mols,dft_mols,lot,bs,name_mol,args,'fullmonte','ani',w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-fullmonte-ani.csv',index=False)

			############

		if os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf'):
			geom_data = get_data(rdkit_mols,None,dft_mols,lot,bs,name_mol,args,'summ',None,w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-summ.csv',index=False)


		if os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') :
			geom_data = get_data(rdkit_mols,None,dft_mols,lot,bs,name_mol,args,'rdkit',None,w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-rdkit.csv',index=False)

		if os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') :
			geom_data = get_data(rdkit_mols,None,dft_mols,lot,bs,name_mol,args,'rdkit',None,w_dir_initial)
			geom_data.to_csv(name_mol+'-all-geom-data-with-fullmonte.csv',index=False)
