#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used for genrating a graph      	    #
#####################################################.

import os
from pyconfort.cmin import rdkit_sdf_read
from rdkit.Chem import AllChem as Chem
import numpy as np
from sklearn.metrics import mean_absolute_error
import statistics as stats

ev_2_kcal_mol = 23.061 #ev to kcal/mol
hartree_to_kcal = 627.509

def stats_calc(y_dft,y):
	y_diff = abs(np.subtract(y, y_dft))
	mae = mean_absolute_error(y, y_dft)
	sd = stats.stdev(y_diff)
	return mae,sd

def get_energy(inmols_min):
	energy_min = []
	for i,mol in enumerate(inmols_min):
		energy_min.append(["_".join(mol.GetProp('_Name').split(' ')),mol.GetProp('Energy')])
	return energy_min

def rename_name(energy,type):
	for i,_ in enumerate(energy):
		energy[i][0] = energy[i][0].split('_'+type)[0]
	return energy

def scaling_with_lowest(energy):
	#scaling arrays
	v = np.array(energy)[:, 1].astype(np.float)
	energy_sc = (v - v.min())
	for i,_ in enumerate(energy):
		energy[i][1] = energy_sc[i]
	return energy

def plot_graph(energy_rdkit,energy_min,energy_min_dft,lot,bs,energy_min_dft_sp,lot_sp,bs_sp,name_mol,args,type_csearch,type_min,w_dir_initial):

	try:
		import matplotlib.path as mpath
		import matplotlib.patches as mpatches
		import matplotlib.pyplot as plt
	except (ModuleNotFoundError,AttributeError):
	    log.write('The Matplotlib module is not installed correctly - the graphing functions are not available')

	def get_cmap(n, name='viridis'):
		'''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
		RGB color; the keyword argument name must be a standard mpl colormap name.'''
		return plt.cm.get_cmap(name, n)

	if len(energy_min_dft_sp)  != 0:
		energy_dft_sp_mae_sd,energy_dft_mae_sd,energy_min_mae_sd,energy_rdkit_mae_sd = [],[],[],[]
		for l,_ in enumerate(energy_min_dft_sp):
			name = energy_min_dft_sp[l][0]
			energy_dft_sp_mae_sd.append(float(energy_min_dft_sp[l][1]))
			for i,_ in enumerate(energy_min_dft):
				if energy_min_dft[i][0] == name:
					energy_dft_mae_sd.append(float(energy_min_dft[i][1]))
			if energy_min is not None:
				for j,_ in enumerate(energy_min):
					if energy_min[j][0] == name:
						energy_min_mae_sd.append(float(energy_min[j][1]))
			if energy_rdkit is not None:
				for k,_ in enumerate(energy_rdkit):
					if energy_rdkit[k][0] == name:
						energy_rdkit_mae_sd.append(float(energy_rdkit[k][1]))

		mae_rdkit,sd_rdkit = stats_calc(energy_dft_sp_mae_sd,energy_rdkit_mae_sd)
		mae_min,sd_min = stats_calc(energy_dft_sp_mae_sd,energy_min_mae_sd)
		mae_dft,sd_dft = stats_calc(energy_dft_sp_mae_sd,energy_dft_mae_sd)

	elif len(energy_min_dft) != 0:
		energy_dft_mae_sd,energy_min_mae_sd,energy_rdkit_mae_sd = [],[],[]
		for i,_ in enumerate(energy_min_dft):
			name = energy_min_dft[i][0]
			energy_dft_mae_sd.append(float(energy_min_dft[i][1]))
			if energy_min is not None:
				for j,_ in enumerate(energy_min):
					if energy_min[j][0] == name:
						energy_min_mae_sd.append(float(energy_min[j][1]))
			if energy_rdkit is not None:
				for k,_ in enumerate(energy_rdkit):
					if energy_rdkit[k][0] == name:
						energy_rdkit_mae_sd.append(float(energy_rdkit[k][1]))

		mae_rdkit,sd_rdkit = stats_calc(energy_dft_mae_sd,energy_rdkit_mae_sd)
		mae_min,sd_min = stats_calc(energy_dft_mae_sd,energy_min_mae_sd)


	x_axis_names = []
	x_axis_names.append(type_csearch)
	if type_min is not None:
		x_axis_names.append(type_min)
	if lot is not None:
		x_axis_names.append(lot+'\n'+bs)
	if lot_sp is not None:
		x_axis_names.append(lot_sp+'\n'+bs_sp)

	x_axis_names_mae = []
	x_axis_names_mae.append(type_csearch)
	if type_min is not None:
		x_axis_names_mae.append(type_min)
	if lot is not None:
		x_axis_names_mae.append(lot+'-'+bs)
	if lot_sp is not None:
		x_axis_names_mae.append(lot_sp+'-'+bs_sp)

	x_axis = [0,1,2,3]
	x_axis_1 = [0,1,2]
	x_axis_2 = [0,1]
	x_axis_3 = [0]


	list_all = []
	name_all = []
	for i,_ in enumerate(energy_rdkit):
		list = []
		name = energy_rdkit[i][0]
		name_all.append(name)
		list.append(energy_rdkit[i][1])
		if energy_min is not None:
			for j,_ in enumerate(energy_min):
				if energy_min[j][0] == name:
					list.append(energy_min[j][1])
		if energy_min_dft is not None:
			for k,_ in enumerate(energy_min_dft):
				if energy_min_dft[k][0] == name:
					list.append(energy_min_dft[k][1])
		if energy_min_dft_sp is not None:
			for l,_ in enumerate(energy_min_dft_sp):
				if energy_min_dft_sp[l][0] == name:
					list.append(energy_min_dft_sp[l][1])
		list_all.append(list)

	fig=plt.figure() #Creates a new figure
	ax1=fig.add_subplot(111) #Plot with: 1 row, 1 column, first subplot.

	cmap = get_cmap(len(list_all))
	Path = mpath.Path
	for i, list in enumerate(list_all):
		if len(list) == 2 or len(list) == 3 or len(list) == 4:
			path_patch_1 = mpatches.PathPatch(
					Path([(x_axis[0], list[0]), (x_axis[0]+0.5, list[1]), (x_axis[1] ,list[1])],
						 [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
					fc="none", transform=ax1.transData, color=cmap(i))
			ax1.add_patch(path_patch_1)
		if len(list) == 3 or len(list) == 4:
			path_patch_2 = mpatches.PathPatch(
					Path([(x_axis[1], list[1]), (x_axis[1]+0.5, list[2]), (x_axis[2] ,list[2])],
						 [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
					fc="none", transform=ax1.transData, color=cmap(i))
			ax1.add_patch(path_patch_2)
		if len(list) == 4:
			path_patch_3 = mpatches.PathPatch(
					Path([(x_axis[2], list[2]), (x_axis[2]+0.5, list[3]), (x_axis[3] ,list[3])],
						 [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
					fc="none", transform=ax1.transData, color=cmap(i))
			ax1.add_patch(path_patch_3)
		if len(list) == 4:
			ax1.scatter(x_axis,list,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)
		elif len(list) == 3:
			ax1.scatter(x_axis_1,list,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)
		elif len(list) == 2:
			ax1.scatter(x_axis_2,list,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)
		elif len(list) == 1:
			ax1.scatter(x_axis_3,list,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)


	plt.xticks(range(0,1), x_axis_names)
	if type_min is not None:
		plt.xticks(range(0,2), x_axis_names)
	if lot is not None:
		plt.xticks(range(0,3), x_axis_names)
	if lot_sp is not None:
		plt.xticks(range(0,4), x_axis_names)

	if len(energy_min_dft_sp) != 0:
		textstr = r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names_mae[0], round(mae_rdkit, 2),round(sd_rdkit, 2))+'\n'
		textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names_mae[1],round(mae_min,2),round(sd_min,2)) +'\n'
		textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names_mae[2],round(mae_dft,2),round(sd_dft,2))
		plt.figtext(0.5, -0.05, textstr ,  ha="center", fontsize=12,bbox=dict(facecolor='grey', alpha=0.25))
	elif len(energy_min_dft) != 0:
		textstr = r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names_mae[0], round(mae_rdkit, 2),round(sd_rdkit, 2))+'\n'
		textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names_mae[1],round(mae_min,2),round(sd_min,2))
		plt.figtext(0.5, -0.03, textstr ,  ha="center", fontsize=12,bbox=dict(facecolor='grey', alpha=0.25))
	#ax1.legend(lines,labels,loc='upper center', prop={'size':4}, bbox_to_anchor=(0.5, -0.13), fancybox=True, shadow=True, ncol=5)
	ax1.set_xlabel('Type of Calculation',fontsize=10)
	ax1.set_ylabel('Relative Energy (kcal/mol)',fontsize=10)
	#plt.setp(ax1.get_xticklabels(), rotation=60, ha="right", visible=True)
	plt.grid(linestyle='--', linewidth=0.75)
	plt.setp(ax1.get_xticklabels(), rotation=0, visible=True)
	title_string=('Energies of Conformers for different Methods : {0}'.format(name_mol))
	ax1.set_title(title_string, fontsize=12)
	fig.tight_layout()
	fig.subplots_adjust(top=0.92,bottom=0.2)
	#creating folder for all molecules to write geom parameter
	folder = w_dir_initial + '/QSTAT/graph'
	try:
		os.makedirs(folder)
		os.chdir(folder)
	except OSError:
		if os.path.isdir(folder):
			os.chdir(folder)
		else:
			raise

	name_to_write = type_csearch
	if type_min is not None:
		name_to_write += '-'+type_min
	if lot is not None:
		name_to_write += '-'+lot+'-'+bs
	if lot_sp is not None:
		name_to_write += '-'+lot_sp+'-'+bs_sp
	plt.savefig(name_mol+'-'+name_to_write+'.png',bbox_inches='tight', format='png', dpi=400)
	plt.close()
	os.chdir(w_dir_initial)

def graph(sdf_rdkit,sdf_xtb,sdf_ani,log_files,sp_files,args,log,lot,bs,lot_sp,bs_sp,name_mol,w_dir_initial,w_dir_sp,w_dir,type):

	try:
		import cclib
	except (ModuleNotFoundError,AttributeError):
		log.write('The cclib module is not installed correctly - the graph funtcion is not available')

	inmols_rdkit = Chem.SDMolSupplier(sdf_rdkit, removeHs=False)
	#get the energy from sdf
	energy_rdkit = get_energy(inmols_rdkit)
	energy_rdkit_sc = scaling_with_lowest(energy_rdkit)

	if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf'):
		#get energy list for all conformers from sdfs of rdkit and minimize
		inmols_xtb =  Chem.SDMolSupplier(sdf_xtb, removeHs=False)
		energy_xtb = get_energy(inmols_xtb)
		energy_xtb = rename_name(energy_xtb,'xtb')
		energy_xtb_sc = scaling_with_lowest(energy_xtb)
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf'):
		#get energy list for all conformers from sdfs of rdkit and minimize
		inmols_ani = Chem.SDMolSupplier(sdf_ani, removeHs=False)
		energy_ani = get_energy(inmols_ani)
		energy_ani = rename_name(energy_ani,'ani')
		energy_ani_sc = scaling_with_lowest(energy_ani)

	energy_rdkit_dft,energy_xtb_dft,energy_ani_dft = [],[],[]
	sp_rdkit_dft,sp_xtb_dft,sp_ani_dft = [],[],[]
	#get energy from log FILES
	if log_files is not None:
		for file in log_files:
			data = cclib.io.ccread(file)
			if len(file.split('_ani.log')) == 2 or len(file.split('_xtb.log')) == 2:
				if len(file.split('_ani.log')) == 2:
					name = file.split('_ani.log')[0]
					energy_ani_dft.append([name,data.scfenergies[0]*ev_2_kcal_mol])
				if len(file.split('_xtb.log')) == 2:
					name = file.split('_xtb.log')[0]
					energy_xtb_dft.append([name,data.scfenergies[0]*ev_2_kcal_mol])
			else:
				name = file.split('.log')[0]
				energy_rdkit_dft.append([name,data.scfenergies[0]*ev_2_kcal_mol])

	if sp_files is not None:
		os.chdir(w_dir_sp)
		if type=='g16':
			for file in sp_files:
				data_sp = cclib.io.ccread(file)
				if len(file.split('_ani_'+args.suffix_sp+'.log')) == 2 or len(file.split('_xtb_'+args.suffix_sp+'.log')) == 2:
					if len(file.split('_ani_'+args.suffix_sp+'.log')) == 2:
						name = file.split('_ani_'+args.suffix_sp+'.log')[0]
						sp_ani_dft.append([name,data_sp.scfenergies[0]*ev_2_kcal_mol])
					if len(file.split('_xtb_'+args.suffix_sp+'.log')) == 2:
						name = file.split('_xtb_'+args.suffix_sp+'.log')[0]
						sp_xtb_dft.append([name,data_sp.scfenergies[0]*ev_2_kcal_mol])
				else:
					name = file.split('.log')[0]
					sp_rdkit_dft.append([name,data_sp.scfenergies[0]*ev_2_kcal_mol])
		if type=='orca':
			for file in sp_files:
				sp_lines = open(file,"r").readlines()
				for i,line in enumerate(sp_lines):
					if sp_lines[i].find('FINAL SINGLE POINT ENERGY') > -1:
						if len(file.split('_ani_'+args.suffix_sp+'.out')) == 2 or len(file.split('_xtb_'+args.suffix_sp+'.out')) == 2:
							if len(file.split('_ani_'+args.suffix_sp+'.out')) == 2:
								name = file.split('_ani_'+args.suffix_sp+'.out')[0]
								sp_ani_dft.append([name,float(sp_lines[i].split()[-1])*hartree_to_kcal])
							if len(file.split('_xtb_'+args.suffix_sp+'.out')) == 2:
								name = file.split('_xtb_'+args.suffix_sp+'.out')[0]
								sp_xtb_dft.append([name,float(sp_lines[i].split()[-1])*hartree_to_kcal])
						else:
							name = file.split('.out')[0]
							sp_rdkit_dft.append([name,float(sp_lines[i].split()[-1])*hartree_to_kcal])
	os.chdir(w_dir)


	energy_ani_dft_sc, energy_xtb_dft_sc, energy_rdkit_dft_sc = [],[],[]
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and log_files is not None:
		energy_ani_dft_sc = scaling_with_lowest(energy_ani_dft)
	if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and log_files is not None:
		energy_xtb_dft_sc = scaling_with_lowest(energy_xtb_dft)
	if not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and log_files is not None:
		energy_rdkit_dft_sc = scaling_with_lowest(energy_rdkit_dft)


	sp_ani_dft_sc, sp_xtb_dft_sc, sp_rdkit_dft_sc = [],[],[]
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and sp_files is not None:
		sp_ani_dft_sc = scaling_with_lowest(sp_ani_dft)
	if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and sp_files is not None:
		sp_xtb_dft_sc = scaling_with_lowest(sp_xtb_dft)
	if not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and sp_files is not None:
		sp_rdkit_dft_sc = scaling_with_lowest(sp_rdkit_dft)


	if  os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf'):
		plot_graph(energy_rdkit_sc,energy_xtb_sc,energy_xtb_dft_sc,lot,bs,sp_xtb_dft_sc,lot_sp,bs_sp,name_mol,args,'RDKit','xTB',w_dir_initial)
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf'):
		plot_graph(energy_rdkit_sc,energy_ani_sc,energy_ani_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'RDKit','ANI',w_dir_initial)

	if  os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf'):
		plot_graph(energy_rdkit_sc,energy_xtb_sc,energy_xtb_dft_sc,lot,bs,sp_xtb_dft_sc,lot_sp,bs_sp,name_mol,args,'SUMM','xTB',w_dir_initial)
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf'):
		plot_graph(energy_rdkit_sc,energy_ani_sc,energy_ani_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'SUMM','ANI',w_dir_initial)

	if  os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf'):
		plot_graph(energy_rdkit_sc,energy_xtb_sc,energy_xtb_dft_sc,lot,bs,sp_xtb_dft_sc,lot_sp,bs_sp,name_mol,args,'Fullmonte','xTB',w_dir_initial)
	if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') and os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf'):
		plot_graph(energy_rdkit_sc,energy_ani_sc,energy_ani_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'Fullmonte','ANI',w_dir_initial)

	if os.path.exists(w_dir_initial+'/CSEARCH/summ/'+name_mol+'_summ.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf'):
		plot_graph(energy_rdkit_sc,None,energy_rdkit_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'SUMM',None,w_dir_initial)
	if os.path.exists(w_dir_initial+'/CSEARCH/rdkit/'+name_mol+'_rdkit.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') :
		plot_graph(energy_rdkit_sc,None,energy_rdkit_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'RDKit',None,w_dir_initial)
	if os.path.exists(w_dir_initial+'/CSEARCH/fullmonte/'+name_mol+'_fullmonte.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf') and not os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf') :
		plot_graph(energy_rdkit_sc,None,energy_rdkit_dft_sc,lot,bs,sp_ani_dft_sc,lot_sp,bs_sp,name_mol,args,'Fullmonte',None,w_dir_initial)
