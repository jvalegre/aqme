#####################################################.
#        This file stores all the functions         #
#           used for genrating a graph              #
#####################################################.

import os
from rdkit.Chem import AllChem as Chem
import numpy as np
from sklearn.metrics import mean_absolute_error
import statistics as stats
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import cclib

ev_2_kcal_mol = 23.061 #ev to kcal/mol
hartree_to_kcal = 627.509

def stats_calc(y_dft,y):
    y_diff = abs(np.subtract(y, y_dft))
    mae = mean_absolute_error(y, y_dft)
    sd = stats.stdev(y_diff)
    return mae,sd

def get_energy(inmols_min):
    energy_min = []
    for _,mol in enumerate(inmols_min):
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

def plot_graph(energy_rdkit,energy_min,energy_min_dft,lot,bs,energy_min_dft_sp,lot_sp,bs_sp,name_mol,type_csearch,cmin_option,w_dir_initial):

    list_all,mae_rdkit,sd_rdkit,mae_min,sd_min,mae_min,sd_min,mae_dft,sd_dft,x_axis_names,name_to_write = params_graph(energy_min_dft_sp,energy_min_dft,energy_min,energy_rdkit,lot,bs,lot_sp,bs_sp,type_csearch,cmin_option)

    plt.xticks(range(0,len(x_axis_names)), x_axis_names)
    fig=plt.figure() #Creates a new figure
    ax1=fig.add_subplot(111) #Plot with: 1 row, 1 column, first subplot.

    cmap = plt.cm.get_cmap('viridis', len(list_all))
    Path = mpath.Path
    x_axis = [0,1,2,3]
    for i,list in enumerate(list_all):
        index_axis = 0
        for j,_ in enumerate(range(len(list))):
            index_axis = add_patch_plot(x_axis,list,ax1,cmap,Path,index_axis,i)
            if j == len(list)-1:
                ax1.scatter(x_axis[:j+1],list,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)
                break

    plt.xticks(range(0,len(x_axis_names)), x_axis_names)

    y_margin = -0.05
    if len(energy_min_dft) != 0:
        textstr = r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[0], round(mae_rdkit, 2),round(sd_rdkit, 2))+'\n'
        textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[1],round(mae_min,2),round(sd_min,2)) +'\n'
        textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[2],round(mae_dft,2),round(sd_dft,2))
        if len(energy_min_dft_sp) != 0:
            y_margin = -0.03
        plt.figtext(0.5, y_margin, textstr, ha="center", fontsize=12,bbox=dict(facecolor='grey', alpha=0.25))

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

    graph_dir = Path(w_dir_initial + '/QSTAT/graph')
    graph_dir.mkdir(exist_ok=True, parents=True)
    os.chdir(graph_dir)

    plt.savefig(name_mol+'-'+name_to_write+'.png',bbox_inches='tight', format='png', dpi=400)
    plt.close()

    os.chdir(w_dir_initial)


def add_patch_plot(x_axis,list,ax1,cmap,Path,index_axis,i):
    '''
    Add patches to plots with multiple columns
    '''
    path_patch = mpatches.PathPatch(
            Path([(x_axis[index_axis], list[index_axis]), (x_axis[index_axis]+0.5, list[index_axis+1]), (x_axis[index_axis+1] ,list[index_axis+1])],
                [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
            fc="none", transform=ax1.transData, color=cmap(i))
    ax1.add_patch(path_patch)
    index_axis += 1

    return index_axis


def get_energy_graph(energy_min,energy_rdkit,name,energy_min_mae_sd,energy_rdkit_mae_sd,list_energies):
    '''
    Get the corresponding energies from CMIN and CSEARCH methods for energy graphs
    '''
    if energy_min is not None:
        for j,_ in enumerate(energy_min):
            if energy_min[j][0] == name:
                energy_min_mae_sd.append(float(energy_min[j][1]))
                list_energies.append(energy_min[j][1])
    if energy_rdkit is not None:
        for k,_ in enumerate(energy_rdkit):
            if energy_rdkit[k][0] == name:
                energy_rdkit_mae_sd.append(float(energy_rdkit[k][1]))
                list_energies.append(energy_min[k][1])

    return energy_min_mae_sd,energy_rdkit_mae_sd,list_energies


def graph(sdf_rdkit,sdf_xtb,sdf_ani,qm_files,sp_files,args,log,lot,bs,lot_sp,bs_sp,name_mol,w_dir_initial,w_dir_sp,w_dir,type):

    inmols_rdkit = Chem.SDMolSupplier(sdf_rdkit, removeHs=False)
    #get the energy from sdf
    energy_rdkit = get_energy(inmols_rdkit)
    energy_rdkit_sc = scaling_with_lowest(energy_rdkit)

    #get energy list for all conformers from SDFs of CMIN
    if os.path.exists(w_dir_initial+'/CSEARCH/xtb/'+name_mol+'_xtb.sdf'):
        sdf_mols = sdf_xtb
        sdf_source = 'xtb'
    if os.path.exists(w_dir_initial+'/CSEARCH/ani/'+name_mol+'_ani.sdf'):
        sdf_mols = sdf_ani
        sdf_source = 'ani'
    inmols_cmin =  Chem.SDMolSupplier(sdf_mols, removeHs=False)
    energy_cmin = get_energy(inmols_cmin)
    energy_cmin = rename_name(energy_cmin,sdf_source)
    energy_cmin_sc = scaling_with_lowest(energy_cmin)

    energy_dft,sp_dft = [],[]
    #get energy from log FILES
    if qm_files is not None:
        energy_dft,csearch_option,cmin_option = get_qm_energy_plot(type,qm_files,energy_dft)

    if sp_files is not None:
        os.chdir(w_dir_sp)
        sp_dft,csearch_option,cmin_option = get_qm_energy_plot(type,sp_files,sp_dft)
        os.chdir(w_dir)

    energy_dft_sc,sp_dft_sc = [],[]
    if qm_files is not None:
        energy_dft_sc = scaling_with_lowest(energy_dft)

    if sp_files is not None:
        sp_dft_sc = scaling_with_lowest(sp_dft)

    plot_graph(energy_rdkit_sc,energy_cmin_sc,energy_dft_sc,lot,bs,sp_dft_sc,lot_sp,bs_sp,name_mol,csearch_option,cmin_option,w_dir_initial)


def params_graph(energy_min_dft_sp,energy_min_dft,energy_min,energy_rdkit,lot,bs,lot_sp,bs_sp,type_csearch,cmin_option):
    if len(energy_min_dft_sp) != 0 or len(energy_min_dft) != 0:
        energy_dft_sp_mae_sd,energy_dft_mae_sd,energy_min_mae_sd,energy_rdkit_mae_sd = [],[],[],[]
        list_all = []
        if len(energy_min_dft_sp) == 0:
            energy_min_dft_sp = [None]
        for l,_ in enumerate(energy_min_dft_sp):
            if energy_min_dft_sp[l] is not None:
                list_energies = []
                name = energy_min_dft_sp[l][0]
                energy_dft_sp_mae_sd.append(float(energy_min_dft_sp[l][1]))
                list_energies.append(float(energy_min_dft_sp[l][1]))
            else:
                name = None
            for i,_ in enumerate(energy_min_dft):
                if name is not None:
                    if energy_min_dft[i][0] == name:
                        energy_dft_mae_sd.append(float(energy_min_dft[i][1]))
                        list_energies.append(float(energy_min_dft[i][1]))
                else:
                    name = energy_min_dft[i][0]
                    energy_dft_mae_sd.append(float(energy_min_dft[i][1]))
                    energy_min_mae_sd,energy_rdkit_mae_sd,list_energies = get_energy_graph(energy_min,energy_rdkit,name,energy_min_mae_sd,energy_rdkit_mae_sd,list_energies)
                    list_all.append(list_energies)
            if energy_min_dft_sp[l] is not None:
                energy_min_mae_sd,energy_rdkit_mae_sd,list_energies = get_energy_graph(energy_min,energy_rdkit,name,energy_min_mae_sd,energy_rdkit_mae_sd,list_energies)
                list_all.append(list_energies)

        if len(energy_min_dft_sp) != 0:
            mae_rdkit,sd_rdkit = stats_calc(energy_dft_sp_mae_sd,energy_rdkit_mae_sd)
            mae_min,sd_min = stats_calc(energy_dft_sp_mae_sd,energy_min_mae_sd)
            mae_dft,sd_dft = stats_calc(energy_dft_sp_mae_sd,energy_dft_mae_sd)
        else:
            mae_rdkit,sd_rdkit = stats_calc(energy_dft_mae_sd,energy_rdkit_mae_sd)
            mae_min,sd_min = stats_calc(energy_dft_mae_sd,energy_min_mae_sd)
            mae_dft,sd_dft = 0,0
        
        # setting some graphing options
        x_axis_names = [type_csearch]

        if cmin_option is not None:
            x_axis_names.append(cmin_option)
        if lot is not None:
            x_axis_names.append(lot+'_'+bs)
        if lot_sp is not None:
            x_axis_names.append(lot_sp+'_'+bs_sp)

        name_to_write = '-'.join([str(int(elem)) for elem in x_axis_names])

        return list_all,mae_rdkit,sd_rdkit,mae_min,sd_min,mae_min,sd_min,mae_dft,sd_dft,x_axis_names,name_to_write


def get_qm_energy_plot(type,plot_files,energy_dft):
    csearch_option,cmin_option = None,None
    for file in plot_files:
        if type == 'g16':
            data_sp = cclib.io.ccread(file)
            energy_qm = data_sp.scfenergies[0]*ev_2_kcal_mol
        elif type == 'orca':
            sp_lines = open(file,"r").readlines()
            for i,_ in reversed(range(sp_lines)):
                if sp_lines[i].find('FINAL SINGLE POINT ENERGY') > -1:
                    energy_qm = float(sp_lines[i].split()[-1])*hartree_to_kcal
                    break
        if len(file.split('_ani.sdf')) == 2 or len(file.split('_xtb.sdf')) == 2:
            name = file.replace('_ani.sdf','_xtb.sdf').split('_xtb.sdf')[0]
            if len(file.split('_ani.sdf')) == 2:
                cmin_option = 'ANI'
            elif len(file.split('_xtb.sdf')) == 2:
                cmin_option = 'xTB'
        else:
            if len(file.split('_summ.sdf')) == 2:
                csearch_option = 'SUMM'
            elif len(file.split('_fullmonte.sdf')) == 2:
                csearch_option = 'Fullmonte'
            elif len(file.split('_rdkit.sdf')) == 2:
                csearch_option = 'RDKit'
            name = file.replace('_summ.sdf','_rdkit.sdf').replace('_fullmonte.sdf','_rdkit.sdf').split('_rdkit.sdf')[0]
        
        energy_dft.append([name,energy_qm])

    return energy_dft,csearch_option,cmin_option