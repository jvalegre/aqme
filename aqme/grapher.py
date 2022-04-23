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

def rename_name(energy,module_type):
    for i,_ in enumerate(energy):
        energy[i][0] = energy[i][0].split('_'+module_type)[0]
    return energy

def scaling_with_lowest(energy):
    #scaling arrays
    v = np.array(energy)[:, 1].astype(np.float)
    energy_sc = (v - v.min())
    for i,_ in enumerate(energy):
        energy[i][1] = energy_sc[i]
    return energy

def plot_graph(qm_data,file_data,csearch_cmin_data):

    list_all,mae_rdkit,sd_rdkit,mae_min,sd_min,mae_min,sd_min,mae_dft,sd_dft,x_axis_names,name_to_write = params_graph(qm_data,csearch_cmin_data)

    plt.xticks(range(0,len(x_axis_names)), x_axis_names)
    fig=plt.figure() #Creates a new figure
    ax1=fig.add_subplot(111) #Plot with: 1 row, 1 column, first subplot.

    cmap = plt.cm.get_cmap('viridis', len(list_all))
    Path = mpath.Path
    x_axis = [0,1,2,3]
    for i,list_level in enumerate(list_all):
        index_axis = 0
        for j,_ in enumerate(range(len(list_level))):
            path_patch = mpatches.PathPatch(
            Path([(x_axis[index_axis], list_level[index_axis]), (x_axis[index_axis]+0.5, list_level[index_axis+1]), (x_axis[index_axis+1] ,list_level[index_axis+1])],
                [Path.MOVETO, Path.CURVE3, Path.CURVE3]), fc="none", transform=ax1.transData, color=cmap(i))
            index_axis += 1
            ax1.add_patch(path_patch)
            if j == len(list_level)-1:
                ax1.scatter(x_axis[:j+1],list_level,color=cmap(i), marker='o',zorder=2,edgecolors= "black",linewidth=0.5)
                break

    plt.xticks(range(0,len(x_axis_names)), x_axis_names)

    y_margin = -0.05
    if len(qm_data['Energy dft']) != 0:
        textstr = r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[0], round(mae_rdkit, 2),round(sd_rdkit, 2))+'\n'
        textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[1],round(mae_min,2),round(sd_min,2)) +'\n'
        textstr += r'{0} = {1} $\pm$ {2} (kcal/mol)'.format(x_axis_names[2],round(mae_dft,2),round(sd_dft,2))
        if len(qm_data['Energy dft SP']) != 0:
            y_margin = -0.03
        plt.figtext(0.5, y_margin, textstr, ha="center", fontsize=12,bbox=dict(facecolor='grey', alpha=0.25))

    #ax1.legend(lines,labels,loc='upper center', prop={'size':4}, bbox_to_anchor=(0.5, -0.13), fancybox=True, shadow=True, ncol=5)
    ax1.set_xlabel('Type of Calculation',fontsize=10)
    ax1.set_ylabel('Relative Energy (kcal/mol)',fontsize=10)
    #plt.setp(ax1.get_xticklabels(), rotation=60, ha="right", visible=True)
    plt.grid(linestyle='--', linewidth=0.75)
    plt.setp(ax1.get_xticklabels(), rotation=0, visible=True)
    title_string=('Energies of Conformers for different Methods : {0}'.format(file_data['Name mol']))
    ax1.set_title(title_string, fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(top=0.92,bottom=0.2)

    graph_dir = Path(file_data['Initial dir'] + '/qstat/graph')
    graph_dir.mkdir(exist_ok=True, parents=True)
    os.chdir(graph_dir)

    plt.savefig(file_data['Name mol']+'-'+name_to_write+'.png',bbox_inches='tight', format='png', dpi=400)
    plt.close()

    os.chdir(file_data['Initial dir'])


def get_energy_graph(csearch_cmin_data,name,list_energies):
    '''
    Get the corresponding energies from CMIN and CSEARCH methods for energy graphs
    '''

    energy_min_mae_sd,energy_rdkit_mae_sd = [],[]
    if csearch_cmin_data['Energy and name CMIN scaled'] is not None:
        for j,_ in enumerate(csearch_cmin_data['Energy and name CMIN scaled']):
            if csearch_cmin_data['Energy and name CMIN scaled'][j][0] == name:
                energy_min_mae_sd.append(float(csearch_cmin_data['Energy and name CMIN scaled'][j][1]))
                list_energies.append(csearch_cmin_data['Energy and name CMIN scaled'][j][1])
    if csearch_cmin_data['Energy and name CSEARCH scaled'] is not None:
        for k,_ in enumerate(csearch_cmin_data['Energy and name CSEARCH scaled']):
            if csearch_cmin_data['Energy and name CSEARCH scaled'][k][0] == name:
                energy_rdkit_mae_sd.append(float(csearch_cmin_data['Energy and name CSEARCH scaled'][k][1]))
                list_energies.append(csearch_cmin_data['Energy and name CSEARCH scaled'][k][1])

    return energy_min_mae_sd,energy_rdkit_mae_sd,list_energies


def graph(qm_data,file_data,csearch_cmin_data):

    inmols_rdkit = Chem.SDMolSupplier(csearch_cmin_data['SDFs CSEARCH'], removeHs=False, sanitize=False)
    #get the energy from sdf
    energy_rdkit = get_energy(inmols_rdkit)
    energy_rdkit_sc = scaling_with_lowest(energy_rdkit)

    #get energy list for all conformers from SDFs of CMIN
    if os.path.exists(file_data['Initial dir']+'/CSEARCH/xtb/'+file_data['Name mol']+'_xtb.sdf'):
        sdf_mols = csearch_cmin_data['SDFs xTB']
        sdf_source = 'xtb'
    if os.path.exists(file_data['Initial dir']+'/CSEARCH/ani/'+file_data['Name mol']+'_ani.sdf'):
        sdf_mols = csearch_cmin_data['SDFs ANI']
        sdf_source = 'ani'
    inmols_cmin =  Chem.SDMolSupplier(sdf_mols, removeHs=False, sanitize=False)
    energy_cmin = get_energy(inmols_cmin)
    energy_cmin = rename_name(energy_cmin,sdf_source)
    energy_cmin_sc = scaling_with_lowest(energy_cmin)

    energy_dft,energy_dft_sp = [],[]
    energy_dft_sc,energy_dft_sp_sc = [],[]
    #get energy from log FILES
    if qm_data['QM files'] is not None:
        energy_dft,type_csearch,type_cmin = get_qm_energy_plot(type,qm_data['QM files'],energy_dft)
        energy_dft_sc = scaling_with_lowest(energy_dft)

    if qm_data['QM files SP'] is not None:
        os.chdir(file_data['Dir SP'])
        energy_dft_sp,type_csearch,type_cmin = get_qm_energy_plot(type,qm_data['QM files SP'],energy_dft_sp)
        energy_dft_sp_sc = scaling_with_lowest(energy_dft_sp)
        os.chdir(file_data['Working dir'])

    qm_data['Energy dft'] = energy_dft_sc
    qm_data['Energy dft SP'] = energy_dft_sp_sc
    csearch_cmin_data['Energy and name CSEARCH scaled'] = energy_rdkit_sc
    csearch_cmin_data['Energy and name CMIN scaled'] = energy_cmin_sc
    csearch_cmin_data['CSEARCH type'] = type_csearch
    csearch_cmin_data['CMIN type'] = type_cmin

    plot_graph(qm_data,file_data,csearch_cmin_data)


def params_graph(qm_data,csearch_cmin_data):
    if len(qm_data['Energy dft SP']) != 0 or len(qm_data['Energy dft']) != 0:
        energy_dft_sp_mae_sd,energy_dft_mae_sd = [],[]
        list_all = []
        if len(qm_data['Energy dft SP']) == 0:
            qm_data['Energy dft SP'] = [None]
        for l,_ in enumerate(qm_data['Energy dft SP']):
            if qm_data['Energy dft SP'][l] is not None:
                list_energies = []
                name = qm_data['Energy dft SP'][l][0]
                energy_dft_sp_mae_sd.append(float(qm_data['Energy dft SP'][l][1]))
                list_energies.append(float(qm_data['Energy dft SP'][l][1]))
            else:
                name = None
            for i,_ in enumerate(qm_data['Energy dft']):
                if name is not None:
                    if qm_data['Energy dft'][i][0] == name:
                        energy_dft_mae_sd.append(float(qm_data['Energy dft'][i][1]))
                        list_energies.append(float(qm_data['Energy dft'][i][1]))
                else:
                    name = qm_data['Energy dft'][i][0]
                    energy_dft_mae_sd.append(float(qm_data['Energy dft'][i][1]))
                    energy_min_mae_sd,energy_rdkit_mae_sd,list_energies = get_energy_graph(csearch_cmin_data,name,list_energies)
                    list_all.append(list_energies)
            if qm_data['Energy dft SP'][l] is not None:
                energy_min_mae_sd,energy_rdkit_mae_sd,list_energies = get_energy_graph(csearch_cmin_data,name,list_energies)
                list_all.append(list_energies)

        if len(qm_data['Energy dft SP']) != 0:
            mae_rdkit,sd_rdkit = stats_calc(energy_dft_sp_mae_sd,energy_rdkit_mae_sd)
            mae_min,sd_min = stats_calc(energy_dft_sp_mae_sd,energy_min_mae_sd)
            mae_dft,sd_dft = stats_calc(energy_dft_sp_mae_sd,energy_dft_mae_sd)
        else:
            mae_rdkit,sd_rdkit = stats_calc(energy_dft_mae_sd,energy_rdkit_mae_sd)
            mae_min,sd_min = stats_calc(energy_dft_mae_sd,energy_min_mae_sd)
            mae_dft,sd_dft = 0,0
        
        # setting some graphing options
        x_axis_names = [csearch_cmin_data['CSEARCH type']]

        if csearch_cmin_data['CMIN type'] is not None:
            x_axis_names.append(csearch_cmin_data['CMIN type'])
        if qm_data['Functional'] is not None:
            x_axis_names.append(qm_data['Functional']+'_'+qm_data['Basis set'])
        if qm_data['Functional SP'] is not None:
            x_axis_names.append(qm_data['Functional SP']+'_'+qm_data['Basis set SP'])

        name_to_write = '-'.join([str(int(elem)) for elem in x_axis_names])

        return list_all,mae_rdkit,sd_rdkit,mae_min,sd_min,mae_min,sd_min,mae_dft,sd_dft,x_axis_names,name_to_write


def get_qm_energy_plot(program,plot_files,energy_dft):
    type_csearch,type_cmin = None,None
    for file in plot_files:
        if program == 'g16':
            data_sp = cclib.io.ccread(file)
            energy_qm = data_sp.scfenergies[0]*ev_2_kcal_mol
        elif program == 'orca':
            sp_lines = open(file,"r").readlines()
            for i,_ in reversed(range(sp_lines)):
                if sp_lines[i].find('FINAL SINGLE POINT ENERGY') > -1:
                    energy_qm = float(sp_lines[i].split()[-1])*hartree_to_kcal
                    break
        if len(file.split('_ani.sdf')) == 2 or len(file.split('_xtb.sdf')) == 2:
            name = file.replace('_ani.sdf','_xtb.sdf').split('_xtb.sdf')[0]
            if len(file.split('_ani.sdf')) == 2:
                type_cmin = 'ANI'
            elif len(file.split('_xtb.sdf')) == 2:
                type_cmin = 'xTB'
        else:
            if len(file.split('_summ.sdf')) == 2:
                type_csearch = 'SUMM'
            elif len(file.split('_fullmonte.sdf')) == 2:
                type_csearch = 'Fullmonte'
            elif len(file.split('_rdkit.sdf')) == 2:
                type_csearch = 'RDKit'
            name = file.replace('_summ.sdf','_rdkit.sdf').replace('_fullmonte.sdf','_rdkit.sdf').split('_rdkit.sdf')[0]
        
        energy_dft.append([name,energy_qm])

    return energy_dft,type_csearch,type_cmin