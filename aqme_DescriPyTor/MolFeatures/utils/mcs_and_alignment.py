import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import openbabel
from networkx import relabel_nodes, Graph
from networkx.algorithms import isomorphism

from enum import Enum

from .file_handlers import *

class TABConstants(Enum):
    """
    Holds constants for tab constants
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic weights
    """
    covalent_radii = {
        'H': 0.31, 'He': 0.28, 'Li': 1.28,
        'Be': 0.96, 'B': 0.84, 'C': 0.76, 
        'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
        'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
        'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
        'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
        'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
        'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
        'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
        'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
        'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
        'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
        'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
        'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
        'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
        'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
        'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
        'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
        'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
        'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
        'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
        'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
        'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
        'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
        'Am': 1.80, 'Cm': 1.69
    }
    
    atomic_weights = {
        'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
        'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
        'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
        'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
        'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
        'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
        'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
        'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
        'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
        'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
        'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
        'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
        'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
        'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
        'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
        'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
        'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
        'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
        'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
        'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
        'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
        'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
        'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
        'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }

def extract_ordered_values(dictionary):
    items=dictionary.items()
    ordered_items=sorted(items, key=lambda x: x[0])
    ordered_values=[item[1] for item in ordered_items]
    return ordered_values

def get_structure_mapping_nx(G_1, G_2, single_mapping=True):
    result=False
    graph_matcher=isomorphism.GraphMatcher(G_1, G_2)
    if graph_matcher.is_isomorphic():
        if single_mapping:
            result=graph_matcher.mapping
        else:
            result=list(graph_matcher.isomorphisms_iter())
    return result

def get_sub_structure_mapping_nx(G_big, G_sub, single_mapping=True):
    result=False
    graph_matcher=isomorphism.GraphMatcher(G_big, G_sub)
    if graph_matcher.subgraph_is_isomorphic():
        if single_mapping:
            result=graph_matcher.mapping
        else:
            result=list(graph_matcher.subgraph_isomorphisms_iter())
    return result

def generate_conformers_rdkit(filename, num_conformers=5, output_filename='conformers.mol'):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    rdkit_mol=Chem.MolFromMolFile(filename)
    AllChem.EmbedMultipleConfs(rdkit_mol, numConfs=num_conformers)
    with open(output_filename, 'w') as outfile:
        for i in range(rdkit_mol.GetNumConformers()):
            conf=rdkit_mol.GetConformer(i)
            Chem.MolToMolFile(rdkit_mol, confId=i, fileName=output_filename, kekulize=False)

def generate_conformers_obabel(filename, num_conformers=5, output_filename='conformers.mol'):
# Load the molecule from the input file
    obConversion=openbabel.OBConversion()
    obConversion.SetInFormat('mol')
    ob_mol=openbabel.OBMol()
    obConversion.ReadFile(ob_mol, filename)

    ff=openbabel.OBForceField.FindForceField('MMFF94')
    ff.Setup(ob_mol)
    ff.DiverseConfGen(num_conformers, 0.5)
    ff.GetConformers(ob_mol)

    # Write the conformers to a new file
    obConversion.SetOutFormat('mol')
    obConversion.WriteFile(ob_mol, output_filename)

def get_low_level_conformers(molecule):
    pass

def rotate_molecule(molecule_1, molecule_2):
    pass

def get_xyz_df_from_file(xyz_filepath):
    """
    Write Me
    """
    xyz_df=pd.read_csv(xyz_filepath,
                       delim_whitespace=True,
                       skiprows=2,
                       names=['element', 'x', 'y', 'z'],
                       error_bad_lines=True) #'skip'
    return xyz_df

def expend_xyz_df(xyz_df, index_atoms=0):
    """
    xyz_df.columns=["element", "x", "y", "z"]
    expended_xyz_df:
       atom1_idx atom2_idx element       x      y       z  weight  cov_radius
    0         C0        C0       C  -1.944  0.247  -0.071  12.011        0.76
    1         C1        C1       C  -3.293 -0.055   0.126  12.011        0.76
    2         C2        C2       C  -4.288  0.537  -0.626  12.011        0.76        
    """
    expended_xyz_df=xyz_df.copy()
    expended_xyz_df['atom1_idx']=["{}{}".format(atm, idx) for atm, idx in zip(expended_xyz_df.element, expended_xyz_df.index.array)] #element + position in xyz C --> C0
    if index_atoms: #first atom has now index 1 in atom name, e.g. C --> C1
            expended_xyz_df['atom1_idx']=["{}{}".format(atm, idx+1) for atm, idx in zip(expended_xyz_df.element, expended_xyz_df.index.array)]
    expended_xyz_df['atom2_idx']=expended_xyz_df['atom1_idx'] #atom 1 is atom 2, duplicate for later use
    #atomic weight dict
    expended_xyz_df['weight']=expended_xyz_df['element'].apply(lambda x: TABConstants.atomic_weights.value[x]) # gemmi - alt: xyz_df['weight'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).weight)
    #covalent radius dict
    expended_xyz_df['cov_radius']=expended_xyz_df['element'].apply(lambda x: TABConstants.covalent_radii.value[x]) # gemmi - alt: xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).covalent_r)
    #reorder data frame
    expended_xyz_df=expended_xyz_df[['atom1_idx','atom2_idx','element','x','y','z','weight','cov_radius']]
    return expended_xyz_df

def dist_matrix_to_series1(dist_matrix): #dm_to_series1
    """
    removes the upper triangle of the distance matrix and zeros
    e.g., from d(A-B) = 1.234 Å = d(B-A) =1.234 Å, d(B-A) will be removed 
    d(A-B) = 0 Å will be removed as well
       C0  C1  C2           C0  C1  C2
    C0 0.0 1.1 2.3   ==\ C0 NaN NaN NaN
    C1 1.1 0.0 1.5   ==/ C1 1.1 NaN NaN
    C2 2.3 1.5 0.0       C2 2.3 1.5 NaN
    """
    dist_matrix=dist_matrix.astype(float) # do not comment this, angle list will be incomplete
    dist_matrix.values[np.triu_indices_from(dist_matrix, k=1)] = np.nan
    dist_matrix=dist_matrix.replace(0, np.nan) #replace zeros with nan
    return dist_matrix.unstack().dropna() #return and drop all nan

def get_distance_df(xyz_df, reduce_distance_matrix=True):
    """
    dist_mat_full:
       C0  C1  C2
    C0 0.0 1.1 2.3
    C1 1.1 0.0 1.5
    C2 2.3 1.5 0.0

    dist_df:
    
    """
    dist_mat_full=pd.DataFrame(squareform(pdist(xyz_df.iloc[:,3:6],'euclid')), #iloc [:,3:6] contains xyz coordinates
                               columns=xyz_df[['atom1_idx','element','cov_radius']],
                               index=xyz_df[['atom2_idx','element','cov_radius']])
    if reduce_distance_matrix:
            dist_mat_red=dist_matrix_to_series1(dist_mat_full) #remove the upper triangle and zeros
            dist_mat_red=dist_mat_red.reset_index(level=[1]) #bring it to a "normal" form
            dist_df=pd.DataFrame(dist_mat_red.to_records()) #bring it to a "normal" form, distance matrix --> disctance data frame
            normal_form_column='index'
            dist_rename_dict={'0':'distance_calc'}
    else:
            dist_mat = dist_mat_full.unstack()
            dist_mat = dist_mat.reset_index()
            dist_df=dist_mat
            normal_form_column='level_0'
            dist_rename_dict={0:'distance_calc'}
    #bring it to a "normal" form ...
    dist_df[['atom1_idx','element1','cov_radius1']]=pd.DataFrame(dist_df[normal_form_column].tolist(), index=dist_df.index)
    dist_df[['atom2_idx','element2','cov_radius2']]=pd.DataFrame(dist_df['level_1'].tolist(), index=dist_df.index)
    dist_df.drop([normal_form_column, 'level_1'], axis=1,inplace=True)
    dist_df.rename(columns=dist_rename_dict,inplace=True)
    dist_df=dist_df[['atom1_idx','element1','cov_radius1','atom2_idx','element2','cov_radius2','distance_calc']] #reorder data frame
    return dist_df, dist_mat_full

def determine_bonds(dist_df, determine_by='radii', radii_ext=0.08, ref_bonds_mask=None, include_bonds=None):
    """
    Write Me
    """
    new_dist_df=dist_df.copy()
    if determine_by=='radii':
            new_dist_df['distance_radii']=(new_dist_df['cov_radius1']+new_dist_df['cov_radius2'])*(1+radii_ext) #column with the sum of the atomic radii from two elements /atoms + x%
            new_dist_df['is_bond']=(new_dist_df['distance_calc']<new_dist_df['distance_radii']) #distance is considered as bond if the calculated distance is smaller than the sum of the atomic radii
    if determine_by=='ref':
            new_dist_df['is_bond']=ref_bonds_mask                
    #include distances from selected atom (pairs) from args include connections sets 'is_bond' to 'True' if 'False'
    if include_bonds:
            #must have atoms 1 and 2 and 'is_bond' must be 'False'
            new_dist_df.loc[(new_dist_df.atom1_idx.isin(include_bonds) & 
                                    new_dist_df.atom2_idx.isin(include_bonds) & 
                               ~new_dist_df.is_bond), 'is_bond'] = True
            #np variant of the above, slightly slower
    ##	dist_df['is_bond'] = np.where((dist_df.atom1_idx.isin(args.includeCon)) &
    ##	(dist_df.atom2_idx.isin(args.includeCon) & (~dist_df.is_bond)), True, dist_df['is_bond'])
    return new_dist_df

def expend_distance_df(dist_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    """
    Write Me
    """
    if bonds_by_ref:
            expended_dist_df=determine_bonds(dist_df, determine_by='ref', ref_bonds_mask=bonds_boolean_mask, include_bonds=include_bonds)
    else:
            expended_dist_df=determine_bonds(dist_df, determine_by='radii', radii_ext=radii_ext, include_bonds=include_bonds)
    return expended_dist_df

def get_expended_dist_df(xyz_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    expended_xyz_df=expend_xyz_df(xyz_df)
    dist_df, _=get_distance_df(expended_xyz_df)
    expended_dist_df=expend_distance_df(dist_df, radii_ext, include_bonds, bonds_by_ref, bonds_boolean_mask)
    return expended_dist_df

def connect_the_dots(xyz_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    expended_dist_df=get_expended_dist_df(xyz_df, radii_ext, include_bonds, bonds_by_ref, bonds_boolean_mask)
    sel_dist=pd.DataFrame(expended_dist_df[expended_dist_df.is_bond])
    return sel_dist

def atom_idx_to_index(atom_idx):
    """
    A function that takes a string of atom index like (C11 or Si13, etc.),
    and returns only the number from the index
    ----------
    Parameters
    ----------
    atom_idx : str.
    Atom index in the form of "element"+"index" like C11 or Si13
    -------
    Returns
    -------
    integer
    The number from the index 
    --------
    Examples
    --------
    atom_idx='C11'
    only_idx=atom_idx_to_index(atom_idx)
    print(only_symbol)
            11
    """
    return int(re.sub(r'\D+', '', atom_idx))

def get_bond_pairs(expended_dist_df):
    indices_1=expended_dist_df['atom1_idx'].apply(atom_idx_to_index).tolist() # remove ".apply(atom_idx_to_idx)"?
    indices_2=expended_dist_df['atom2_idx'].apply(atom_idx_to_index).tolist()
    bond_pairs=zip(indices_1, indices_2)
    return bond_pairs

def get_bond_pairs_from_xyz_df(xyz_df):
    dist_df=connect_the_dots(xyz_df.copy())
    bond_pairs=list(get_bond_pairs(dist_df))
    return bond_pairs

def xyz_df_to_nx_graph(xyz_df):
    atom_types=xyz_df['element'].to_numpy()
    bonded_pairs=get_bond_pairs_from_xyz_df(xyz_df)
    G=Graph()
    for atom_index, atom_type in enumerate(atom_types):
        G.add_node(atom_index,
                   element=atom_type)
    for bond in bonded_pairs:
        G.add_edge(bond[0], bond[1])
    return G

def align_xyz_df_by_mapping(xyz_df, mapping):
    actual_mapping=extract_ordered_values(mapping)            
    aligned_df=xyz_df.copy().reindex(actual_mapping).reset_index(drop=True)
    return aligned_df

def align_nx_graph_by_mapping(G, mapping):
    aligned_G=relabel_nodes(G, mapping)
    return aligned_G

def get_aligned_xyz_name(xyz_filepath):
    actual_name=xyz_filepath.split('.')[-2]
    aligned_name=actual_name+'_aligned'+FileExtensions.XYZ.value

def align_xyz_by_substructure(xyz_filepath, substructure_xyz, comment_line='comment_line'):
    # part 1 - renumbering
    output_filename=get_aligned_xyz_name(xyz_filepath)
    big_xyz_df=get_xyz_df_from_file(xyz_filepath)
    mcs_xyz_df=get_xyz_df_from_file(substructure_xyz)
    G_big=xyz_df_to_nx_graph(big_xyz_df)
    G_mcs=xyz_df_to_nx_graph(mcs_xyz_df)
    sub_mapping=get_sub_structure_mapping_nx(G_big, G_mcs)
    renumbered_big_xyz_df=align_xyz_df_by_mapping(big_xyz_df, sub_mapping)
    save_single_xyz_file(symbols=renumbered_big_xyz_df['element'].to_numpy(),
                         coordinates=renumbered_big_xyz_df[['x', 'y', 'z']].to_numpy(),
                         output_filename=output_filename,
                         comment_line=comment_line)
    # G_big_renumbered=align_nx_graph_by_mapping(G_big, sub_mapping)
    # part 2 - alignment via conformer search
    #get RMSD/align according to those common points
    #if no low RMSD match was found - run conformer search
    #build and renumber the sub_structure
    #re-construct the rest of the molecule
    #return the molecule
    return output_filename

def get_aligned_xyz_by_mcs_xyz(mcs_xyz, xyz_filepaths_to_align):
    aligned_xyz_filepaths=[align_xyz_by_substructure (xyz_filepath, mcs_xyz) for xyz_filepath in xyz_filepaths_to_align]
    return aligned_xyz_filepaths