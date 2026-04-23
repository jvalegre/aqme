import numpy as np
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from ..utils import help_functions
except:
    from utils import help_functions

import igraph as ig
import pandas as pd
def count_0(x):
    """Count the number of leading zeros in a pandas Series before the first non-zero entry."""
    return (x.cumsum() == 0).sum()

def get_transfomed_plane_for_sterimol(plane,degree):
    """
    a function that gets a plane and rotates it by a given degree
    in the case of sterimol the plane is the x,z plane.
    Parameters:
    ----------
    plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
    degree : float
    """
    # print(degree,plane)
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    transformed_plane=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(3)
    return transformed_plane

def calc_B1(transformed_plane,avs,edited_coordinates_df,column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    
    ## get the index of the min value in the column compared to the avs.min
    idx=np.where(np.isclose(np.abs(transformed_plane[:,column_index]),(avs.min()).round(4)))[0][0]
    if transformed_plane[idx,column_index]<0:
        new_idx=np.where(np.isclose(transformed_plane[:,column_index],transformed_plane[:,column_index].min()))[0][0]
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[new_idx,column_index],
                                 transformed_plane[:,column_index]<=transformed_plane[new_idx,column_index]+1)
        
        transformed_plane[:,column_index]=-transformed_plane[:,column_index]
    else:
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[idx,column_index]-1,
                                 transformed_plane[:,column_index]<=transformed_plane[idx,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]
    for i in range(1,transformed_plane.shape[0]): 
        if bool_list[i]:
            against.append(np.array(transformed_plane[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        if len(against)>0:
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
            
        else:
            B1.append(np.abs(transformed_plane[idx,column_index]+edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])
            
    # print(f'B1: {B1}, B1_loc: {B1_loc}')      
    return [B1,B1_loc]

def b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane):
    """
    a function that gets a plane transform it and calculate the b1s for each degree.
    checks if the plane is in the x or z axis and calculates the b1s accordingly.
    Parameters:
    ----------
    extended_df : pd.DataFrame
    b1s : list
    b1s_loc : list
    degree_list : list
    plane : np.array
    """
    degree=[]
    for degree in degree_list:
        transformed_plane=get_transfomed_plane_for_sterimol(plane, degree)
        
        avs=np.abs([max(transformed_plane[:,0]),min(transformed_plane[:,0]), 
                    max(transformed_plane[:,1]),min(transformed_plane[:,1])])
        
        if min(avs) == 0:
            min_avs_indices = np.where(avs == min(avs))[0]
            if any(index in [0, 1] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                continue  # Skip the rest of the loop

            elif any(index in [2, 3] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                continue

        if np.where(avs==avs.min())[0][0] in [0,1]:
            B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,0)
            
  
        elif np.where(avs==avs.min())[0][0] in [2,3]:
            B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,1)
             
        
        b1s.append(np.unique(np.vstack(B1)).max())####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())

def get_b1s_list(extended_df, scans=90//5):
    
    b1s,b1s_loc=[],[]
    scans=scans
    degree_list=list(range(18,108,scans))
    plane=np.array(extended_df[['x','z']].astype(float))
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    
    if b1s:
        try:
            back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
            front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
            degree_list=range(back_ang,front_ang+1)
        except:
            
            back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
            front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
            degree_list=range(back_ang,front_ang+1)
    else:
        print('no b1s found')
        return [np.array(b1s),np.array(b1s_loc)]
    # print(f'specific degree list: {degree_list}')
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    # print(f'b1 arrays: {[np.array(b1s),np.array(b1s_loc)]}')
    return [np.array(b1s),np.array(b1s_loc)]

def get_molecule_connections(bonds_df,source,direction):
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    longest_path=np.unique(help_functions.flatten_list(with_direction))
    return longest_path


from scipy.special import cbrt

def calc_angle(p1, p2, degrees: bool=False) -> float: ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_p1=np.linalg.norm(p1)
    norm_p2=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_p1*norm_p2))
    if degrees:
        thetha=np.degrees(thetha)   
    return thetha

def vec_organizer(block, density, blocks, x_steps, y_steps, z_steps):
    
    start_index = blocks[block] 
    end_index = start_index + int(np.ceil(z_steps / 6)) -1
    vec = density.iloc[start_index:end_index + 1].values.flatten()
    # Calculate the required vector length
    target_length = int(np.ceil(cbrt(x_steps * y_steps * z_steps)))-1
    # Append zeros if the vector is shorter than the required length
    if len(vec) < target_length:
        vec = np.append(vec, np.zeros(target_length - len(vec)))
    return vec

def pt_space_block(non_zero,x_blocks):

    df = x_blocks[non_zero]
    # Get indices of non-zero elements, np.argwhere returns a 2D array of indices
    non_zero_indices = np.argwhere(df.to_numpy() != 0)
    # Create a DataFrame from non-zero indices
    non_zero_df = pd.DataFrame(non_zero_indices, columns=['Row', 'Column'])
    ## rearrage the columns by the order of Column
    
    # Concatenate original block with its non-zero indices
    non_zero=pd.DataFrame([non_zero]*len(non_zero_df))
    return pd.concat([non_zero, non_zero_df], axis=1).sort_values(by='Column')

def pt_space_block_binder(non_zero,x_blocks):
    # Apply pt_space_block to each DataFrame in the list and concatenate them
    
    results = [pt_space_block(non_zero[i], x_blocks) for i in range(len(non_zero))]
    return pd.concat(results, ignore_index=True)


def dens_to_pt(point, dense_points, x_origin, y_origin, z_origin, x_size, y_size, z_size):
    # Extract the row corresponding to 'point'
    num_pt = (dense_points.iloc[point].values)+1
    
    
    # Calculate new coordinates
    x_coord = x_origin + num_pt[0] * x_size
    y_coord = y_origin + num_pt[1] * y_size
    z_coord = z_origin + num_pt[2] * z_size
    
    return x_coord, y_coord, z_coord

from scipy.spatial.distance import pdist, squareform

def extract_connectivity(xyz_df, threshhold_distance=1.82):
    coordinates=np.array(xyz_df[['x','y','z']].values)
    atoms_symbol=np.array(xyz_df['atom'].values)
    # compute the pairwise distances between the points
    distances = pdist(coordinates)
    # convert the flat array of distances into a distance matrix
    dist_matrix = squareform(distances)
    dist_df=pd.DataFrame(dist_matrix).stack().reset_index()
    dist_df.columns = ['a1', 'a2', 'value']
    dist_df['first_atom']=[atoms_symbol[i] for i in dist_df['a1']]
    dist_df['second_atom']=[atoms_symbol[i] for i in dist_df['a2']]
    remove_list=[]
    dist_array=np.array(dist_df)
    remove_list = []
    for idx, row in enumerate(dist_array):
        remove_flag = False
      
        if row[0] == row[1]:
            remove_flag = True
          
        if ((row[3] == 'H') & (row[4] not in help_functions.XYZConstants.NOF_ATOMS.value)):
            remove_flag = True
           
        if ((row[3] == 'H') & (row[4] == 'H')):
            remove_flag = True
          
        if (((row[3] == 'H') | (row[4] == 'H')) & (row[2] >= 1.5)):
            remove_flag = True
           
        if ((row[2] >= threshhold_distance) | (row[2] == 0)):
            remove_flag = True
      

        if remove_flag:
            remove_list.append(idx)

    dist_df=dist_df.drop(remove_list)
    dist_df[['min_col', 'max_col']] = pd.DataFrame(np.sort(dist_df[['a1', 'a2']], axis=1), index=dist_df.index)
    dist_df = dist_df.drop(columns=['a1', 'a2']).rename(columns={'min_col': 0, 'max_col': 1})
    dist_df = dist_df.drop_duplicates(subset=[0, 1])
    return pd.DataFrame(dist_df[[0,1]]+1)

def direction_atoms_for_sterimol(bonds_df,base_atoms)->list: #help function for sterinol
    """
    a function that return the base atom indices for coordination transformation according to the bonded atoms.
    you can insert two atom indicess-[1,2] output [1,2,8] or the second bonded atom
    if the first one repeats-[1,2,1] output [1,2,3]
    """
    
    base_atoms_copy=base_atoms[0:2]
    origin,direction=base_atoms[0],base_atoms[1]
    bonds_df = bonds_df[~((bonds_df[0] == origin) & (bonds_df[1] == direction)) & 
                              ~((bonds_df[0] == direction) & (bonds_df[1] == origin))]
    
    try :
        base_atoms[2]==origin
        if(any(bonds_df[0]==direction)):
            # take the second atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[0]==direction)][1].iloc[1])
        else:
            # take the first atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[1]==direction)][0].iloc[1])
    except: 
        
        for _, row in bonds_df.iterrows():
            if row[0] == direction:
                base_atoms_copy.append(row[1])
                break
            elif row[1] == direction:
                base_atoms_copy.append(row[0])
                break
    return base_atoms_copy



def get_molecule_connections(bonds_df,source,direction):
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    longest_path=np.unique(help_functions.flatten_list(with_direction))
    return longest_path

def calc_new_base_atoms(coordinates_array, atom_indices):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices is 4, the origin will be the middle of the first two atoms.
    """
    new_origin=coordinates_array[atom_indices[0], 0:3]
    
    if (len(atom_indices)==4):
        new_origin=(new_origin+coordinates_array[atom_indices[1]])/2
    new_y=(coordinates_array[atom_indices[-2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-2]]-new_origin))
    coplane=((coordinates_array[atom_indices[-1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-1]]-new_origin)+0.00000001))
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results


def transform_row(row_array, new_basis, new_origin, round_digits):
    translocated_row = row_array - new_origin
    return np.dot(new_basis, translocated_row).round(round_digits)

def mag(x):
    return np.linalg.norm(x, axis=1)

def mag_2d(x):
    return np.sqrt(x[0] ** 2 + x[2] ** 2)
    

def close_atom(dens_pt, trans_co, n_atoms):
    # Calculate vector differences
    differences = trans_co[dens_pt, :3] - trans_co[:n_atoms, :3]
    # Apply the magnitude function
    distances = mag(differences)
    # Return the index of the minimum distance (1-based index)
    return np.argmin(distances) + 1

ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}


class cube():

    def __init__(self, fname=None, base_atoms=[1,2]):
        self.fname=fname
        self.base_atoms=base_atoms
        self.xyz_from_cube()
        self.process_cube_data()
        self.get_sterimol()
        
    def xyz_from_cube(self):
        # Remove the existing .xyz file if it exists
        xyz_file_path = os.path.splitext(self.fname)[0] + ".xyz"
        if os.path.exists(xyz_file_path):
            os.unlink(xyz_file_path)

        # Read data from the cube file
        with open(self.fname, 'r') as file:
            self.cube_data = file.readlines()

        # Processing the header and extracting the number of atoms
        self.n_atoms = int(self.cube_data[2].split()[0])
        atom_data = self.cube_data[6:6 + self.n_atoms]

        # Creating a DataFrame from the atom data
        
        atoms = []
        for line in atom_data:
            parts = line.split()
            atom_type = int(parts[0])  # Assuming the atomic number is given
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            # Convert from Bohr to Angstrom
            atoms.append([atom_type, x * 0.529177249, y * 0.529177249, z * 0.529177249])
            # atoms_reg.append([atom_type, x , y , z ])

        # Create a DataFrame
        df = pd.DataFrame(atoms, columns=['atom', 'x', 'y', 'z'])
        atoms=df['atom'].astype(str).map(ATOMIC_NUMBERS)
        df['atom']=atoms
        self.xyz_df=df
        return None
    
    def process_cube_data_to_array(self):
        # Determine the maximum length of rows
        max_length = max(len(row) for row in self.cube_data)

        # Pad rows with NaN to make them all the same length
        padded_data = []
        for row in self.cube_data[2:]:
            float_row = [float(x) for x in row.split()]
            # if len(float_row) < max_length:
            #     float_row.extend([np.nan] * (max_length - len(float_row)))
            padded_data.append(float_row)

        # Convert the padded data to a numpy array
        array_data = np.array(padded_data)
        return array_data

    def pad_cube_data(self,data, max_length):
        processed_data = []
        for line in data:
            float_values = list(map(float, line.split()))
            # Pad the row with NaN to ensure it has exactly 6 columns
            while len(float_values) < max_length:
                float_values.append(np.nan)
            processed_data.append(float_values)
        return pd.DataFrame(processed_data)
    
    def process_cube_data(self,isovalue=0.003):

        self.structure = [self.cube_data[2:][i] for i in range(4 + self.n_atoms)]
        self.density=self.cube_data[(6 + self.n_atoms):]
        self.structure_df=self.pad_cube_data(self.structure, 5)
        self.density_df=self.pad_cube_data(self.density, 6)
        x_origin, y_origin, z_origin = self.structure_df.iloc[0,1], self.structure_df.iloc[0,2], self.structure_df.iloc[0,3]
        x_step, y_step, z_step = self.structure_df.iloc[1,0], self.structure_df.iloc[2,0], self.structure_df.iloc[3,0]
        x_size, y_size, z_size = self.structure_df.iloc[1,1], self.structure_df.iloc[2,2], self.structure_df.iloc[3,3]

        self.bonds_df=extract_connectivity(self.xyz_df)
        rows_per_block = int(np.ceil(z_step / 6)) + 1 # 12
        target_length = int(np.ceil(cbrt(x_step * y_step * z_step)))-1 # 80
        blocks = list(range(0, len(self.density_df) , rows_per_block-1))
        self.density_df=self.density_df.fillna(0)
        row_dens_data = [vec_organizer(block, self.density_df, blocks, x_step, y_step, z_step) for block in range(len(blocks))]

        row_dens_df=pd.DataFrame()
        for i in range(target_length):
            row=[]
            for j in range(len(blocks)):
                row.append(row_dens_data[j][i])
            row_dens_df[i]=row
        row_dens_df[(row_dens_df > isovalue * 1.1) | (row_dens_df < isovalue)] = 0
        
        x_blocks = []
        for i in range(0, len(row_dens_df), int(y_step)):
            x_blocks.append(row_dens_df.iloc[i:i+int(y_step)])
        # x_blocks_df = pd.concat(x_blocks, axis=1)
        non_zero_regions = [index for index, df in enumerate(x_blocks) if df.sum().sum() != 0]
   
        dense_points=pt_space_block_binder(non_zero_regions, x_blocks)
        x=[]
        y=[]
        z=[]
        for point in range(len(dense_points)):
            x_coord, y_coord, z_coord = dens_to_pt(point, dense_points, x_origin, y_origin, z_origin, x_size, y_size, z_size)
            x.append(x_coord)
            y.append(y_coord)
            z.append(z_coord)
        coordinates_space_df=pd.DataFrame({'x':x,'y':y,'z':z})
        df=self.xyz_df[['x','y','z']]/0.529177249
        self.xyz=pd.concat([df,coordinates_space_df]).reset_index(drop=True)

    def get_sterimol(self):

        if isinstance(self.base_atoms[0], list):
            sterimol_list=[]
            for base_atom in self.base_atoms:
                sterimol_list.append(self.calc_sterimol(base_atom))
            self.sterimol_df=pd.concat(sterimol_list, axis=1)
            print(self.sterimol_df)
        else:
            self.sterimol_df=self.calc_sterimol(self.base_atoms)

    def calc_sterimol(self,base_atoms_indices):      
       
        base_atoms=np.array(direction_atoms_for_sterimol(self.bonds_df,base_atoms_indices))-1
       
        tag=None
        
        if count_0(self.xyz.iloc[:self.n_atoms, :3].sum()) >= 2:
            
            col_1 = count_0(self.xyz.iloc[:, 0])
            col_2 = count_0(self.xyz.iloc[:, 1])
            col_3 = count_0(self.xyz.iloc[:, 2])
            # Find the column with the maximum count of leading zeros
            place_num = np.argmax([col_1, col_2, col_3]) + 1  # +1 because Python uses zero-based indexing
            # Create a new row where the max place is 1 and others are 0
            new_row = pd.DataFrame([[0, 0, 0]], columns=self.xyz.columns)
            new_row.iloc[0, place_num - 1] = 1  # -1 to adjust for zero-based index
            # Insert the new row at position n_atoms
            self.xyz = pd.concat([self.xyz.iloc[:self.n_atoms], new_row, self.xyz.iloc[self.n_atoms:]]).reset_index(drop=True)
            # Suppose numeric_atoms is a list that holds some values where the 3rd element is to be updated
            numeric_atoms = [None, None, None]  # Example initialization
            numeric_atoms[2] = self.n_atoms + 1  # Update the third element
            tag = 1  # Set tag to 1 as per the R script logic
 
        new_origin,new_y,coplane=calc_new_base_atoms(self.xyz.values,base_atoms)
        cross_y_coplane=np.cross(coplane,new_y)
        coef_mat=np.stack([new_y,coplane,cross_y_coplane])
        angle_new_y_coplane=calc_angle(coplane,new_y)
        cop_ang_x=angle_new_y_coplane-(np.pi/2)
        result_vector=[0,np.cos(cop_ang_x),0]
        #result_vector=[np.cos(cop_ang_x), 0, 0]
        new_x,_,_,_=np.linalg.lstsq(coef_mat,result_vector,rcond=None)
        new_basis=np_cross_and_vstack(new_x, new_y)
        coordinates_array=self.xyz.values
        transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, 6), 1,
                                                  coordinates_array)
        nan_atoms=[np.nan]*self.n_atoms
        close_atoms=[]
        for i in range(self.n_atoms , len(transformed_coordinates)):
            close_atoms.append(close_atom(i, transformed_coordinates, self.n_atoms))
        ## add a column to the transformed coordinates of close atoms, put Nan in length of the original atoms first
        close_atoms=nan_atoms+close_atoms
        ## add the close atoms to the transformed coordinates
        transformed_coordinates=np.column_stack([transformed_coordinates,close_atoms])
        rlev=get_molecule_connections(self.bonds_df,base_atoms_indices[0],base_atoms_indices[1])
        
        filtered_indices = np.isin(transformed_coordinates[:, 3], rlev)
        tcs = transformed_coordinates[filtered_indices, :3]
        ## now make rest_tcs from the rest of the transformed coordinates
        rest_tcs = transformed_coordinates[~filtered_indices, :3]
        
        mutate=[mag_2d(tcs[i])*0.529177249 for i in range(len(tcs))]
        tcs=np.column_stack([tcs,mutate])
        L=max(tcs[:,1]) * 0.529177249
        B5=max(tcs[:,3])
        # make plane variable form column 0 and 2
        plane=tcs[:,[0,2]]
        df=pd.DataFrame(columns=['b1', 'b1_loc', 'x_b1', 'z_b1'])
        for i in range(1,91):
            tc_plane=get_transfomed_plane_for_sterimol(plane,i)
            avs=np.abs([max(tc_plane[:,0]),min(tc_plane[:,0]), 
                    max(tc_plane[:,1]),min(tc_plane[:,1])])
            
            tc_plane=tc_plane.round(3)
            
            if np.where(avs==avs.min())[0][0] in [0,1]: 

                idx=np.where(np.isclose(np.abs(tc_plane[:,0]),(avs.min()).round(3)))[0][0]
                value_from_tcs = tcs[idx, 1]
                B1_loc=  value_from_tcs * 0.529177249
                x_b1 = tcs[idx, 0]
                z_b1 = tcs[idx, 2]
            elif np.where(avs==avs.min())[0][0] in [2,3]:

                idx=np.where(np.isclose(np.abs(tc_plane[:,1]),(avs.min()).round(3)))[0][0]
                value_from_tcs = tcs[idx, 1]
                B1_loc=  value_from_tcs * 0.529177249
                x_b1 = tcs[idx, 0]
                z_b1 = tcs[idx, 2]
            b1=np.abs(min(avs))* 0.529177249
            b_df=pd.DataFrame({'b1': b1, 'b1_loc': B1_loc, 'x_b1': x_b1, 'z_b1': z_b1}, index=[0])
            df=pd.concat([df,b_df])
        
        df=df.reset_index(drop=True)
        
        B1=df['b1'].min()
        x_b1=df['x_b1'].iloc[df['b1'].idxmin()]
        z_b1=df['z_b1'].iloc[df['b1'].idxmin()]
        loc_b1=df['b1_loc'].iloc[df['b1'].idxmin()]
        loc_b5=tcs[np.argmax(tcs[:,3]),1]*0.529177249

        idx_b5 = np.argmax(tcs[:, 3])
        x_b5 = tcs[idx_b5, 0]
        z_b5 = tcs[idx_b5, 2]
        v1 = np.array([x_b1, z_b1])
        v2 = np.array([x_b5, z_b5])
        
        # Compute the angle (in radians) and convert to degrees
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        angle_rad = np.arccos(dot_product / (norm_v1 * norm_v2))
        angle_deg = np.degrees(angle_rad)

        name=self.fname.split('.')[0]
        sterimol_df = pd.DataFrame({
            f'B1_{base_atoms_indices}': [B1],
            f'B5_{base_atoms_indices}': [B5],
            f'L_{base_atoms_indices}': [L],
            f'loc_b5_{base_atoms_indices}': [loc_b5],
            f'B1_B5_angle_{base_atoms_indices}': [angle_deg]
        }, index=[name])
        
       
        return sterimol_df

class cube_many():

    def __init__(self, path, base_atoms):
        self.path=path
        os.chdir(path)
        self.file_names=[file for file in os.listdir(path) if file.endswith('.cube')]
        self.base_atoms=base_atoms
        self.get_sterimol_many()

    def get_sterimol_many(self):
        sterimol_list = []
        
        for file in self.file_names:
            # try:
                cube_file = cube(file, self.base_atoms)
                df = cube_file.sterimol_df.copy()
                df.index = [file]  # Set the index to the file name
                sterimol_list.append(df)
            # except Exception as e:
            #     print(f'File {file} failed: {e}')

        # Concatenate all DataFrames, indexed by file name
        self.sterimol_df = pd.concat(sterimol_list, axis=0)
        print(self.sterimol_df)
        # save as csv
        self.sterimol_df.to_csv('cube_sterimol')
        return self.sterimol_df  # Return for debugging if needed

if __name__ == "__main__":
    pass