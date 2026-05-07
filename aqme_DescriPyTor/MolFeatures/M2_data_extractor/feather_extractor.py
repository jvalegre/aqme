import re
import pandas as pd
import numpy as np
import os
from enum import Enum
from typing import List, Optional
import re
import pandas as pd
from typing import List, Pattern, Union
class ReExpressions(Enum):
    FLOAT = r'[-+]?[0-9]+\.[0-9]+'
    FLOAT_POL= r'-?\d*\.\d*'
    FLOATS_ONLY= "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
    BONDS= r'R\(\d+.\d+'
    FREQUENCY= r'\-{19}'
    CHARGES=r'^-?[0-9]+([.][0-9]+)?$'

class FileFlags(Enum):
    DIPOLE_START='Dipole moment'
    DIPOLE_END='Quadrupole moment'
    MAIN_CUTOFF= r'\-{69}'
    STANDARD_ORIENTATION_START='Standard orientation:'
    POL_START='iso'
    POL_END='xx'
    CHARGE_START='Summary of'
    CHARGE_END='====='
    FREQUENCY_START='Harmonic frequencies'
    FREQUENCY_END='Thermochemistry'
    HIRSH_CHARGE_START='Hirshfeld'
    HIRSH_CHARGE_END='Hirshfeld charges with hydrogens '

class Names(Enum):
    
    DIPOLE_COLUMNS=['dip_x','dip_y','dip_z','total_dipole']
    STANDARD_ORIENTATION_COLUMNS=['atom','x','y','z']
    DF_LIST=['standard_orientation_df', 'dipole_df', 'pol_df','energy_value' ,'atype_df','charge_df', 'bonds_df', 'info_df']


class GeneralConstants(Enum):    
    ATOMIC_NUMBERS ={
     '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
              '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
     

    


def search_phrase_in_text(text_lines: str, key_phrase: str) : 
        """
        Searches for a key phrase in a list of text lines.
        
        Parameters:
            text_lines: text string to search through.
            key_phrase: The phrase to search for.
            
        Returns:
            The index of the first line where the key phrase is found, or None if not found.
        """   
        search_result=re.compile(key_phrase).search(text_lines)
        # print(f"Found key phrase '{key_phrase}' at line {search_result.start()}") if search_result else print(f"Key phrase '{key_phrase}' not found")
        return search_result

def extract_lines_from_text(
    text: str,
    re_expression: Union[str, Pattern[str]]
) -> List[str]:
    """
    Extract all substrings from `text` that match the given regular expression.

    Parameters:
        text (str): The full text to search.
        re_expression (str | Pattern[str]): A regex pattern (as a string or compiled Pattern).

    Returns:
        List[str]: A list of all non-overlapping matches.
    """
    pattern = re_expression if isinstance(re_expression, re.Pattern) else re.compile(re_expression)
    return pattern.findall(text)


def find_all_matches(
    text: str,
    key_phrase: Union[str, Pattern[str]]
) -> List[re.Match[str]]:
    """
    Find all occurrences of `key_phrase` in `text`, returning the Match objects.

    Parameters:
        text (str): The text to search.
        key_phrase (str | Pattern[str]): A substring or regex pattern to locate.

    Returns:
        List[re.Match[str]]: A list of Match objects, in order of appearance.
    """
    pattern = key_phrase if isinstance(key_phrase, re.Pattern) else re.compile(key_phrase)
    return list(pattern.finditer(text))


def process_gaussian_charge_text(log_text: str) -> pd.DataFrame:
    """
    Parse a Gaussian log file’s text to extract Mulliken charges and return them in a DataFrame.

    This function locates the last occurrence of the start and end markers
    defined by `FileFlags.CHARGE_START.value` and `FileFlags.CHARGE_END.value`,
    slices out that block, extracts all floating-point numbers, then picks
    every 6th value starting from the second one (index 1) to form the charge list.

    Parameters:
        log_text (str): Entire contents of the Gaussian log file.

    Returns:
        pandas.DataFrame: A single-column DataFrame with header 'charge'.

    Raises:
        ValueError: If start/end markers aren’t found or no floats are extracted.
    """
    # Locate charge block
    starts = find_all_matches(log_text, FileFlags.CHARGE_START.value)
    ends   = find_all_matches(log_text, FileFlags.CHARGE_END.value)
    if not starts or not ends:
        raise ValueError("Could not locate charge start/end markers in log text.")

    start_idx = starts[-1].end()
    end_idx   = ends[-1].start()
    section   = log_text[start_idx:end_idx]

    # Extract all floats in that block
    floats = extract_lines_from_text(section, ReExpressions.FLOATS_ONLY.value)
    if not floats:
        raise ValueError("No floating‐point numbers found in charge section.")

    # Convert to float and pick every 6th entry starting with index 1
    arr = np.array(floats, dtype=float)
    charge_values = arr[1::6]

    return pd.DataFrame({"nbo_charge": charge_values})


def process_gaussian_dipole_text(log_text: str) -> pd.DataFrame:
    """
    Extract the dipole moment components from a Gaussian log file.

    The function finds the last block between the DIPOLE_START and DIPOLE_END markers,
    extracts exactly four floating-point numbers, and returns them in a one-row DataFrame
    using the column names defined in Names.DIPOLE_COLUMNS.value.

    Parameters:
        log_text (str): Full text of the Gaussian log file.

    Returns:
        pd.DataFrame: A DataFrame with a single row and columns = Names.DIPOLE_COLUMNS.value.

    Raises:
        ValueError: If start/end markers are not found or fewer than four floats are extracted.
    """
    # Find all occurrences of the start and end markers
    starts = find_all_matches(log_text, FileFlags.DIPOLE_START.value)
    ends   = find_all_matches(log_text, FileFlags.DIPOLE_END.value)
    if not starts or not ends:
        raise ValueError("Could not locate dipole start/end markers in log text.")

    # Use the last occurrence
    start_idx = starts[-1].end()
    end_idx   = ends[-1].start()
    block     = log_text[start_idx:end_idx]

    # Extract floats
    floats = extract_lines_from_text(block, ReExpressions.FLOAT.value)
    if len(floats) < 4:
        raise ValueError(f"Expected 4 dipole components, found {len(floats)}.")

    # Take exactly the first 4 values
    values = floats[:4]
    # Build DataFrame
    df = pd.DataFrame([values], columns=Names.DIPOLE_COLUMNS.value)
    return df


def gauss_first_split(log_text: str) -> List[str]:
    """
    Perform an initial two-stage split of Gaussian log text:
      1. Split by the STANDARD_ORIENTATION_START marker and take the last chunk.
      2. Split that chunk by the MAIN_CUTOFF marker.

    Parameters:
        log_text (str): Full Gaussian log file text.

    Returns:
        List[str]: A list of text sections after the two splits.
    """
    # Split off everything before the final standard orientation block
    after_std = re.split(FileFlags.STANDARD_ORIENTATION_START.value, log_text)[-1]
    # Then split that by the main cutoff
    parts = re.split(FileFlags.MAIN_CUTOFF.value, after_std)
    return parts


def process_gaussian_standard_orientation_text(log_text: str) -> pd.DataFrame:
    """
    Extract the Standard Orientation coordinates table from a Gaussian log file.

    This function:
     1. Finds all floating-point numbers in the text.
     2. Reshapes them into rows of six values.
     3. Removes the first two columns (atomic index and atomic number).
     4. Builds a DataFrame with columns = Names.STANDARD_ORIENTATION_COLUMNS.value.
     5. Converts the 'atom' column from atomic number to symbol via GeneralConstants.ATOMIC_NUMBERS.value.
     6. Ensures x, y, z columns are floats.
     7. Prepends two header rows: a NaN row, and a row whose 'atom' cell holds the number of atoms.

    Parameters:
        log_text (str): Full Gaussian log file text.

    Returns:
        pd.DataFrame: A DataFrame representing the standard orientation table,
                      with two extra header rows at the top.
    """
    # 1) Extract all floats
    float_lines = extract_lines_from_text(log_text, ReExpressions.FLOATS_ONLY.value)
    arr = np.array(float_lines, dtype=float)

    # 2) Reshape into N×6 matrix
    if arr.size % 6 != 0:
        raise ValueError("Floating-point count is not a multiple of 6.")
    mat = arr.reshape(-1, 6)

    # 3) Drop columns 0 and 2 (index and atomic number)
    coords = np.delete(mat, (0, 2), axis=1)

    # 4) Build DataFrame
    df = pd.DataFrame(coords, columns=Names.STANDARD_ORIENTATION_COLUMNS.value)

    # 5) Convert atomic numbers to symbols
    
    df['atom'] = df['atom'].astype(int).astype(str)  # Ensure 'atom' is string type for mapping
    mapping = GeneralConstants.ATOMIC_NUMBERS.value  # e.g. {1:'H', 6:'C', 7:'N', 8:'O', ...}
    df['atom'] = df['atom'].map(mapping)
    
    # 6) Ensure x,y,z are floats
    for col in ['x', 'y', 'z']:
        df[col] = df[col].astype(float)

    # 7) Prepend two header rows
    n_atoms = df.shape[0]
    headers = pd.DataFrame(
        [[np.nan] * df.shape[1],  # first header row of NaNs
         [n_atoms] + [np.nan] * (df.shape[1] - 1)],  # second row: atom-count in 'atom' column
        columns=df.columns
    )
    final_df = pd.concat([headers, df], ignore_index=True)
    
    return final_df


def search_phrase_in_text_pol(text: str, key_phrase: str) -> int:
    """
    Search for a key phrase in a text string and return the start position of the match.

    Parameters:
        text (str): The text to search within.
        key_phrase (str): The phrase to search for.

    Returns:
        int: The start index of the first match or None if not found.
    """
    pattern = re.compile(rf"{re.escape(str(key_phrase))}\s")
    match = pattern.search(text)
    if match:
        return match.start()
    return None

def process_hirshfeld_charges(log_file_lines):
   
    charges_start = find_all_matches(log_file_lines, key_phrase=FileFlags.HIRSH_CHARGE_START.value)
    charges_end = find_all_matches(log_file_lines, key_phrase=FileFlags.HIRSH_CHARGE_END.value)

    if charges_start and charges_end:
        charges_start = charges_start[-2].end()
        charges_end = charges_end[-1].start()
        
        lines=log_file_lines[charges_start:charges_end]
    
    selected_lines = (extract_lines_from_text(lines,
                                             re_expression=ReExpressions.FLOAT.value))
   
    charge_array = np.array(selected_lines).reshape(-1,6)
    
  
    hirsh_charge_array=charge_array[:,0]
    cm5_charge_array=charge_array[:,5]

    

    return pd.DataFrame(hirsh_charge_array, columns=['hirshfeld_charge']),pd.DataFrame(cm5_charge_array, columns=['cm5_charge'])



def process_gaussian_pol_text(log_file_lines: List[str]) -> Optional[pd.DataFrame]:
    """
    Processes Gaussian polarization data and returns a DataFrame.
    
    Parameters:
        log_file_lines: List of lines from the Gaussian log file.
        
    Returns:
        A DataFrame containing polarization data or None if not found.
    """
    # If log_file_lines is a list, convert it to a single string
    if isinstance(log_file_lines, list):
        log_file_text = '\n'.join(log_file_lines)
    else:
        log_file_text = log_file_lines

    pol_start = search_phrase_in_text_pol(log_file_text, key_phrase=FileFlags.POL_START.value)
    pol_end = search_phrase_in_text_pol(log_file_text, key_phrase=FileFlags.POL_END.value)
    
    if pol_start is not None and pol_end is not None:
        pol = extract_lines_from_text(log_file_lines[pol_start:pol_end], re_expression=ReExpressions.FLOAT_POL.value)
        pol_df = pd.DataFrame([float(pol[0])*1000, float(pol[5])*1000], index=['iso', 'aniso'], dtype=float).T
        return pol_df
    else:
        # print("Failed to create.")
        pol_df = pd.DataFrame([100, 100], index=['iso', 'aniso'], dtype=float).T
        return pol_df




def process_gaussian_bonds(log_file_lines):
    bonds=extract_lines_from_text(log_file_lines, re_expression=ReExpressions.BONDS.value)
    bonds_text=[re.sub(r'R\(','',line).split(',') for line in bonds]
    return bonds_text


def remove_floats_until_first_int(input_list):
    output_list = []
    encountered_integer = False
    for item in input_list:
        try:
            # Try converting the string to an integer
            int_item = int(item)
            encountered_integer = True
            output_list.append(item)
        except ValueError:
            # If conversion to integer fails, try float
            try:
                float_item = float(item)
                if encountered_integer:
                    output_list.append(item)
            except ValueError:
                # If conversion to both integer and float fails, keep the item
                output_list.append(item)
    return output_list


def process_gaussian_vibs_string(log_file_lines):
    pattern = re.compile(rf'{FileFlags.FREQUENCY_START.value}([\s\S]*?){FileFlags.FREQUENCY_END.value}')
    match = pattern.search(log_file_lines)
    if match:
        vibration_section = match.group(1).strip()

        # Further slicing the obtained text based on 'Frequencies --'
        frequencies_blocks = re.split(r'Frequencies --', vibration_section)

        # Removing the last part that contains '------'
        final_blocks = [block.split('-------------------')[0].strip() for block in frequencies_blocks[1:]]

        return final_blocks
    
def process_gaussian_info(frequency_string):
    ir,frequency=[],[]
    for data in frequency_string:
        match=extract_lines_from_text(data, re_expression=ReExpressions.FLOATS_ONLY.value)
        frequency.append((match[0:3]))
        ir.append(match[6:9])
    info_df=pd.DataFrame()
    info_df['Frequency']=[(item) for sublist in frequency for item in sublist]
    info_df['IR']=[(item) for sublist in ir for item in sublist] 
    info_df=info_df.astype(float)

    return info_df

def vib_array_list_to_df(array_list):
    array_list_df = []
    col_names = []

    for i, array in enumerate(array_list, 1):
        new_array = np.delete(array, [0, 1], axis=1).reshape(-1, 3)
        new_df = pd.DataFrame(new_array)
        # Assign mode-specific x, y, z columns
        new_df.columns = [f'mode_{i}_x', f'mode_{i}_y', f'mode_{i}_z']
        array_list_df.append(new_df)

    # Concatenate horizontally (each mode's x/y/z in its own block)
    vibs_df = pd.concat(array_list_df, axis=1)
    vibs_df = vibs_df.astype(float)
    return vibs_df



def process_gaussian_frequency_string(final_blocks):
    """
    Processes Gaussian frequency output blocks into a DataFrame of vibrational data.
    """
    vibs_list = []
    lengths = []
    for i, data in enumerate(final_blocks):
        # Extract all floats
        match = re.findall(ReExpressions.FLOATS_ONLY.value, data)
        # Remove the first 12 floats (headers or irrelevant values)
        match = match[12:]
        match = np.array(match)
        
        # Special cleaning for the second block if needed
        if i == 1:
            match = remove_floats_until_first_int(match)
        
        lengths.append(len(match))
        
        # Try to reshape, handle various "planar" or extra data issues
        reshaped = None
        for n_remove in (0, 1, 3):  # Try removing 0, 1, or 3 elements at the end
            try:
                if n_remove > 0:
                    match_try = np.delete(match, [-k for k in range(1, n_remove + 1)])
                else:
                    match_try = match
                
                reshaped = np.array(match_try).reshape(-1, 11)
               
                vibs_list.append(reshaped)
                break
            except ValueError:
                continue
        else:
            raise ValueError(f"Could not reshape vibrational data for block {i}: shape {match.shape}")
    
    # Stack all vibrational data
    vibs = np.vstack(vibs_list)
    np.set_printoptions(threshold=np.inf)
    
    # Organize by atom index (assumes first column is atom index)
    try:
        final_atom = int(float(vibs[-1][0]))
    except Exception as e:
        raise ValueError(f"Could not determine final atom index: {e}")

    ordered_vibs = [vibs[vibs[:,0] == str(i)] for i in range(1, final_atom+1)]
    vibs_df = vib_array_list_to_df(ordered_vibs) if ordered_vibs else None
    return vibs_df


def df_list_to_dict(df_list):
    my_dict={}
    for name,df in zip(Names.DF_LIST.value,df_list):
        my_dict[name]=df
    return my_dict



def process_gaussian_energy_text(energy_string):
    cut=re.split('SCF Done', energy_string)[1]
    energy=extract_lines_from_text(cut, re_expression=ReExpressions.FLOAT.value)
    data = np.array([[energy[0]]])
    df = pd.DataFrame(data, columns=['energy'])
    df=df.astype(float)

    return df





def gauss_file_handler(gauss_filename, export=False):
    string_report=''
    
    
    with open(os.path.abspath(gauss_filename)) as f:
        log_file_lines = f.read()

    try:
        charge_df = process_gaussian_charge_text(log_file_lines)
        charge_df=charge_df.astype(float)
    except Exception as e:
        charge_df = pd.DataFrame()
        print(f"{gauss_filename}: Error processing charge: {e}")
        string_report+=f"{gauss_filename}: Error processing charge: {e}\n"
    try:
        dipole_df = process_gaussian_dipole_text(log_file_lines)
        dipole_df=dipole_df.astype(float)
    except Exception as e:
        dipole_df = pd.DataFrame()
        print(f"{gauss_filename}: Error processing dipole: {e}")
        string_report+=f"{gauss_filename}: Error processing dipole: {e}\n"
    try:
        pol_df = process_gaussian_pol_text(log_file_lines)
        pol_df=pol_df.astype(float)
    except Exception as e:
        pol_df = pd.DataFrame()
        print(f"{gauss_filename}: Error processing polarization: {e}")
        string_report+=f"{gauss_filename}: Error processing polarization: {e}\n"
    try:
        gauss_data = gauss_first_split(log_file_lines)
        standard_orientation_df = process_gaussian_standard_orientation_text(gauss_data[2])
        
        
    except Exception as e:
        standard_orientation_df = pd.DataFrame()
        print(f"{gauss_filename}: Error processing standard orientation: {e}")
        string_report+=f"{gauss_filename}: Error processing standard orientation: {e}\n"
    try:
        frequency_str = process_gaussian_vibs_string(log_file_lines)
    except Exception as e:
        frequency_str = pd.DataFrame()
        print(f"{gauss_filename}: Error processing vibrations: {e}")
        string_report+=f"{gauss_filename}: Error processing vibrations: {e}\n"
    try:
        info_df = process_gaussian_info(frequency_str)
        info_df=info_df.astype(float)
        
    except Exception as e:
        info_df = pd.DataFrame()  # or some default DataFrame
        print(f"{gauss_filename}: Error processing info: {e}")
        string_report+=f"{gauss_filename}: Error processing info: {e}\n"

    try:
        vibs_df = process_gaussian_frequency_string(frequency_str)
    except Exception as e:
        vibs_df = pd.DataFrame() # or some default DataFrame
        print(f"{gauss_filename}: Error processing frequency: {e}")
        string_report+=f"{gauss_filename}: Error processing frequency: {e}\n"
    
    try:
        hirsh_charge_df, cm5_charge_df = process_hirshfeld_charges(log_file_lines)
        hirsh_charge_df=hirsh_charge_df.astype(float)
        cm5_charge_df=cm5_charge_df.astype(float)
        print(f"{gauss_filename}: NBO, Hirshfeld and CM5 charges")
        # charge_df={'nbo_charge':charge_df,'hirsh_charge':hirsh_charge_df,'cm5_charge':cm5_charge_df} ### need to edit this -- working

    except Exception as e:
        hirsh_charge_df = pd.DataFrame({'hirshfeld_charge': [np.nan] * standard_orientation_df.shape[0]})
        cm5_charge_df = pd.DataFrame({'cm5_charge': [np.nan] * standard_orientation_df.shape[0]})
        print(f"{gauss_filename}: NBO charge only")
        string_report+=f"{gauss_filename}: NBO charge only\n"

    try:
        ev=extract_homo_lumo_gap(gauss_filename=gauss_filename)
        ev_df=pd.DataFrame({'energy': [ev]})

    except:
        ev_df=pd.DataFrame()
        print(f'{gauss_filename}: eV Gap was not found')
    
    try:
        concatenated_df = pd.concat([standard_orientation_df, dipole_df, pol_df, ev_df ,charge_df, hirsh_charge_df,cm5_charge_df, info_df, vibs_df], axis=1)
  
    except Exception as e:
        concatenated_df = pd.DataFrame()  # or some default DataFrame
        print(f"{gauss_filename}: Error concatenating data: {e}")
        string_report+=f"{gauss_filename}: Error concatenating data: {e}\n"

    return concatenated_df, string_report


def save_to_feather(df, filename):
    # Create a DataFrame from the list of strings
    # Each string becomes a column. Since all columns must be of the same length,
    # you may need to handle this if your strings are of different lengths.
    # Define column names that correspond to the data in string_list
    # column_names = ['Standard_Orientation', 'Dipole', 'Polarizability', 'Frequency', 'Charge', 'Energy']

    df['atom'] = df['atom'].astype(str)
    feather_filename = filename + '.feather'
    # df.columns = range(df.shape[1]) # [str(i) for i in range(df.shape[1])]
    # df.columns = df.columns.map(str)
    df.to_feather(feather_filename)

    print(f"Data saved to {feather_filename}")
    string_report=f"Data saved to {feather_filename}\n"
    return string_report

def logs_to_feather(dir_path, feather_dir = 'feather_files'):
    """
    Processes Gaussian log files in the given directory and saves them as Feather files.
    Default place to save feather files is a subdirectory named 'feather_files'.
    """
    string_report=''
    failed_files_string=None
    os.chdir(dir_path)
    if not os.path.exists(feather_dir):
        os.mkdir(feather_dir)

    for file in os.listdir(dir_path):
        if file.endswith(".log"):
            try:
                df, gauss_string_report = gauss_file_handler(file)
                string_report+=gauss_string_report
            except Exception as e:
                print(f"Error processing file {file}: {e}")
                string_report+=f"Error processing file {file}: {e}\n"
                failed_files_string+=f"{file}\n"
                continue  # Skip to the next file

            os.chdir('feather_files')
            
            string_report+=save_to_feather(df, file.split('.')[0])  # Assuming you want to remove the .log extension
            os.chdir('..')
        else:
            continue
    if failed_files_string is not None:
        string_report = failed_files_string + 'Check the log files and reported errors for more information.'
        
    os.chdir(dir_path)
    
    print('Done!')
    return string_report


def extract_homo_lumo_gap(gauss_filename):
    with open(gauss_filename, 'r') as file:
        lines = file.readlines()

    homo, lumo = None, None
    last_occ_line_idx = None

    # Find index of last 'Alpha occ.' line
    for i, line in enumerate(lines):
        if "Alpha  occ. eigenvalues" in line:
            last_occ_line_idx = i

    if last_occ_line_idx is None:
        return None

    # Get HOMO from the last occ line
    homo_values = re.findall(r"-?\d+\.\d+", lines[last_occ_line_idx])
    homo = float(homo_values[-1])

    # LUMO is the first number from the first 'Alpha virt.' line after last occ
    for j in range(last_occ_line_idx + 1, len(lines)):
        if "Alpha virt. eigenvalues" in lines[j]:
            lumo_values = re.findall(r"-?\d+\.\d+", lines[j])
            lumo = float(lumo_values[0])
            break

    if homo is not None and lumo is not None:
      
        gap_hartree = lumo - homo
        gap_ev = gap_hartree * 27.2114
        return gap_ev
    
def logs_to_dict(dir_path):
    """
    Processes Gaussian log files in the given directory, returning a dictionary of DataFrames.

    Parameters:
    dir_path (str): Path to the directory containing Gaussian log files.

    Returns:
    dict: A dictionary where keys are file names (without .log) and values are DataFrames.
    """
    file_dict = {}
    error_log = {}

    os.chdir(dir_path)

    for file in os.listdir(dir_path):
        if file.endswith(".log"):
            try:
                df, gauss_string_report = gauss_file_handler(file)
                file_dict[file.split('.')[0]] = df
                dict_list=df_file_handler(df)
            except Exception as e:
                print(f"Error processing file {file}: {e}")
                error_log[file] = str(e)  # Store error message for debugging

    return dict_list, error_log


def df_list_to_dict(df_list):
    my_dict={}
    for name,df in zip(['standard_orientation_df', 'dipole_df', 'pol_df','energy', 'info_df'],df_list):
        my_dict[name]=df
    return my_dict

def df_file_handler(data):
    # Read the feather file
    
    data.columns=range(len(data.columns))
    xyz = data.iloc[:, 0:4].dropna()
    dipole_df = data.iloc[:, 4:8].dropna()
    pol_df = data.iloc[:, 8:10].dropna()
    nbo_charge_df = data.iloc[:, 10:11].dropna()
    ## if the whole column is NaN, it will be removed
    hirshfeld_charge_df = data.iloc[:,11:12].dropna().reset_index(drop=True)
    cm5_charge_df = data.iloc[:,12:13].dropna().reset_index(drop=True)
    info_df = data.iloc[:, 13:15].dropna()
    vectors = data.iloc[:, 15:].dropna() 
    xyz.rename(columns={xyz.columns[0]: 'atom', xyz.columns[1]: 'x', xyz.columns[2]: 'y', xyz.columns[3]: 'z'}, inplace=True)
    xyz = xyz.reset_index(drop=True)
    xyz[['x', 'y', 'z']] = xyz[['x', 'y', 'z']].astype(float)
    xyz=xyz.dropna()
    dipole_df.rename(columns={dipole_df.columns[0]: 'dip_x', dipole_df.columns[1]: 'dip_y', dipole_df.columns[2]: 'dip_z', dipole_df.columns[3]: 'total_dipole'}, inplace=True)
    dipole_df=dipole_df.astype(float)
    dipole_df=dipole_df.dropna()
    dipole_df=dipole_df.reset_index(drop=True)
    nbo_charge_df.rename(columns={nbo_charge_df.columns[0]: 'charge'}, inplace=True)
    nbo_charge_df=nbo_charge_df.astype(float)
    nbo_charge_df=nbo_charge_df.dropna()
    nbo_charge_df=nbo_charge_df.reset_index(drop=True)
    charge_dict={'nbo':nbo_charge_df, 'hirshfeld':hirshfeld_charge_df, 'cm5':cm5_charge_df}
    pol_df.rename(columns={pol_df.columns[0]: 'aniso', pol_df.columns[1]: 'iso'}, inplace=True)
    pol_df=pol_df.astype(float)
    pol_df=pol_df.dropna()
    pol_df=pol_df.reset_index(drop=True)
    info_df.rename(columns={info_df.columns[0]: 'Frequency', info_df.columns[1]: 'IR'}, inplace=True)
    info_df=info_df.astype(float)
    info_df=info_df.dropna()
    info_df=info_df.reset_index(drop=True)
    
    def split_to_dict(dataframe):
        
        num_columns = dataframe.shape[1]
        num_dfs = num_columns // 3  # Integer division to get the number of 3-column dataframes
        dfs_dict = {}
        for i in range(num_dfs):
            start_col = i * 3
            end_col = start_col + 3
            key = f'vibration_atom_{i + 1}'
            dfs_dict[key] = np.array(dataframe.iloc[:, start_col:end_col].values.astype(float))  # Storing as a NumPy array
            
        return dfs_dict
    
    vib_dict=split_to_dict(vectors)
    df_list=[xyz ,dipole_df, pol_df, info_df]
    df_list=[df.dropna(how='all') for df in df_list if not df.empty]
    df_dict=df_list_to_dict(df_list)
    dict_list=[df_dict,vib_dict,charge_dict]

    return dict_list

