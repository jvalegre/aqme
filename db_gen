import glob,os
import numpy as np
import pandas as pd
from pandas import DataFrame
from rdkit import Chem,DataStructs
from rdkit.Chem import rdChemReactions,AllChem,Lipinski,Descriptors
import openbabel as ob

# From a single sdf file containing all the molecules, this generates one individual
# sdf file for each molecule. This is useful to generate RDKit-optimized 3D structures

# set working directory
w_dir = 'C:\Google Drive\Rob Paton CSU\Project Pd Lily\Arene auto generation'

# Set True if you want to generate new sdf files
new_sdf = True

# Method to generate conformers and get their E to remove duplicates and generate files.
# Options: RDKit and xtb
conform_gen_method = RDKit

# Number of conformers you want to optimize and calculate E, creating extra SDF files
n_opt_conf = 10

# Energy threshold to remove duplicates in the conformer generation
threshold = 0.2

# Set a basename (for the first sdf, it doesn't matter)
basename = 'arenes'

# Especify where are the sdf files containing the molecules
os.chdir(w_dir+'\SDF_input_files')
sdf_files = glob.glob('*.sdf')

# If more types of atoms are included, add them here. This is to include the atoms in the GENECP option
possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I']

if new_sdf == True:
    i = 1
    for file in sdf_files:
        os.chdir(w_dir+'\SDF_input_files')
        suppl = Chem.SDMolSupplier(file)
        for mol in suppl:
            valid_structure = True
            # First filter: number of rotatable bonds
            if Chem.Lipinski.NumRotatableBonds(mol) < 3:
                # Second filter: molecular weight
                if Chem.Descriptors.MolWt(mol) < 500:
                    # Third filter: this filters salts off (2 separated components)
                    if len(Chem.MolToSmiles(mol).split('.')) == 1:
                        for atom in mol.GetAtoms():
                            #Fourth filter: atoms outside the scope chosen in 'possible_atoms'
                            if atom.GetSymbol() not in possible_atoms:
                                valid_structure = False
                    else: valid_structure = False
                else: valid_structure = False
            else: valid_structure = False
            if valid_structure == True:
                mol = Chem.AddHs(mol)
                cids = AllChem.EmbedMultipleConfs(mol, numConfs=30, numThreads=0)
                # res stores the E values of all the conformers
                res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
                # This part creates a list containing the energies and a list that contains
                # the indexes of those energies from res
                energy_list, res_index = [],[]
                for j in range(len(res)):
                    energy_list.append(res[j][1])
                    res_index.append(j)
                # This part removes duplicates
                print(res_index)
                print(energy_list)
                # This part sorts all the energies and res indexes starting from the lowest energy value
                energy_list, res_index = (list(t) for t in zip(*sorted(zip(energy_list, res_index))))
                for j in reversed(range(1,len(res))):
                    if energy_list[j] > energy_list[j-1]-threshold\
                    and energy_list[j] < energy_list[j-1]+threshold:
                        energy_list.remove(energy_list[j])
                        res_index.remove(res_index[j])
                print(res_index)
                print(energy_list)

                rmslist = []
                AllChem.AlignMolConformers(mol, RMSlist=rmslist)
                # Especify where to separated sdf files with RDKit-optimized geometries will be saved
                os.chdir(w_dir+'\SDF_output_files')

                for conf in range(n_opt_conf):
                    if conf < len(res_index):
                        print(conf)
                        w = Chem.SDWriter(basename+str(i)+'_CONF'+str(conf)+'.sdf')
                        w.write(mol, confId=res_index[conf])
                        w.close()
                i = i+1
else: print('No new structures were generated.')

# Generate an excel file with the info of the original molecules (this will be useful
# for future comparisons to the final structures)

# set working directory
w_dir = 'C:\Google Drive\Rob Paton CSU\Project Pd Lily\Arene auto generation'

# If more types of atoms are included, add them here. This is to include the atoms in the GENECP option
possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I']

os.chdir(w_dir+'\SDF_output_files')
sdf_to_com_files = glob.glob('*.sdf')

IDs, names, smiles, fingerprints = [],[],[],[]
for file in sdf_to_com_files:
    os.chdir(w_dir+'\SDF_output_files')
    ID,name = '',''
    f = open(file,"r")
    readlines = f.readlines()
    for i in range(len(readlines)):
        if readlines[i].find('>  <ID>') > -1:
            ID = readlines[i+1].split()[0]
        if readlines[i].find('>  <NAME>') > -1:
            name = readlines[i+1].split()[0]
    suppl = Chem.SDMolSupplier(file)
    smile = Chem.MolToSmiles(suppl[0])
    fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile), radius=2, nBits = 40)
    arr = np.zeros((0,), dtype=np.int8)
    fp_np = DataStructs.ConvertToNumpyArray(fp, arr)
    f.close()
    IDs.append(ID)
    names.append(name)
    smiles.append(smile)
    fingerprints.append(arr)

os.chdir(w_dir+'\Excel_database')
df_dict = {'IDs': IDs, 'Names': names, 'Smiles': smiles, 'Fingerprints': fingerprints}
df_components = DataFrame.from_dict(df_dict, orient='index')
df_components = df_components.transpose()
export_excel = df_components.to_csv('Initial_structures.csv', index = None, header=True)


# Generate com files with Link1 from the previous sdf files

# Set True if you want to generate new com files
new_sdf = False

# set working directory
w_dir = 'C:\Google Drive\Rob Paton CSU\Project Pd Lily\Arene auto generation'

# If more types of atoms are included, add them here. This is to include the atoms in the GENECP option
possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I']

obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("sdf", "com")
mol = ob.OBMol()

os.chdir(w_dir+'\SDF_output_files')
sdf_to_com_files = glob.glob('*.sdf')
if new_sdf == True:
    for file in sdf_to_com_files:
        os.chdir(w_dir+'\SDF_output_files')
        ID,name = '',''
        f = open(file,"r")
        readlines = f.readlines()
        for i in range(len(readlines)):
            if readlines[i].find('>  <ID>') > -1:
                ID = readlines[i+1].split()[0]
            if readlines[i].find('>  <NAME>') > -1:
                name = readlines[i+1].split()[0]
        f.close()
        obConversion.ReadFile(mol, file)
        os.chdir(w_dir+'\COM_output_files')
        obConversion.WriteFile(mol, ID+'.com')

        ecp_list,ecp_I = [],False
        read_lines = open(ID+'.com',"r").readlines()
        fileout = open(ID+'.com', "w")
        # Detect if there are I atoms to use genecp or not (to use gen)
        for i in range(4,len(read_lines)):
            if read_lines[i].split(' ')[0] not in ecp_list and read_lines[i].split(' ')[0] in possible_atoms:
                ecp_list.append(read_lines[i].split(' ')[0])
            if read_lines[i].split(' ')[0] == 'I':
               ecp_I = True
        if ecp_I == False:
            genecp = 'gen'
        if ecp_I == True:
            genecp = 'genecp'
        # Modify these parameters for the Gaussian input files
        n_of_processors = '36'
        memory = '30GB'
        functional = 'wb97xd'

        basis_set = '6-31g(d)'
        basis_set_I = 'LANL2DZ'
        functional_nmr = 'mpw1pw91'
        basis_set_nmr = '6-311+G(d,p)'
        basis_set_I_nmr = 'def2tzvpp'
        # Leave blank if gas phase is used, if not use standard keywords (i.e. scrf=(smd,solvent=toluene))
        solvent = ''
        solvent_nmr = ''
        keywords_opt = functional+'/'+genecp+' '+solvent+' opt freq=noraman '
        keywords_nmr = functional_nmr+'/'+genecp+' '+solvent_nmr+' nmr=GIAO geom=check '
        keywords_nbo = functional+'/'+genecp+' '+solvent+' pop=(savenbos,nbo7read) geom=check'
        nbo_sterics = '$nbo steric bndidx $end'

        #Write input for part 1 (opt + freq)
    #     fileout.write("%chk="+ID+".chk"+"\n")
        fileout.write("%mem="+memory+"\n")
        fileout.write("%nprocshared="+n_of_processors+"\n")
        fileout.write("# "+keywords_opt+"\n")
        fileout.write("\n")
        fileout.write(name+"\n")
        fileout.write("\n")
        for i in range(4,len(read_lines)):
            fileout.write(read_lines[i])
        for i in range(len(ecp_list)):
            if ecp_list[i] != 'I':
                fileout.write(ecp_list[i]+' ')
        fileout.write('0\n')
        fileout.write(basis_set+'\n')
        fileout.write('****\n')
        if ecp_I == False:
            fileout.write('\n')
        else:
            fileout.write('I     0\n')
            fileout.write(basis_set_I+'\n')
            fileout.write('****\n\n')
            fileout.write('I 0\n')
            fileout.write(basis_set_I+'\n\n')
    #     #Write input for part 2 (NMR)
    #     fileout.write("--Link1--"+"\n")
    #     fileout.write("%chk="+ID+".chk"+"\n")
    #     fileout.write("%mem="+memory+"\n")
    #     fileout.write("%nprocshared="+n_of_processors+"\n")
    #     fileout.write("# "+keywords_nmr+"\n")
    #     fileout.write("\n")
    #     fileout.write(name+"\n")
    #     fileout.write("\n")
    #     fileout.write(read_lines[4])
    #     fileout.write("\n")
    #     for i in range(len(ecp_list)):
    #         if ecp_list[i] != 'I':
    #             fileout.write(ecp_list[i]+' ')
    #     fileout.write('0\n')
    #     fileout.write(basis_set_nmr+'\n')
    #     if ecp_I == False:
    #         fileout.write('\n')
    #     else:
    #         fileout.write('****\n')
    #         fileout.write('I     0\n')
    #         fileout.write(basis_set_I_nmr+'\n')
    #         fileout.write('****\n\n')
    #         fileout.write('I 0\n')
    #         fileout.write(basis_set_I_nmr+'\n\n')
    #     #Write input for part 3 (NBO)
    #     fileout.write("--Link1--"+"\n")
    #     fileout.write("%chk="+ID+".chk"+"\n")
    #     fileout.write("%mem="+memory+"\n")
    #     fileout.write("%nprocshared="+n_of_processors+"\n")
    #     fileout.write("# "+keywords_nbo+"\n")
    #     fileout.write("\n")
    #     fileout.write(name+"\n")
    #     fileout.write("\n")
    #     fileout.write(read_lines[4])
    #     fileout.write("\n")
    #     for i in range(len(ecp_list)):
    #         if ecp_list[i] != 'I':
    #             fileout.write(ecp_list[i]+' ')
    #     fileout.write('0\n')
    #     fileout.write(basis_set+'\n')
    #     if ecp_I == False:
    #         fileout.write('\n')
    #     else:
    #         fileout.write('****\n')
    #         fileout.write('I     0\n')
    #         fileout.write(basis_set_I+'\n')
    #         fileout.write('****\n\n')
    #         fileout.write('I 0\n')
    #         fileout.write(basis_set_I+'\n\n')
    #     fileout.write(nbo_sterics+"\n")
    #     fileout.write("\n")
        fileout.close()
else: print('No new com files were generated.')

print('test')
print('test-shree')
