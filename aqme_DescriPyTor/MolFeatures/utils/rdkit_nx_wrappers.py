import io
from PIL import Image
from networkx import relabel_nodes, Graph

from rdkit.Chem import rdFMCS, MolFromSmiles, MolFromSmarts, Draw, RWMol
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.rdmolfiles import MolToXYZFile
from rdkit.Chem.rdmolops import AddHs
from rdkit.Chem.rdchem import Atom

from .constants import FileExtensions

def draw_rdkit_molecule(molecule, x_dim=512, y_dim=512, atom_indices=True):
    drawer=Draw.MolDraw2DCairo(x_dim, y_dim)
    drawer.drawOptions().useBWAtomPalette()
    drawer.drawOptions().addAtomIndices=atom_indices
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()
    img=Image.open(io.BytesIO(img_bytes))
    img.show()

def rdkit_molecule_to_nx_graph(molecule):
    G=Graph()
    for atom in molecule.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   atom_symbol=atom.GetSymbol(),
                   ring=atom.IsInRing(),
                   charge=atom.GetFormalCharge())
    for bond in molecule.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType(),
                   ring=bond.IsInRing())
    return G
    
def nx_graph_to_rdkit_mol(G):
    dummy_molecule=MolFromSmiles('')
    molecule_writer=RWMol(dummy_molecule)
    for _, node_data in G.nodes(data=True):
        new_atom=Atom(node_data['atomic_num'])
        new_atom.SetFormalCharge(node_data['charge'])
        molecule_writer.AddAtom(new_atom)
    for edge in G.edges(data=True):
        molecule_writer.AddBond(edge[0], edge[1], edge[2]['bond_type'])
    actual_molecule=molecule_writer.GetMol()
    return actual_molecule

def rdkit_mol_to_xyz(molecule, molecule_name, output_filename=None):
    output_filename=molecule_name+FileExtensions.XYZ.value if output_filename==None else output_filename
    molecule.UpdatePropertyCache(strict=False)    
    molecule=AddHs(molecule)
    EmbedMolecule(molecule)
    MMFFOptimizeMolecule(molecule)
    MolToXYZFile(molecule, output_filename)

def get_mcs_mol_from_smiles_rdkit(smiles_list):
    molecules=tuple([MolFromSmiles(smile_str) for smile_str in smiles_list])
    mcs_search_results=rdFMCS.FindMCS(molecules, bondCompare=rdFMCS.BondCompare.CompareAny)
    mcs_smarts=mcs_search_results.smartsString
    mcs_mol=MolFromSmarts(mcs_smarts) #mcs_search_results.queryMol
    return mcs_mol


