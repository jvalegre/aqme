#####################################################.
#       This file stores functions related to       #
#          metal templates used in CSEARCH          #
#####################################################.

import sys
from importlib.resources import files
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDistGeom, rdMolAlign
from aqme.utils import get_conf_RMS, load_sdf

TEMPLATES_PATH = files("aqme").joinpath("templates")


def template_embed(self, mol, complex_type, metal_idx, maxsteps, heavyonly, maxmatches, name):
    """Apply template-based embedding to metal complexes.
    
    Automatically selects and applies the appropriate embedding function based on
    metal coordination number. Handles geometry optimization and filtering of
    similar conformers.

    Args:
        self: Class instance containing args with logger
        mol (rdkit.Chem.rdchem.Mol): Molecule to embed
        complex_type (str): Type of complex geometry (e.g. "squareplanar")
        metal_idx (list): Indices of metal atoms in molecule
        maxsteps (int): Maximum optimization steps
        heavyonly (bool): If True, only consider heavy atoms for RMSD
        maxmatches (int): Maximum number of atom matches for RMSD
        name (str): Molecule identifier for logging
        
    Returns:
        list: Items containing embedded conformers and associated data
        
    Note:
        The embedding function is chosen based on the number of metal neighbors:
        - 2: Linear geometry
        - 3: Trigonal planar
        - 4: Square planar
        - 5: Square pyramidal
    """
    # Map coordination number to embedding function
    embed_functions = {
        2: two_embed,
        3: three_embed,
        4: four_embed,
        5: five_embed
    }

    # Load appropriate geometry template
    template = load_template(complex_type, self.args.log)

    # Generate embeddings based on metal coordination
    neighbours = calc_neighbours(mol, metal_idx)
    embed_func = embed_functions[len(neighbours)]
    results = embed_func(
        mol, template, metal_idx, neighbours, 
        name, maxsteps, self.args.log
    )

    # Filter similar conformers if multiple were generated
    molecules = results[0]
    if len(molecules) > 1:
        # Identify conformers to ignore based on RMSD
        ignored_indices = []
        for idx, mol_filter in enumerate(molecules):
            if filter_template_mol(mol_filter, molecules, heavyonly, maxmatches):
                ignored_indices.append(idx)
        
        # Remove ignored conformers from results
        results = [
            item for idx, item in enumerate(results) 
            if idx not in ignored_indices
        ]

    return results


def template_embed_optimize(target, template, metal_idx, maxsteps, log, template_n=None, cumulative_algMap=None):
    """Embed and optimize a molecular conformation using template constraints.

    Embeds a new conformation into a molecule, applies template-based distance
    constraints, performs UFF optimization, and aligns the result to the template.
    Handles special cases for metal coordination environments.

    Args:
        target (rdkit.Chem.rdchem.Mol): Target molecule for embedding
        template (rdkit.Chem.rdchem.Mol): Template molecule for geometry constraints
        metal_idx (list): Indices of metal atoms
        maxsteps (int): Maximum optimization steps
        log: Logger object for status messages
        template_n (int, optional): Template number for metal coordination
        cumulative_algMap (list, optional): Accumulated alignment mappings

    Returns:
        tuple: (molecule, coord_map, alg_map, conf_id, cumulative_algMap) where:
            - molecule: Optimized molecule
            - coord_map: Coordinate mapping dictionary
            - alg_map: Alignment atom mapping
            - conf_id: Conformer ID (-1 if embedding failed)
            - cumulative_algMap: Updated alignment mappings

    Note:
        Uses distance constraints from template and UFF optimization
        with a force constant of 10000 for template constraints.
    """
    # Initialize parameters
    SEED = -1
    FORCE_CONSTANT = 10000
    MAX_ALIGN_ITERS = 100
    
    if cumulative_algMap is None:
        cumulative_algMap = []

    # Get atom mappings between target and template
    coord_map, alg_map, cumulative_algMap = get_mappings(
        target, template, metal_idx, 
        cumulative_algMap, template_n, 
        conformer_id=-1
    )

    # Prepare molecule and embed conformation
    molecule = Chem.AddHs(target)
    conf_id = rdDistGeom.EmbedMolecule(
        molecule, 
        coordMap=coord_map, 
        randomSeed=SEED
    )

    if conf_id < 0:
        log.write("Could not embed molecule.")
        return molecule, None, None, conf_id, cumulative_algMap

    # Setup force field with template constraints
    forcefield = Chem.UFFGetMoleculeForceField(molecule, confId=conf_id)
    
    # Add distance constraints from template
    constraints = get_distance_constrains(coord_map)
    for atom1_idx, atom2_idx, target_dist in constraints:
        forcefield.AddDistanceConstraint(
            atom1_idx, atom2_idx,
            target_dist, target_dist,
            FORCE_CONSTANT
        )

    # Optimize geometry
    forcefield.Initialize()
    forcefield.Minimize(maxIts=maxsteps)

    # Align optimized structure to template
    rdMolAlign.AlignMol(
        molecule, template,
        atomMap=alg_map,
        reflect=True,
        maxIters=MAX_ALIGN_ITERS
    )

    return molecule, coord_map, alg_map, conf_id, cumulative_algMap


def filter_template_mol(molecule_new, mol_objects, heavyonly, max_matches):
    """Filter molecules based on RMSD similarity to existing conformers.

    Determines if a new molecular conformer should be kept by comparing its RMSD
    to a list of existing conformers. Conformers that are too similar to existing
    ones are filtered out.

    Args:
        molecule_new (rdkit.Chem.rdchem.Mol): New conformer to evaluate
        mol_objects (list): List of existing conformers to compare against
        heavyonly (bool): If True, only consider non-H atoms for RMSD
        max_matches (int): Maximum number of atom matches for RMSD calculation

    Returns:
        bool: True if conformer should be kept (unique), False if too similar

    Note:
        Uses an RMSD threshold of 0.5 Å. Conformers with RMSD < 0.5 Å to any
        existing conformer are considered duplicates and filtered out.
    """
    RMSD_THRESHOLD = 0.5  # Angstroms
    
    # Keep if this is the first conformer
    if not mol_objects:
        return True

    # Compare against existing conformers
    for mol in mol_objects:
        rms = get_conf_RMS(
            mol, molecule_new,
            -1,  # Use first conformer
            -1,  # Use first conformer 
            heavyonly,
            max_matches
        )
        # Filter out if too similar to any existing conformer
        if rms < RMSD_THRESHOLD:
            return False
            
    return True


def get_mappings(molecule, template, metal_idx, cumulative_algMap, template_n, conformer_id=-1):
    """Get atom mappings between molecule and template.

    Creates coordinate and atom mappings between a molecule and its template,
    with special handling for metal complexes where ligands may have identical
    atom types but different positions.

    Args:
        molecule (rdkit.Chem.rdchem.Mol): Molecule to map
        template (rdkit.Chem.rdchem.Mol): Template molecule
        metal_idx (list): Indices of metal atoms
        cumulative_algMap (list): Previous alignment mappings
        template_n (int, optional): Template variant number (1 or 2)
        conformer_id (int, default=-1): Conformer ID to use from template

    Returns:
        tuple: (coordMap, algMap, cumulative_algMap) where:
            - coordMap (dict): Maps atom indices to 3D coordinates
            - algMap (list): Atom index pairs for alignment
            - cumulative_algMap (list): Updated alignment mappings

    Note:
        Special handling for metal complexes with identical ligands to ensure
        correct mapping of ligand positions based on template variant.
    """
    # Get initial mappings
    match = molecule.GetSubstructMatch(template)
    conformer = template.GetConformer(conformer_id)
    
    # Create initial coordinate and alignment maps
    coord_map = {}
    alg_map = []
    for template_idx, atom_idx in enumerate(match):
        alg_map.append((atom_idx, template_idx))
        coord_map[atom_idx] = conformer.GetAtomPosition(template_idx)

    # Handle special case of identical ligands
    if alg_map in cumulative_algMap and template_n is not None:
        # Separate metal and ligand mappings
        neigh_coord_map, neigh_alg_map = {}, []
        metal_coord_map, metal_alg_map = {}, []
        
        for idx, match in alg_map:
            if idx == metal_idx[0]:
                metal_alg_map.append((idx, match))
                metal_coord_map[idx] = coord_map[idx]
            else:
                neigh_alg_map.append((idx, match))
                neigh_coord_map[idx] = coord_map[idx]

        # Define ligand reordering based on template variant
        rel_idx = [2, 3, 1] if template_n == 1 else [3, 1, 2]
        
        # Rebuild coordinate map with reordered ligands
        new_coord_map = {}
        coord_keys = list(neigh_coord_map.keys())
        
        # Map ligand coordinates
        new_coord_map[coord_keys[0]] = neigh_coord_map[coord_keys[0]]
        for i, idx in enumerate(rel_idx):
            new_coord_map[coord_keys[idx]] = neigh_coord_map[coord_keys[i + 1]]
        
        # Add metal coordinates
        new_coord_map[metal_idx[0]] = metal_coord_map[metal_idx[0]]
        
        # Rebuild alignment map with reordered ligands
        new_alg_map = [
            neigh_alg_map[0],  # First ligand
            (neigh_alg_map[rel_idx[0]][0], neigh_alg_map[1][1]),
            (neigh_alg_map[rel_idx[1]][0], neigh_alg_map[2][1]),
            (neigh_alg_map[rel_idx[2]][0], neigh_alg_map[3][1]),
            metal_alg_map[0]   # Metal center
        ]
        
        # Update mappings
        coord_map = new_coord_map
        alg_map = new_alg_map

    # Track this mapping
    cumulative_algMap.append(alg_map)

    return coord_map, alg_map, cumulative_algMap


def load_template(complex_type, log):
    """Load molecular template for a given complex geometry type.

    Loads the appropriate geometric template for a given metal complex type
    from the template library. Verifies template accessibility and handles
    errors appropriately.

    Args:
        complex_type (str): Type of complex geometry. One of:
            - "squareplanar": Square planar geometry
            - "squarepyramidal": Square pyramidal geometry
            - "linear": Linear geometry
            - "trigonalplanar": Trigonal planar geometry
        log: Logger object for status messages

    Returns:
        rdkit.Chem.rdchem.Mol: Template molecule for the specified geometry

    Raises:
        SystemExit: If template folder is not found
        
    Note:
        Some templates are shared between different geometries:
        - template-4-and-5.sdf: Used for both square planar and pyramidal
        - template-2.sdf: Used for linear geometries
        - template-3.sdf: Used for trigonal planar geometries
    """
    # Map complex types to template files
    TEMPLATE_FILES = {
        "squareplanar": "template-4-and-5.sdf",
        "squarepyramidal": "template-4-and-5.sdf",
        "linear": "template-2.sdf",
        "trigonalplanar": "template-3.sdf"
    }

    # Check template folder exists
    if not TEMPLATES_PATH.exists():
        log.write(
            "x  The templates folder was not found, probably due to a "
            "problem while installing AQME"
        )
        log.finalize()
        sys.exit()

    # Load appropriate template file
    template_file = TEMPLATES_PATH.joinpath(TEMPLATE_FILES[complex_type])
    templates = load_sdf(str(template_file))
    
    # Return last template in file (most relevant one)
    return templates[-1]


def calc_neighbours(molecule, metals_idx):
    """Process metal centers and return their neighboring atoms.

    Identifies the first metal atom in the molecule, replaces it with an
    appropriate surrogate atom (Si or I) based on coordination number,
    and returns its neighboring atoms. Handles charge adjustment for
    higher coordination numbers.

    Args:
        molecule (rdkit.Chem.rdchem.Mol): Molecule containing metal centers
        metals_idx (list): List of metal atom indices

    Returns:
        list: Neighboring atoms of the first metal center found

    Note:
        Uses surrogate atoms to fit templates:
        - Si (atomic num 14) for 4 and 5 coordination
        - I (atomic num 53) for 2 and 3 coordination
        
        For 5-coordinate centers, adjusts formal charge to maintain
        valid valence state.
    """
    # Map coordination number to surrogate atomic numbers
    COORD_TO_ATOMIC_NUM = {
        5: 14,  # Silicon for 5-coordinate
        4: 14,  # Silicon for 4-coordinate
        3: 53,  # Iodine for 3-coordinate
        2: 53   # Iodine for 2-coordinate
    }
    
    # Maximum charge adjustment attempts for 5-coordinate
    MAX_CHARGE_ATTEMPTS = 5

    # Find first metal atom and process it
    for atom in molecule.GetAtoms():
        if atom.GetIdx() in metals_idx:
            # Get neighbors and bond count
            neighbours = atom.GetNeighbors()
            n_bonds = len(atom.GetBonds())
            
            # Replace metal with appropriate surrogate
            atom.SetAtomicNum(COORD_TO_ATOMIC_NUM[n_bonds])
            
            # Handle special case of 5-coordinate centers
            if n_bonds == 5:
                # Try increasing charges until valid
                for charge_adj in range(MAX_CHARGE_ATTEMPTS):
                    atom.SetFormalCharge(1 + charge_adj)
                    try:
                        # Test if molecule is valid with this charge
                        mol_test = Chem.Mol(molecule)
                        Chem.SanitizeMol(mol_test)
                        break  # Valid charge found
                    except Chem.AtomValenceException:
                        continue  # Try next charge
            
            return neighbours
    
    # Return empty list if no metal found
    return []


def check_metal_neigh(mol, complex_type, metal_idx_ind, log, valid_template):
    """Validate metal center coordination number against template.
    
    Checks if the metal center's coordination number matches the expected
    number of ligands for the chosen template geometry type.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule containing metal center
        complex_type (str): Type of geometry template:
            - "linear": 2-coordinate linear geometry
            - "trigonalplanar": 3-coordinate trigonal planar
            - "squareplanar": 4-coordinate square planar
            - "squarepyramidal": 5-coordinate square pyramidal
        metal_idx_ind (int): Index of the metal atom to check
        log: Logger object for status messages
        valid_template (bool): Current template validity state
    
    Returns:
        bool: Updated template validity state. False if coordination
        numbers don't match, otherwise maintains input value.
        
    Note:
        Logs an error message if coordination numbers don't match.
    """
    # Define expected coordination numbers
    COORD_NUMS = {
        "linear": 2,
        "trigonalplanar": 3, 
        "squareplanar": 4,
        "squarepyramidal": 5
    }
    
    # Get expected coordination for this geometry
    expected_coord = COORD_NUMS[complex_type]
    
    # Count actual coordination number
    metal_atom = mol.GetAtoms()[metal_idx_ind]
    actual_coord = len(metal_atom.GetNeighbors())
    
    # Check for mismatch
    if actual_coord != expected_coord:
        log.write(
            f"x  The number of neighbours of the metal ({actual_coord}) "
            f"does not match the number of expected neighbours for the "
            f"template selected ({complex_type}). No templates will be "
            f"applied to this system."
        )
        valid_template = False

    return valid_template


def get_distance_constrains(coordMap):
    """Generate distance constraints between mapped atoms.
    
    Creates a list of distance constraints between all pairs of atoms
    in the coordinate mapping, calculating their 3D distances.

    Args:
        coordMap (dict): Maps atom indices to Point3D coordinates

    Returns:
        list: List of tuples (atom1_idx, atom2_idx, distance) for all
            atom pairs in the mapping

    Note:
        Distance between atoms is calculated using Point3D.Distance()
        method from RDKit's geometry utilities.
    """
    # Get list of atom indices
    atom_idxs = list(coordMap.keys())
    
    # Generate constraints for all pairs
    constrains = []
    for k, i in enumerate(atom_idxs):
        # Only process each pair once
        for j in atom_idxs[k + 1:]:
            # Calculate 3D distance
            d = coordMap[i].Distance(coordMap[j])
            constrains.append((i, j, d))
            
    return constrains

# Embedding functions
def two_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log):
    """Generate linear coordination geometry embeddings.
    
    Sets up and performs template-based embedding for linear metal complexes
    with two coordinating atoms. Uses iodine (I) as a surrogate metal center.
    
    Args:
        molecule: Target molecule to embed
        template: Linear geometry template 
        metal_idx: Index of metal atom
        neighbours: List of neighboring atoms to metal
        name: Identifier for generated conformer
        maxsteps: Maximum optimization steps
        log: Logger for messages
    
    Returns:
        Multiple lists containing generated conformer data:
        - List of embedded molecules
        - List of conformer names
        - List of coordinate mappings
        - List of alignment mappings
        - List of templates used
        - List of any atom number overrides
        
    Note:
        Uses template-2.sdf with I (atomic num 53) as surrogate metal.
        Returns empty lists if embedding fails.
    """
    # Set template atom types to match neighbors
    template.GetAtomWithIdx(0).SetAtomicNum(neighbours[0].GetAtomicNum())
    template.GetAtomWithIdx(1).SetAtomicNum(neighbours[1].GetAtomicNum())
    template.GetAtomWithIdx(2).SetAtomicNum(53)  # I as surrogate metal

    # Embed and optimize
    mol_obj, coord_map, alg_map, ci, _ = template_embed_optimize(
        molecule, template, metal_idx, maxsteps, log
    )

    # Track any atom number overrides (none for linear case)
    original_atn_list = [None]

    # Return results if embedding succeeded
    if ci >= 0:
        return [mol_obj], [name], [coord_map], [alg_map], [template], original_atn_list

    # Return empty lists on failure
    return [], [], [], [], [], original_atn_list


def three_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log):
    """Generate trigonal planar coordination geometry embeddings.
    
    Sets up and performs template-based embedding for trigonal planar metal
    complexes with three coordinating atoms. Uses iodine (I) as a surrogate
    metal center.
    
    Args:
        molecule: Target molecule to embed
        template: Trigonal planar geometry template
        metal_idx: Index of metal atom
        neighbours: List of neighboring atoms to metal
        name: Identifier for generated conformer
        maxsteps: Maximum optimization steps
        log: Logger for messages
    
    Returns:
        Multiple lists containing generated conformer data:
        - List of embedded molecules
        - List of conformer names
        - List of coordinate mappings 
        - List of alignment mappings
        - List of templates used
        - List of any atom number overrides
        
    Note:
        Uses template-3.sdf with I (atomic num 53) as surrogate metal.
        Returns empty lists if embedding fails.
    """
    # Set template atom types to match neighbors
    template.GetAtomWithIdx(0).SetAtomicNum(53)  # I as surrogate metal
    template.GetAtomWithIdx(1).SetAtomicNum(neighbours[0].GetAtomicNum())
    template.GetAtomWithIdx(2).SetAtomicNum(neighbours[1].GetAtomicNum())
    template.GetAtomWithIdx(3).SetAtomicNum(neighbours[2].GetAtomicNum())

    # Embed and optimize
    mol_obj, coord_map, alg_map, conf_id, _ = template_embed_optimize(
        molecule, template, metal_idx, maxsteps, log
    )
    
    # Track any atom number overrides (none for trigonal planar case)
    original_atn_list = [None]

    # Return results if embedding succeeded
    if conf_id >= 0:
        return [mol_obj], [name], [coord_map], [alg_map], [template], original_atn_list

    # Return empty lists on failure
    return [], [], [], [], [], original_atn_list

def four_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log):
    """Generate square planar coordination geometry embeddings.
    
    Sets up and performs template-based embedding for square planar metal
    complexes with four coordinating atoms. Uses template-4-and-5.sdf with
    special handling for problematic atoms like arsenic ligands and carbenes.
    Attempts multiple embeddings with different ligand orderings to find
    optimal geometry.
    
    Args:
        molecule (rdkit.Chem.rdchem.Mol): Molecule containing metal center 
        template (rdkit.Chem.rdchem.Mol): Square planar template molecule
        metal_idx (list): Index of metal atom to transform
        neighbours (list): Metal's neighboring atoms 
        name (str): Base name for generated conformers
        maxsteps (int): Maximum optimization iterations
        log: Logger for status messages

    Returns:
        tuple: Lists of generated conformer data:
            - mol_objects: Embedded molecules 
            - name_return: Conformer names
            - coord_maps: Coordinate mappings
            - alg_maps: Alignment mappings
            - mol_templates: Template molecules used
            - original_atn_list: Original atomic numbers (for restoring As)
            
    Note:
        Special handling includes:
        - Arsenic ligands are temporarily replaced with phosphorus
        - N-C-N carbenes are converted to N-C=N+ form
        - Attempts 3 different ligand orderings around metal
        Uses silicon (atomic number 14) as surrogate metal center
    """
    # Initialize result containers
    mol_objects = []
    name_return = []
    coord_maps = []
    alg_maps = []
    mol_templates = []

    # Remove F atoms from the template
    template = Chem.RWMol(template)
    for atom in template.GetAtoms():
        if atom.GetSymbol() == "F":
            template = Chem.RWMol(template)
            template.RemoveAtom(atom.GetIdx())
    template = template.GetMol()

    # Process neighbors and handle special cases
    neighbor_atomic_nums = []
    original_atn = None
    
    # Process each neighbor atom
    for i, neighbor in enumerate(neighbours):
        # Get atomic number
        atom_num = neighbor.GetAtomicNum()
        atom_neigh = molecule.GetAtomWithIdx(neighbor.GetIdx())
        
        # Handle problematic atoms
        if atom_num == 33:  # Replace As with P
            original_atn = [atom_num, neighbor.GetIdx()]
            atom_neigh.SetAtomicNum(15)  # P atomic number
            neighbor_atomic_nums.append(15)
        elif (atom_num == 6 and  # Handle carbene carbons
              atom_neigh.GetTotalValence() == 3):
            # Check for N-C-N carbene pattern
            c_neighbors = atom_neigh.GetNeighbors()
            n_neighbors = [n for n in c_neighbors 
                         if n.GetAtomicNum() == 7]
            
            if len(c_neighbors) == 3 and len(n_neighbors) == 2:
                # Convert N-C-N to N-C=N+ for RDKit compatibility
                c_idx = atom_neigh.GetIdx()
                n_idx = n_neighbors[0].GetIdx()
                
                # Update bond and charges
                bond = molecule.GetBondBetweenAtoms(c_idx, n_idx)
                if bond:
                    bond.SetBondType(Chem.BondType.DOUBLE)
                molecule.GetAtomWithIdx(c_idx).SetFormalCharge(0)
                molecule.GetAtomWithIdx(n_idx).SetFormalCharge(1)
                Chem.SanitizeMol(molecule)
            neighbor_atomic_nums.append(atom_num)
        else:
            neighbor_atomic_nums.append(atom_num)

    # Define possible ligand orderings for square planar geometry
    # Each tuple represents: (L1, L2, L3, L4, M) where M is surrogate metal (Si)
    replacements = [
        (neighbor_atomic_nums[0], neighbor_atomic_nums[1], 
         neighbor_atomic_nums[2], neighbor_atomic_nums[3], 14),
        (neighbor_atomic_nums[0], neighbor_atomic_nums[2],
         neighbor_atomic_nums[3], neighbor_atomic_nums[1], 14),
        (neighbor_atomic_nums[0], neighbor_atomic_nums[3],
         neighbor_atomic_nums[1], neighbor_atomic_nums[2], 14)
    ]
    
    # Try each ligand ordering
    cumulative_algMap = []
    original_atn_list = []
    
    for i, atom_nums in enumerate(replacements):
        # Create template copy and set atomic numbers
        template_mol = Chem.Mol(template)
        for idx, atom_num in enumerate(atom_nums):
            template_mol.GetAtomWithIdx(idx).SetAtomicNum(atom_num)

        # Try embedding with this configuration
        mol_obj, coord_map, alg_map, conf_id, cumulative_algMap = template_embed_optimize(
            molecule, template_mol, metal_idx, maxsteps, log,
            template_n=i, cumulative_algMap=cumulative_algMap
        )

        # Store successful embedding
        if conf_id >= 0:
            name_out = f"{name.split()[0]}_{i}"
            mol_objects.append(mol_obj)
            name_return.append(name_out)
            coord_maps.append(coord_map)
            alg_maps.append(alg_map)
            mol_templates.append(template_mol)
            original_atn_list.append(original_atn)

    return mol_objects, name_return, coord_maps, alg_maps, mol_templates, original_atn_list

def five_embed(molecule, template, metal_idx, neighbours, name, maxsteps, log):
    """Generate square pyramidal coordination geometry embeddings.
    
    Sets up and performs template-based embedding for square pyramidal metal
    complexes with five coordinating atoms. Uses a modified template-4-and-5.sdf
    template and attempts multiple embeddings with different ligand orderings.
    
    Args:
        molecule: Target molecule to embed
        template: Square pyramidal geometry template
        metal_idx: Index of metal atom
        neighbours: List of neighboring atoms to metal
        name: Identifier for generated conformer
        maxsteps: Maximum optimization steps
        log: Logger for messages
    
    Returns:
        Multiple lists containing generated conformer data:
        - List of embedded molecules
        - List of conformer names
        - List of coordinate mappings
        - List of alignment mappings
        - List of templates used
        - List of any atom number overrides
        
    Note:
        Makes 15 embedding attempts with different ligand orderings.
        Uses Si (atomic num 14) as the surrogate metal center.
        Silicon atom is given a +1 formal charge for valence.
    """
    # Initialize result containers
    mol_objects = []
    name_return = []
    coord_maps = []
    alg_maps = []
    mol_templates = []
    counter = 0

    # Get atomic numbers of neighboring atoms
    atomic_numbers = [atom.GetAtomicNum() for atom in neighbours]

    # Define ligand permutations for all possible orientations
    replacements = [
        # Ligand ordering variations (15 total)
        [4, 0, 1, 2, 3],  # Base ordering
        [4, 0, 2, 3, 1],  # Rotate base ligands
        [4, 0, 3, 1, 2],
        [3, 0, 1, 2, 4],  # Move axial ligand to equatorial
        [3, 0, 2, 4, 1],
        [3, 0, 4, 1, 2],
        [2, 0, 1, 4, 3],
        [2, 0, 4, 3, 1],
        [2, 0, 4, 1, 3],
        [1, 0, 4, 2, 3],
        [1, 0, 2, 3, 4],
        [1, 0, 3, 4, 2],
        [0, 4, 1, 2, 3],
        [0, 4, 2, 3, 1],
        [0, 4, 3, 1, 2],
    ]

    # Try each ligand ordering
    original_atn_list = []
    for replacement in replacements:
        # Map ligand atomic numbers according to ordering
        at0, at1, at2, at3, at4 = [atomic_numbers[r] for r in replacement]
        
        # Set template atom types
        template.GetAtomWithIdx(0).SetAtomicNum(at0)
        template.GetAtomWithIdx(1).SetAtomicNum(at1)
        template.GetAtomWithIdx(2).SetAtomicNum(at2)
        template.GetAtomWithIdx(3).SetAtomicNum(at3)
        template.GetAtomWithIdx(4).SetAtomicNum(at4)
        template.GetAtomWithIdx(5).SetAtomicNum(14)  # Silicon
        template.GetAtomWithIdx(5).SetFormalCharge(1)  # +1 charge for valence

        # Attempt embedding with this configuration
        mol_obj, coord_map, alg_map, conf_id, _ = template_embed_optimize(
            molecule, template, metal_idx, maxsteps, log
        )

        # Store successful embedding
        if conf_id >= 0:
            name_out = f"{name}_{counter}"
            mol_objects.append(mol_obj)
            name_return.append(name_out)
            coord_maps.append(coord_map)
            alg_maps.append(alg_map)
            mol_templates.append(template)
            original_atn_list.append(None)
            counter += 1
            
    return mol_objects, name_return, coord_maps, alg_maps, mol_templates, original_atn_list
