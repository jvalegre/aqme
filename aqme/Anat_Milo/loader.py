import os
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

class PairedMoleculeLoader:
    """Manages file ingestion, structural data generation, and property matching."""
    def __init__(self, base_name: str, xyz_path: str, json_path: str, cube_path: str = None):
        self.name = base_name
        self.elements = []
        self.coordinates = []
        self.raw_json = {}
        self.cube_path = cube_path if cube_path and os.path.exists(cube_path) else None
        
        self._load_xyz(xyz_path)
        self._load_json(json_path)
        self._build_connectivity_table()

    def _load_xyz(self, path: str):
        with open(path, 'r') as f:
            lines = f.readlines()
        start_idx = 2 if len(lines) > 2 and lines[0].strip().isdigit() else 0
        for line in lines[start_idx:]:
            tokens = line.strip().split()
            if len(tokens) >= 4:
                self.elements.append(tokens[0])
                self.coordinates.append([float(x) for x in tokens[1:4]])
        self.coordinates = np.array(self.coordinates, dtype=float)
        
        self.xyz_df = pd.DataFrame({
            'atom': self.elements,
            'x': self.coordinates[:, 0],
            'y': self.coordinates[:, 1],
            'z': self.coordinates[:, 2]
        })

    def _load_json(self, path: str):
        with open(path, 'r') as f:
            self.raw_json = json.load(f)

    def _build_connectivity_table(self, threshold_distance: float = 1.82, metal_threshold: float = 2.8):
        """Builds an undirected connectivity table using native VdW cutoff matrices."""
        metal_elements = {'Li', 'Na', 'Mg', 'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Pd', 'Ag', 'Pt', 'Au', 'Hg'}
        num_atoms = len(self.elements)
        if num_atoms < 2:
            self.bonds_df = pd.DataFrame(columns=[0, 1])
            return

        dist_matrix = squareform(pdist(self.coordinates))
        bonds = []
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = dist_matrix[i, j]
                if dist == 0: continue
                
                is_metal = self.elements[i] in metal_elements or self.elements[j] in metal_elements
                limit = metal_threshold if is_metal else threshold_distance
                
                if dist <= limit:
                    if self.elements[i] == 'H' and self.elements[j] == 'H': continue
                    if (self.elements[i] == 'H' or self.elements[j] == 'H') and dist >= 1.5: continue
                    bonds.append([i + 1, j + 1]) # Convert to 1-based indexing for internal utility engines
                    
        self.bonds_df = pd.DataFrame(bonds, columns=[0, 1])

    def get_prop(self, key: str, default=None):
        return self.raw_json.get(key, default)