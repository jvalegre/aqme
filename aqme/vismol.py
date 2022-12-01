import py3Dmol
from rdkit import Chem
from ipywidgets import interact
import ipywidgets
import os
import subprocess
import time


from aqme.utils import (
    load_variables,
    mol_from_sdf_or_mol_or_mol2,
)


class vismol:
    """
    Class to visualize the molecules from SDF files
    """

    def __init__(self, **kwargs):

        start_time_overall = time.time()
        # load default and user-specified variables
        self.args = load_variables(kwargs, "vismol")

        # write input files
        for file in self.args.files:
            name = file.replace("/", "\\").split("\\")[-1].split(".")[0]
            if file.split(".")[1].lower() in ["sdf", "xyz", "pdb"]:
                sdf_files = []
                if file.split(".")[1].lower() == "xyz":
                    # separate the parent XYZ file into individual XYZ files
                    command_xyz = [
                        "obabel",
                        "-ixyz",
                        file,
                        "-osdf",
                        f"-O{name}.sdf",
                    ]
                    subprocess.run(
                        command_xyz,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    sdf_files.append(f"{name}.sdf")

                elif file.split(".")[1].lower() == "pdb":
                    command_pdb = ["obabel", "-ipdb", file, "-osdf", f"-O{name}.sdf"]
                    subprocess.run(
                        command_pdb,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    sdf_files.append(f"{name}.sdf")
                else:
                    sdf_files.append(file)

            for sdf_file in sdf_files:
                self.confs = mol_from_sdf_or_mol_or_mol2(sdf_file, "qprep")
                interact(
                    self.style_selector,
                    idx=ipywidgets.IntSlider(min=0, max=len(self.confs) - 1, step=1),
                    s=ipywidgets.Dropdown(
                        options=["line", "stick", "sphere"],
                        value="line",
                        description="Style:",
                    ),
                )
                if file.split(".")[1].lower() in ["xyz", "pdb"]:
                    try:
                        # delete SDF files when the input was an XYZ/PDB file
                        os.remove(sdf_file)
                    except PermissionError:
                        pass # in Windows, the generated SDF file cannot be removed

        if self.args.verbose:
            elapsed_time = round(time.time() - start_time_overall, 2)
            self.args.log.write(f"\nTime VIZMOL: {elapsed_time} seconds\n")
        self.args.log.finalize()

    def style_selector(self, idx, s):
        conf = self.confs[idx]
        return self.MolTo3DView(conf, style=s).show()

    def MolTo3DView(
        self, mol, size=(300, 300), style="sphere", surface=False, opacity=0.5
    ):
        """
        Draw molecule in 3D

        Args:
        ----
            mol: rdMol, molecule to show
            size: tuple(int, int), canvas size
            style: str, type of drawing molecule
                        style can be 'line', 'stick', 'sphere', 'carton'
            surface, bool, display SAS
            opacity, float, opacity of surface, range 0.0-1.0
        Return:
        ----
            viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
        """

        assert style in ("line", "stick", "sphere", "carton")
        mblock = Chem.MolToMolBlock(mol)
        viewer = py3Dmol.view(width=size[0], height=size[1])
        viewer.addModel(mblock, "mol")
        viewer.setStyle({style: {}})
        if surface:
            viewer.addSurface(py3Dmol.SAS, {"opacity": opacity})
        viewer.zoomTo()
        return viewer
