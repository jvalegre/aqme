"""
Parameters
----------

General
+++++++

   files : str or list of str, default=None
     Input files. Formats accepted: SDF. Also, lists can be used
     (i.e. [FILE1.sdf, FILE2.sdf] or \*.sdf).
   method : str, default=None
     FAMEX backend used for geometry optimization. If omitted, the default
     xTB/tblite backend is used.
     Current options: 'xtb', 'aimnet2', 'mace', 'orb', 'so3lr', 'uma'
   w_dir_main : str, default=os.getcwd()
     Working directory
   destination : str, default=None
     Directory to create the output file(s)
   varfile : str, default=None
     Option to parse the variables using a yaml file (specify the filename)
   charge : int, default=None
     Charge for FAMEX calculations. If not set, read from SDF 'Real charge'
     property; if missing, defaults to 0.
   mult : int, default=None
     Multiplicity for FAMEX calculations. If not set, read from SDF 'Mult'
     property; if missing, defaults to 1.
   ewin_cmin : float, default=5.0
     Energy window in kcal/mol to discard conformers after optimization.
   initial_energy_threshold : float, default=0.0001
     Energy difference in kcal/mol between unique conformers (first filter).
   energy_threshold : float, default=0.25
     Energy difference in kcal/mol between unique conformers (second filter).
   rms_threshold : float, default=0.25
     RMS difference threshold for the second filter.
   opt_fmax : float, default=0.05
     Convergence criterion (max force, eV/Å) for the BFGS optimizer in FAMEX.
   opt_steps : int, default=1000
     Maximum number of BFGS steps.
   constraints_dist : list of lists, default=[]
     Distance constraints for FAMEX FixInternals as [[AT1,AT2,DIST], ...].
   constraints_angle : list of lists, default=[]
     Angle constraints for FAMEX FixInternals as [[AT1,AT2,AT3,ANGLE], ...].
     Angles are specified in degrees.
   constraints_dihedral : list of lists, default=[]
     Dihedral constraints for FAMEX FixInternals as [[AT1,AT2,AT3,AT4,DIHEDRAL], ...].
     Dihedrals are specified in degrees.
   target : str, default='minima'
     Optimization target for FAMEX. Accepted values: 'minima' and 'ts'.
   freq : bool, default=False
     If True, calculates vibrational frequencies for conformers surviving filters.
   prefix : str, default=''
     Prefix added to all output names.
   suffix : str, default=''
     Suffix added to all output names.
"""

#####################################################.
#          This file stores the CMIN class          #
#        Conformer refinement via FAMEX backend       #
#####################################################.

import os
import sys
import time
import ast
import contextlib
import threading
import numpy as np
from pathlib import Path
import concurrent.futures

from rdkit.Chem import AllChem as Chem
from rdkit.Chem.PropertyMol import PropertyMol
from progress.bar import IncrementalBar
from rdkit.Geometry import Point3D

from aqme.utils import (
    load_variables,
    mol_from_sdf_or_mol_or_mol2,
    add_prefix_suffix,
    check_dependencies,
    set_destination,
)
from aqme.csearch.utils import _translate_constraint_indices
from aqme.filter import conformer_filters, cluster_conformers

SUPPORTED_FAMEX_BACKENDS = {"xtb", "tblite", "aimnet2", "mace", "orb", "so3lr", "uma"}
EV_TO_KCAL = 23.0609  # 1 eV = 23.0609 kcal/mol


def normalize_cmin_backend(args):
    """Normalize the requested CMIN backend and keep xtb mapped to tblite."""
    backend = (
        getattr(args, "method", None)
        or getattr(args, "model", None)
        or getattr(args, "program", None)
        or "xtb"
    )
    backend = str(backend).lower()

    if backend not in SUPPORTED_FAMEX_BACKENDS:
        raise ValueError(
            f"Method '{backend}' is not supported. "
            f"Choose one of: {', '.join(sorted(SUPPORTED_FAMEX_BACKENDS))}"
        )

    if backend == "xtb":
        backend = "tblite"

    args.method = backend
    args.model = backend
    args.program = backend
    return args


class cmin:
    _tblite_lock = threading.Lock()

    """
    Conformer refinement using a FAMEX backend (default: xTB/tblite).

    Reads conformers from SDF files, optimizes them using FAMEX, applies 
    energy/RMSD filters, and writes the results to SDF format.

    Parameters
    ----------
    kwargs : keyword arguments
        Refer to the module docstring for the full list of options (excluding 'method').
    """

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def __init__(self, **kwargs):
        start_time = time.time()

        self.args = load_variables(kwargs, "cmin")
        # check_dependencies(self)
        self._validate_target()

        # Validate backend choice
        self._validate_program()

        # Check we have input files
        if not self.args.files:
            self.args.log.write(
                "\nx  No files were found! Make sure you use quotation marks "
                "if you are using * (i.e. --files \"*.sdf\")"
            )
            self.args.log.finalize()
            sys.exit()

        # Only SDF input is supported now
        file_format = Path(self.args.files[0]).suffix.lower().lstrip(".")
        if file_format != "sdf":
            self.args.log.write(
                f"\nx  Only SDF input is supported in this version of CMIN "
                f"(got '.{file_format}'). Convert your files to SDF first."
            )
            self.args.log.finalize()
            sys.exit()

        bar = IncrementalBar(
            "\no  Number of finished jobs from CMIN", max=len(self.args.files)
        )

        for sdf_file in self.args.files:
            self.mols, self.name = self._load_sdf(sdf_file)
            self.name = add_prefix_suffix(self.name, self.args)
            self.args.log.write(f"\n\n   ----- {self.name} -----")

            self._setup_output_directories()
            self.compute_cmin(sdf_file)

            bar.next()

        bar.finish()

        elapsed = round(time.time() - start_time, 2)
        self.args.log.write(f"\nTime CMIN: {elapsed} seconds\n")
        self.args.log.finalize()

        # Return to the original directory (important for Jupyter)
        os.chdir(self.args.initial_dir)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_program(self):
        """Check that the requested FAMEX backend is importable."""
        try:
            self.args = normalize_cmin_backend(self.args)
        except ValueError as exc:
            self.args.log.write(f"\nx  {exc}")
            self.args.log.finalize()
            sys.exit()

        try:
            import famex 
        except ImportError:
            self.args.log.write(
                "\nx  FAMEX is not installed. Install it with: pip install famex"
            )
            self.args.log.finalize()
            sys.exit()

    def _validate_target(self):
        """Validate and normalize the requested FAMEX optimization target."""
        target = getattr(self.args, "target", None)
        if target in [None, ""]:
            target = "minima"

        target = str(target).lower()
        if target not in {"minima", "ts"}:
            self.args.log.write(
                "\nx  Invalid target value. Choose 'minima' or 'ts'."
            )
            self.args.log.finalize()
            sys.exit()

        self.args.target = target

    def _get_famex_target(self):
        """Return the FAMEX Explorer target for the current optimization mode."""
        return "ts" if getattr(self.args, "target", "minima") == "ts" else "minima"

    def _log_famex_citation(self):
        """Log the recommended FAMEX citation the first time it is used."""
        if getattr(self, "_famex_citation_logged", False):
            return

        self.args.log.write(
            "\n   Please cite FAMEX as:\n"
            "   FAMEX Development Team, FAMEX: Fast Mechanistic Explorer (2026).\n"
            "   Available at https://github.com/rlaplaza-lab/famex"
        )
        self._famex_citation_logged = True

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    def _load_sdf(self, sdf_file):
        """Load molecules from an SDF file.

        Returns
        -------
        tuple : (list[Mol], str)
            Molecule objects and the base name (stem of the file).
        """
        try:
            mols = mol_from_sdf_or_mol_or_mol2(sdf_file, "cmin", self.args)
        except OSError:
            path = Path(self.args.initial_dir) / sdf_file
            mols = mol_from_sdf_or_mol_or_mol2(str(path), "cmin", self.args)

        name = Path(sdf_file).stem
        return mols, name

    def _setup_output_directories(self):
        """Create CMIN output folders and open SDWriters."""
        self.cmin_folder = set_destination(self, "CMIN")
        self.cmin_folder.mkdir(exist_ok=True, parents=True)
        self.cmin_folder.joinpath("All_confs").mkdir(exist_ok=True, parents=True)

        all_confs_path = self.cmin_folder / "All_confs" / (
            f"{self.name}_all_confs{self.args.output}"
        )
        self.sdwriterall = Chem.SDWriter(str(all_confs_path))

        filtered_path = self.cmin_folder / (
            f"{self.name}{self.args.output}"
        )
        self.sdwriter = Chem.SDWriter(str(filtered_path))

    # ------------------------------------------------------------------
    # Charge / multiplicity
    # ------------------------------------------------------------------

    def _read_charge_mult_from_sdf(self, sdf_file):
        """Read 'Real charge' and 'Mult' properties written by CSEARCH."""
        charge, mult = None, None
        try:
            def _parse_int_like(raw_value):
                """Parse numeric strings like '1', '1.0', '+1', '-1.000' into int."""
                token = str(raw_value).strip().split()[0]
                return int(round(float(token)))

            with open(sdf_file, "r") as fh:
                lines = fh.readlines()
            charge_found = mult_found = False
            for i, line in enumerate(lines):
                # Support common SDF header variants, e.g.:
                # >  <Real charge>  (1)
                # > <Real charge>
                # >  <Mult>  (1)
                if "<Real charge>" in line:
                    charge = _parse_int_like(lines[i + 1])
                    charge_found = True
                if "<Mult>" in line:
                    mult = _parse_int_like(lines[i + 1])
                    mult_found = True
                if charge_found and mult_found:
                    break
        except Exception:
            pass
        return charge, mult

    def _determine_charge_mult(self, sdf_file):
        """Return (charge, mult) to pass to FAMEX.

        Priority: user arg > SDF property > default (0, 1).
        """
        file_charge, file_mult = self._read_charge_mult_from_sdf(sdf_file)

        charge = self.args.charge if self.args.charge is not None else (
            file_charge if file_charge is not None else 0
        )
        mult = self.args.mult if self.args.mult is not None else (
            file_mult if file_mult is not None else 1
        )

        if self.args.charge is None and file_charge is None:
            self.args.log.write(
                "\n   No charge found - defaulting to 0. "
                "Set charge= or add 'Real charge' to SDF."
            )
        if self.args.mult is None and file_mult is None:
            self.args.log.write(
                "\n   No multiplicity found - defaulting to 1. "
                "Set mult= or add 'Mult' to SDF."
            )

        return int(charge), int(mult)
    # ------------------------------------------------------------------
    # Core FAMEX optimisation
    # ------------------------------------------------------------------

    def _mol_to_ase_atoms(self, mol, charge, mult):
        """Convert an RDKit Mol (single conformer) to an ASE Atoms object."""
        try:
            from ase import Atoms
        except ImportError:
            self.args.log.write(
                "\nx  ASE is not installed. Install it with: pip install ase"
            )
            self.args.log.finalize()
            sys.exit()

        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        positions = mol.GetConformer().GetPositions()  # Å
        atoms = Atoms(symbols=symbols, positions=positions)
        atoms.info["charge"] = int(charge)
        atoms.info["spin"] = int(mult)
        return atoms

    def _update_mol_from_ase(self, mol, ase_atoms):
        """Write optimised ASE coordinates back into the RDKit conformer (in-place)."""
        positions = ase_atoms.get_positions()  # Å
        conf = mol.GetConformer()
        for i, (x, y, z) in enumerate(positions):
            conf.SetAtomPosition(i, Point3D(x, y, z))

    def _get_atom_map_to_idx(self, mol):
        """Return atom-map-number to RDKit-index mapping, if present."""
        map_to_idx = {}
        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()
            if atom_map:
                map_to_idx[atom_map] = atom.GetIdx()
        return map_to_idx

    def _prepare_constraint_indices(self, mol, constraints, n_indices):
        """Translate atom-map numbers when available, otherwise keep direct indices."""
        constraints = self._as_constraint_list(constraints)
        if len(constraints) == 0:
            return []

        map_to_idx = self._get_atom_map_to_idx(mol)
        if map_to_idx:
            return _translate_constraint_indices(
                constraints,
                map_to_idx,
                n_indices=n_indices,
                log=self.args.log,
            )

        return [list(constraint) for constraint in constraints]

    def _as_constraint_list(self, constraints):
        """Return constraints as a regular list without numpy truth-value ambiguity."""
        if constraints is None:
            return []
        if isinstance(constraints, str):
            return ast.literal_eval(constraints)
        if isinstance(constraints, np.ndarray):
            return constraints.tolist()
        return constraints

    def _add_fixinternals_specs(
        self, specs, constraints, n_indices, famex_type, arg_name, n_atoms
    ):
        """Append validated FAMEX FixInternals specs from AQME nested constraints."""
        for constraint in constraints:
            if not isinstance(constraint, (list, tuple, np.ndarray)):
                raise ValueError(
                    f"\nx  {arg_name} must be a list of lists. Invalid entry: {constraint}"
                )

            if len(constraint) != n_indices + 1:
                raise ValueError(
                    f"\nx  {arg_name} entries must contain {n_indices} atom index(es) "
                    f"and one target value. Invalid entry: {constraint}"
                )

            try:
                atom_indices = [int(idx) for idx in constraint[:n_indices]]
                target_value = float(constraint[n_indices])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"\nx  Invalid numeric value in {arg_name}: {constraint}"
                ) from exc

            for atom_idx in atom_indices:
                if atom_idx < 0 or atom_idx >= n_atoms:
                    raise ValueError(
                        f"\nx  Atom index {atom_idx} in {arg_name} is out of range "
                        f"for a molecule with {n_atoms} atoms."
                    )

            indices_str = ",".join(str(idx) for idx in atom_indices)
            specs.append(f"fixinternals_{famex_type} {indices_str} value={target_value}")

    def _build_famex_constraints(self, mol):
        """Builds FAMEX constraints using local copies to prevent accumulation."""
        
        # Get user-defined constraints from the command line/args
        user_dist = self._as_constraint_list(getattr(self.args, "constraints_dist", []))
        
        # Combine with automatic haptic constraints (returned as a new list)
        constraints_dist = self._apply_automatic_metal_constraints(mol, user_dist)
        
        constraints_angle = self._as_constraint_list(getattr(self.args, "constraints_angle", []))
        constraints_dihedral = self._as_constraint_list(getattr(self.args, "constraints_dihedral", []))

        if not any(len(c) > 0 for c in [constraints_dist, constraints_angle, constraints_dihedral]):
            return None

        n_atoms = mol.GetNumAtoms()
        # Translate the combined list of map numbers to RDKit indices
        constraints_dist = self._prepare_constraint_indices(mol, constraints_dist, 2)
        constraints_angle = self._prepare_constraint_indices(mol, constraints_angle, 3)
        constraints_dihedral = self._prepare_constraint_indices(mol, constraints_dihedral, 4)

        specs = []
        self._add_fixinternals_specs(specs, constraints_dist, 2, "bond", "constraints_dist", n_atoms)
        self._add_fixinternals_specs(specs, constraints_angle, 3, "angle", "constraints_angle", n_atoms)
        self._add_fixinternals_specs(specs, constraints_dihedral, 4, "dihedral", "constraints_dihedral", n_atoms)

        constraints_str = "; ".join(specs)
        self.args.log.write(f"\n   Applying FAMEX FixInternals constraints: {constraints_str}")
        return constraints_str

    def _optimize_with_famex(self, mol, conf_name, charge, mult, constraints=None):
        """Run a single FAMEX local minimisation.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            Molecule with exactly one conformer.
        conf_name : str
            Label for logging.
        charge : int
            Total molecular charge.
        mult : int
            Spin multiplicity.

        Returns
        -------
        tuple : (mol, energy_kcal, success)
        """
        import famex

        self.args.log.write(f"\no  FAMEX optimisation [{self.args.program}] ({conf_name})")

        try:
            ase_atoms = self._mol_to_ase_atoms(mol, charge, mult)
            target = self._get_famex_target()

            fmax = getattr(self.args, "opt_fmax", 0.05)
            steps = getattr(self.args, "opt_steps", 1000)

            if self.args.program == "tblite":
                with cmin._tblite_lock:
                    # Redirect to os.devnull to avoid buffering large strings.
                    with open(os.devnull, "w", encoding="utf-8") as devnull:
                        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                            explorer = famex.Explorer(
                                atoms=ase_atoms,
                                backend=self.args.program,
                                target=target,
                                strategy="local",
                                default_charge=charge,
                                default_spin=mult,
                                constraints=constraints,
                                verbose=0
                            )
                            result = explorer.run(fmax=fmax, steps=steps)
            else:
                explorer = famex.Explorer(
                    atoms=ase_atoms,
                    backend=self.args.program,
                    target=target,
                    strategy="local",
                    default_charge=charge,
                    default_spin=mult,
                    constraints=constraints,
                    verbose=0
                )
                result = explorer.run(fmax=fmax, steps=steps)

            optimised_atoms = result["optimized_atoms"]
            energy_ev = optimised_atoms.get_potential_energy()  # eV
            energy_kcal = energy_ev * EV_TO_KCAL

            self._update_mol_from_ase(mol, optimised_atoms)
            self.args.log.write(
                f"\n   Converged. E = {energy_ev:.6f} eV  ({energy_kcal:.4f} kcal/mol)"
            )
            return mol, energy_kcal, True

        except Exception as exc:
            self.args.log.write(
                f"\nx  FAMEX optimisation failed for {conf_name}: {exc}"
            )
            return mol, 0.0, False

    def _calculate_frequencies(self, mol, conf_label, charge, mult):
        """Calculate vibrational frequencies for a given conformer using the FAMEX analysis module.
        
        Args:
            mol (rdkit.Chem.PropertyMol): Conformer molecule object.
            conf_label (str): Label used for naming the output file.
            charge (int): Molecular charge.
            mult (int): Spin multiplicity.
        """
        import famex
        from famex.analysis.frequency import FrequencyAnalysis

        self.args.log.write(f"\no  FAMEX frequency calculation [{self.args.program}] ({conf_label})")

        try:
            ase_atoms = self._mol_to_ase_atoms(mol, charge, mult)
            target = self._get_famex_target()

            # Suppress internal FAMEX logs during calculator initialization steps
            with open(os.devnull, "w", encoding="utf-8") as devnull:
                with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                    if self.args.program == "tblite":
                        with cmin._tblite_lock:
                            explorer = famex.Explorer(
                                atoms=ase_atoms,
                                backend=self.args.program,
                                target=target,
                                strategy="local",
                                default_charge=charge,
                                default_spin=mult,
                                verbose=0
                            )
                            # Run 1 step to ensure the result dictionary populates "optimized_atoms" with the calculator
                            res = explorer.run(steps=1)
                            active_atoms = res["optimized_atoms"]
                    else:
                        explorer = famex.Explorer(
                            atoms=ase_atoms,
                            backend=self.args.program,
                            target=target,
                            strategy="local",
                            default_charge=charge,
                            default_spin=mult,
                            verbose=0
                        )
                        res = explorer.run(steps=1)
                        active_atoms = res["optimized_atoms"]

            # Perform frequency and normal mode analysis
            calc = active_atoms.calc
            freq_analyzer = FrequencyAnalysis(active_atoms, calc)
            frequencies = freq_analyzer.get_frequencies()
            normal_modes = freq_analyzer.get_normal_modes() # Returns (3N x M) matrix
            freq_values = [float(freq) for freq in np.asarray(frequencies).ravel()]
            imaginary_analysis = []
            for i, freq in enumerate(freq_values):
                if freq < 0.0:
                    mode_vectors = normal_modes[:, i].reshape(-1, 3)
                    displacements = np.linalg.norm(mode_vectors, axis=1)
                    imaginary_analysis.append(
                        {
                            "frequency": float(freq),
                            "top_atoms": [
                                {
                                    "atom_idx": int(atom_idx),
                                    "symbol": active_atoms.get_chemical_symbols()[atom_idx],
                                    "displacement": float(displacements[atom_idx]),
                                }
                                for atom_idx in np.argsort(displacements)[::-1][:5]
                            ],
                        }
                    )
            return {
                "conf_label": conf_label,
                "frequencies": freq_values,
                "n_negative": int(sum(1 for freq in freq_values if freq < 0.0)),
                "imaginary_analysis": imaginary_analysis,
            }

        except Exception as exc:
            self.args.log.write(
                f"\nx  FAMEX frequency calculation failed for {conf_label}: {exc}"
            )
            return None

    def _write_frequencies_report(self, freq_results):
        """Write the combined frequency report to a single file."""
        freq_file = self.cmin_folder / "frecuencies.dat"
        target = getattr(self.args, "target", "minima")

        one_negative = [item["conf_label"] for item in freq_results if item["n_negative"] == 1]
        two_or_more = [
            item["conf_label"] for item in freq_results if item["n_negative"] >= 2
        ]
        none = [item["conf_label"] for item in freq_results if item["n_negative"] == 0]

        with open(freq_file, "w", encoding="utf-8") as fh:
            fh.write(f"CMIN frequencies report for {self.name}\n")
            fh.write(f"Target: {target}\n")
            fh.write(f"Program: {self.args.program}\n\n")

            if target == "ts":
                fh.write("TS frequency summary\n")
                fh.write(
                    f"One negative frequency: {', '.join(one_negative) if one_negative else 'none'}\n"
                )
                fh.write(
                    f"Two or more negative frequencies: {', '.join(two_or_more) if two_or_more else 'none'}\n"
                )
                fh.write(
                    f"No negative frequencies: {', '.join(none) if none else 'none'}\n\n"
                )
                self.args.log.write("\no  TS frequency summary")
                self.args.log.write(
                    f"\n   One negative frequency: "
                    f"{', '.join(one_negative) if one_negative else 'none'}"
                )
                self.args.log.write(
                    f"\n   Two or more negative frequencies: "
                    f"{', '.join(two_or_more) if two_or_more else 'none'}"
                )
                self.args.log.write(
                    f"\n   No negative frequencies: "
                    f"{', '.join(none) if none else 'none'}"
                )

            if not freq_results:
                fh.write("No frequency calculations were performed.\n")
                self.args.log.write(
                    f"\no  Frequencies report saved to {freq_file}"
                )
                return freq_file

            for item in freq_results:
                fh.write(f"Conformer: {item['conf_label']}\n")
                fh.write(f"Negative frequencies: {item['n_negative']}\n")
                fh.write(f"Frequencies (cm^-1): {item['frequencies']}\n")
                if item["imaginary_analysis"]:
                    for analysis in item["imaginary_analysis"]:
                        fh.write(
                            f"  Imaginary frequency: {analysis['frequency']:.2f} cm^-1\n"
                        )
                        fh.write("  Top moving atoms:\n")
                        for atom in analysis["top_atoms"]:
                            fh.write(
                                f"    Atom {atom['atom_idx']} ({atom['symbol']}) : "
                                f"{atom['displacement']:.4f}\n"
                            )
                fh.write("\n")

        self.args.log.write(f"\no  Frequencies report saved to {freq_file}")
        return freq_file

    # ------------------------------------------------------------------
    # Main compute loop
    # ------------------------------------------------------------------

    def compute_cmin(self, sdf_file):
        """Optimise all conformers from *sdf_file* and write results."""
        charge, mult = self._determine_charge_mult(sdf_file)
        self.args.log.write(
            f"\n   charge={charge}  mult={mult}  program={self.args.program}"
        )

        cenergy, outmols = [], []
        nprocs = getattr(self.args, "nprocs", None)
        if not nprocs:
            nprocs = 1
        else:
            nprocs = max(1, int(nprocs))

        # Prepare the list of tasks (skipping None values)
        tasks = []
        for i, mol in enumerate(self.mols):
            if mol is not None:
                conf_name = f"{self.name}_conf_{i}"
                tasks.append((mol, conf_name, charge, mult))
        if not tasks:
            self.args.log.write(
                f"\nx  No valid initial conformers found for {self.name}."
            )
            return

        self._log_famex_citation()
        constraints = self._build_famex_constraints(tasks[0][0])
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
            futures_list = [
                executor.submit(
                    self._optimize_with_famex,
                    task[0],
                    task[1],
                    task[2],
                    task[3],
                    constraints,
                )
                for task in tasks
            ]

            # Gather results as they complete
            for future in concurrent.futures.as_completed(futures_list):
                try:
                    result = future.result()
                    if result is not None:
                        mol, energy, ok = result
                        if ok:
                            pmol = PropertyMol(mol)
                            outmols.append(pmol)
                            cenergy.append(energy)
                except Exception as exc:
                    self.args.log.write(f"\nx  A parallel worker crashed: {exc}")

        if not cenergy:
            self.args.log.write(
                f"\nx  No conformers were successfully optimised for {self.name}."
            )
            self.sdwriterall.close()
            self.sdwriter.close()
            return

        # Sort by energy
        cids = list(range(len(outmols)))
        sorted_cids = sorted(cids, key=lambda cid: cenergy[cid])

        # Add properties
        for cid in sorted_cids:
            outmols[cid].SetProp("Energy", str(cenergy[cid]))
            outmols[cid].SetProp("Real charge", str(charge))
            outmols[cid].SetProp("Mult", str(mult))

        # Write all conformers before filtering
        for cid in sorted_cids:
            self.sdwriterall.write(outmols[cid])
        self.sdwriterall.close()

        # Apply energy + RMSD filters
        self.args.log.write(
            f"\o  Applying filters after {self.args.program} minimisation"
        )
        selected_cids = conformer_filters(self, sorted_cids, cenergy, outmols)

        # Write filtered conformers
        filtered_mols = [outmols[cid] for cid in selected_cids]

        freq_results = []
        # Calculate frequencies ONLY for structures that survived the filtering criteria
        if getattr(self.args, "freq", False):
            for i, mol in enumerate(filtered_mols):
                conf_label = f"{self.name}_filtered_conf_{i}"
                freq_result = self._calculate_frequencies(mol, conf_label, charge, mult)
                if freq_result is not None:
                    freq_results.append(freq_result)
            self._write_frequencies_report(freq_results)

        for mol in filtered_mols:
            self.sdwriter.write(mol)
        self.sdwriter.close()

        final_count = len(filtered_mols)
        if (
            self.args.auto_cluster
            and hasattr(self.args, "sample")
            and len(filtered_mols) > int(self.args.sample)
        ):
            self.args.log.write(
                f"\no  Applying final Butina clustering after minimisation "
                f"filters ({self.name})"
            )
            output_sdf = self.cmin_folder / (self.name + self.args.output)
            clustered_mols = cluster_conformers(
                self, filtered_mols, "rdkit", output_sdf, self.name, int(self.args.sample)
            )
            final_count = len(clustered_mols)

        self.args.log.write(
            f"\no  {final_count} conformer(s) written to "
            f"{self.cmin_folder / (self.name + self.args.output)}"
        )


# ------EXTRA THINGS------
# List of transition metals for detection
    TRANSITION_METALS = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
    }

    def _apply_automatic_metal_constraints(self, mol, current_dist_constraints):
        """
        Creates haptic constraints based on map numbers (1000, 2000, etc.)
        without modifying the persistent self.args object.
        """
        # 1. Ensure map numbers are set from SDF properties if missing
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                if atom.HasProp('molAtomMapNumber'):
                    atom.SetAtomMapNum(int(atom.GetProp('molAtomMapNumber')))

        # 2. Create a local copy of the constraints for THIS molecule only
        total_dist = [list(c) for c in current_dist_constraints]
        
        haptic_groups = {}
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num >= 1000:
                group_prefix = (map_num // 1000) * 1000
                if group_prefix not in haptic_groups:
                    haptic_groups[group_prefix] = []
                haptic_groups[group_prefix].append(atom)

        if not haptic_groups:
            return total_dist

        # 3. Process each ring group independently
        for prefix, atoms in haptic_groups.items():
            ring_size = len(atoms)
            # Standard distances: 2.1 A (Cp) or 2.25 A (Benzene)
            dist_val = 2.1 if ring_size == 5 else 2.25
            
            ring_pos = np.mean([mol.GetConformer().GetAtomPosition(a.GetIdx()) for a in atoms], axis=0)
            closest_metal = None
            min_dist = float('inf')
            
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in self.TRANSITION_METALS:
                    metal_pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                    dist = np.linalg.norm(ring_pos - metal_pos)
                    if dist < min_dist:
                        min_dist = dist
                        closest_metal = atom

            if closest_metal is None:
                continue

            # Ensure the metal has a map number below 1000
            metal_map = closest_metal.GetAtomMapNum()
            if metal_map == 0 or metal_map >= 1000:
                existing_maps = [a.GetAtomMapNum() for a in mol.GetAtoms() if 0 < a.GetAtomMapNum() < 1000]
                metal_map = max(existing_maps) + 1 if existing_maps else 1
                closest_metal.SetAtomMapNum(metal_map)

            # Target atoms 1, 3, and 5 of the ring
            target_maps = [prefix, prefix + 2, prefix + 4]
            for t_map in target_maps:
                if any(a.GetAtomMapNum() == t_map for a in atoms):
                    total_dist.append([t_map, metal_map, dist_val])

            self.args.log.write(f"\n   o  Group {prefix}: Applying {ring_size}-member ring constraints ({dist_val} A).")

        return total_dist
