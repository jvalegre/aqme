#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   CMIN module                      #
######################################################.

import os
import glob
import sys
import types
from pathlib import Path

import pytest
import rdkit
import shutil

from aqme.cmin import cmin

# Mapeo de rutas relativas a la ubicación del archivo del test
tests_dir = os.path.dirname(os.path.abspath(__file__))
w_dir_main = os.path.dirname(tests_dir) # Raíz del repositorio (aqme)

cmin_methods_dir = os.path.join(tests_dir, "cmin_methods")
cmin_empty_dir = os.path.join(tests_dir, "cmin_empty")
cmin_ts_dir = os.path.join(tests_dir, "cmin_TS")

if not os.path.exists(cmin_methods_dir):
    os.mkdir(cmin_methods_dir)


def _repo_path(*parts):
    return Path(w_dir_main).joinpath(*parts)


def _copy_fixture_to_local(source_relative, target_relative):
    source = Path(w_dir_main) / source_relative
    target = _repo_path(*target_relative)
    target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return target


def _write_sdf_without_properties(source_relative, target_relative):
    source = Path(w_dir_main) / source_relative
    lines = source.read_text(encoding="utf-8").splitlines()
    end_idx = next(i for i, line in enumerate(lines) if line.strip() == "M  END")
    target = _repo_path(*target_relative)
    target.write_text("\n".join(lines[: end_idx + 1]) + "\n$$$$\n", encoding="utf-8")
    return target


def _read_first_mol_properties(sdf_path):
    mols = rdkit.Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    return next(mol for mol in mols if mol is not None)


class _FakeAtoms:
    def __init__(self, mol, energy=0.0):
        self._positions = mol.GetConformer().GetPositions()
        self._symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        self._energy = energy
        self.calc = object()

    def get_positions(self):
        return self._positions

    def get_potential_energy(self):
        return self._energy

    def get_chemical_symbols(self):
        return self._symbols


def _make_fake_famex(monkeypatch, explorer_calls):
    fake_famex = types.ModuleType("famex")
    fake_analysis = types.ModuleType("famex.analysis")
    fake_frequency = types.ModuleType("famex.analysis.frequency")
    fake_frequency.FrequencyAnalysis = object

    class FakeExplorer:
        def __init__(
            self,
            atoms,
            backend,
            target,
            strategy,
            default_charge,
            default_spin,
            constraints=None,
            verbose=0,
        ):
            explorer_calls.append(
                {
                    "backend": backend,
                    "target": target,
                    "charge": default_charge,
                    "spin": default_spin,
                    "constraints": constraints,
                    "strategy": strategy,
                }
            )
            self.atoms = atoms

        def run(self, fmax=None, steps=None):
            return {"optimized_atoms": self.atoms}

    fake_famex.Explorer = FakeExplorer
    fake_analysis.frequency = fake_frequency
    fake_famex.analysis = fake_analysis
    monkeypatch.setitem(sys.modules, "famex", fake_famex)
    monkeypatch.setitem(sys.modules, "famex.analysis", fake_analysis)
    monkeypatch.setitem(sys.modules, "famex.analysis.frequency", fake_frequency)
    return fake_famex


def _fake_optimize_success(self, mol, conf_name, charge, mult, constraints=None):
    return mol, 0.0, True


def _fake_frequency_result(conf_label, symbols):
    return {
        "conf_label": conf_label,
        "frequencies": [-123.4, 25.0, 80.1],
        "n_negative": 1,
        "imaginary_analysis": [
            {
                "frequency": -123.4,
                "top_atoms": [
                    {"atom_idx": 0, "symbol": symbols[0], "displacement": 1.0},
                    {"atom_idx": 4, "symbol": symbols[4], "displacement": 0.8},
                    {"atom_idx": 5, "symbol": symbols[5], "displacement": 0.7},
                ],
            }
        ],
    }


# tests for invalid input handling
@pytest.mark.parametrize(
    "program",
    ["invalid_method"],
)
def test_cmin_invalid_backend_raises_system_exit(monkeypatch, capsys, program):
    monkeypatch.chdir(cmin_methods_dir)
    sdf_path = _repo_path("tests", "cmin_methods", "pentane_rdkit_methods.sdf")

    with pytest.raises(SystemExit):
        cmin(program=program, files=str(sdf_path))

    captured = capsys.readouterr()
    assert "not supported" in captured.out


@pytest.mark.parametrize("extension", [".mol2", ".pdb"])
def test_cmin_rejects_non_sdf_extensions(monkeypatch, capsys, extension):
    monkeypatch.chdir(cmin_methods_dir)
    input_file = _repo_path("tests", "cmin_methods", f"input{extension}")
    input_file.write_text("dummy", encoding="utf-8")
    _make_fake_famex(monkeypatch, [])

    with pytest.raises(SystemExit):
        cmin(program="xtb", files=str(input_file))

    captured = capsys.readouterr()
    assert "Only SDF input is supported" in captured.out
    input_file.unlink(missing_ok=True)


@pytest.mark.parametrize("target", ["unknown", "minimum", "ts-like"])
def test_cmin_rejects_invalid_target(monkeypatch, capsys, target):
    monkeypatch.chdir(cmin_methods_dir)
    sdf_path = _repo_path("tests", "cmin_methods", "pentane_rdkit_methods.sdf")

    with pytest.raises(SystemExit):
        cmin(program="xtb", files=str(sdf_path), target=target)

    captured = capsys.readouterr()
    assert "Invalid target value" in captured.out


# tests for charge and multiplicity handling
def test_cmin_reads_charge_and_mult_from_sdf(monkeypatch):
    monkeypatch.chdir(cmin_ts_dir)
    sdf_path = _repo_path("tests", "cmin_TS", "methyl_Cl_Br.sdf")
    _make_fake_famex(monkeypatch, [])

    monkeypatch.setattr("aqme.cmin.cmin._optimize_with_famex", _fake_optimize_success)
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="xtb", files=str(sdf_path))

    output_file = _repo_path("tests", "cmin_TS", "CMIN", "methyl_Cl_Br.sdf")
    assert output_file.exists()
    out_mol = _read_first_mol_properties(output_file)

    assert out_mol.GetProp("Real charge") == "-1"
    assert out_mol.GetProp("Mult") == "1"
    output_file.unlink(missing_ok=True)


def test_cmin_explicit_charge_and_mult_override_sdf(monkeypatch):
    monkeypatch.chdir(cmin_ts_dir)
    sdf_path = _repo_path("tests", "cmin_TS", "methyl_Cl_Br.sdf")
    _make_fake_famex(monkeypatch, [])

    monkeypatch.setattr("aqme.cmin.cmin._optimize_with_famex", _fake_optimize_success)
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="xtb", files=str(sdf_path), charge=2, mult=3)

    output_file = _repo_path("tests", "cmin_TS", "CMIN", "methyl_Cl_Br.sdf")
    assert output_file.exists()
    out_mol = _read_first_mol_properties(output_file)

    assert out_mol.GetProp("Real charge") == "2"
    assert out_mol.GetProp("Mult") == "3"
    output_file.unlink(missing_ok=True)


def test_cmin_defaults_charge_and_mult_when_missing(monkeypatch):
    monkeypatch.chdir(cmin_methods_dir)
    sdf_path = _write_sdf_without_properties(
        "tests/cmin_methods/pentane_rdkit_methods.sdf",
        ("tests", "cmin_methods", "pentane_no_charge_mult.sdf"),
    )
    _make_fake_famex(monkeypatch, [])

    monkeypatch.setattr("aqme.cmin.cmin._optimize_with_famex", _fake_optimize_success)
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="xtb", files=str(sdf_path))

    output_file = _repo_path("tests", "cmin_methods", "CMIN", "pentane_no_charge_mult.sdf")
    assert output_file.exists()
    out_mol = _read_first_mol_properties(output_file)

    assert out_mol.GetProp("Real charge") == "0"
    assert out_mol.GetProp("Mult") == "1"
    output_file.unlink(missing_ok=True)
    sdf_path.unlink(missing_ok=True)


# tests for target handling and frequency output
@pytest.mark.parametrize(
    "target, expected_famex_target",
    [("minima", "minima"), ("ts", "ts")],
)
def test_cmin_forwards_target_to_famex(monkeypatch, target, expected_famex_target):
    monkeypatch.chdir(cmin_ts_dir)
    sdf_path = _repo_path("tests", "cmin_TS", "methyl_Cl_Br.sdf")

    explorer_calls = []
    _make_fake_famex(monkeypatch, explorer_calls)

    def fake_mol_to_ase_atoms(self, mol, charge, mult):
        return _FakeAtoms(mol)

    monkeypatch.setattr("aqme.cmin.cmin._mol_to_ase_atoms", fake_mol_to_ase_atoms)
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="tblite", files=str(sdf_path), target=target)

    assert explorer_calls
    assert explorer_calls[0]["target"] == expected_famex_target


def test_cmin_writes_frequencies_report_when_freq_enabled(monkeypatch):
    monkeypatch.chdir(cmin_ts_dir)
    sdf_path = _repo_path("tests", "cmin_TS", "methyl_Cl_Br.sdf")
    _make_fake_famex(monkeypatch, [])

    monkeypatch.setattr("aqme.cmin.cmin._optimize_with_famex", _fake_optimize_success)
    monkeypatch.setattr(
        "aqme.cmin.cmin._calculate_frequencies",
        lambda self, mol, conf_label, charge, mult: _fake_frequency_result(
            conf_label, [atom.GetSymbol() for atom in mol.GetAtoms()]
        ),
    )
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="xtb", files=str(sdf_path), freq=True)

    freq_file = _repo_path("tests", "cmin_TS", "CMIN", "frecuencies.dat")
    assert freq_file.exists()
    assert "Negative frequencies: 1" in freq_file.read_text(encoding="utf-8")
    freq_file.unlink(missing_ok=True)


def test_cmin_ts_report_marks_negative_frequency_and_top_atoms(monkeypatch):
    monkeypatch.chdir(cmin_ts_dir)
    sdf_path = _repo_path("tests", "cmin_TS", "methyl_Cl_Br.sdf")
    _make_fake_famex(monkeypatch, [])

    monkeypatch.setattr("aqme.cmin.cmin._optimize_with_famex", _fake_optimize_success)
    monkeypatch.setattr(
        "aqme.cmin.cmin._calculate_frequencies",
        lambda self, mol, conf_label, charge, mult: {
            "conf_label": conf_label,
            "frequencies": [-215.4, 11.0, 34.0],
            "n_negative": 1,
            "imaginary_analysis": [
                {
                    "frequency": -215.4,
                    "top_atoms": [
                        {"atom_idx": 0, "symbol": "C", "displacement": 1.4},
                        {"atom_idx": 4, "symbol": "Br", "displacement": 1.1},
                        {"atom_idx": 5, "symbol": "Cl", "displacement": 0.9},
                    ],
                }
            ],
        },
    )
    monkeypatch.setattr(
        "aqme.cmin.conformer_filters",
        lambda self, sorted_cids, cenergy, outmols: sorted_cids,
    )

    cmin(program="xtb", files=str(sdf_path), target="ts", freq=True)

    freq_file = _repo_path("tests", "cmin_TS", "CMIN", "frecuencies.dat")
    assert freq_file.exists()
    report = freq_file.read_text(encoding="utf-8")
    assert "Target: ts" in report
    assert "TS frequency summary" in report
    assert "Imaginary frequency: -215.40 cm^-1" in report
    assert any(symbol in report for symbol in ["(C)", "(Br)", "(Cl)"])
    freq_file.unlink(missing_ok=True)

# tests of basic QME optimizations (xtb, mace, aimnet2)
@pytest.mark.parametrize(
    "path, program, sdf, output_nummols",
    [
        ("complete", "xtb", "pentane_rdkit_methods.sdf", 4),
        ("complete", "mace", "pentane_rdkit_methods.sdf", 4), 
        ("complete", "aimnet2", "pentane_rdkit_methods.sdf", 4),
        ("partial", "xtb", "tests/cmin_methods/pentane_rdkit_methods.sdf", 4), 
        ("name", "xtb", "pentane_rdkit_methods.sdf", 4), 
    ],
)
def test_cmin_methods(
    path, program, sdf, output_nummols
):
    # runs the program with the different tests
    os.chdir(cmin_methods_dir)
    if path == 'complete':
        file_path = os.path.normpath(f'{cmin_methods_dir}/{sdf}')
        cmin(program=program, files=file_path)
        os.chdir(w_dir_main)
    elif path == 'partial':
        os.chdir(w_dir_main)
        file_path = os.path.normpath(sdf)
        cmin(program=program, files=file_path)
        os.chdir(cmin_methods_dir)
        sdf = 'pentane_rdkit_methods.sdf'
    elif path == 'name':
        os.chdir(cmin_methods_dir) 
        cmin(program=program, files=sdf)

    file = os.path.normpath(f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}.sdf')
    file2 = os.path.normpath(f'{cmin_methods_dir}/CMIN/All_confs/{sdf.split(".")[0]}_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)

    mols = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    mols_all = rdkit.Chem.SDMolSupplier(file2, removeHs=False, sanitize=False)
    
    assert len(mols_all) == output_nummols
    assert len(mols) == output_nummols

    # check that the optimizations work (different geometry than initial)
    outfile = open(file, "r")
    outlines = outfile.readlines()
    outfile.close()
    coords = ['2.5236','0.0073','0.1777']
    for coord in coords:
        assert coord not in outlines[4]

    # Force the release of the RDKit file handle for Windows
    del mols
    del mols_all
    os.chdir(w_dir_main)

# Test to check the application of constraints (distance, angle, and dihedral)
# CMIN is executed applying the three types of constraints supported by QME FixInternals
def test_cmin_constraints():
    os.chdir(cmin_methods_dir)
    sdf = "pentane_rdkit_methods.sdf"
    file_path = os.path.normpath(f'{cmin_methods_dir}/{sdf}')

    cmin(
        program="xtb",
        files=file_path,
        constraints_dist=[[0, 1, 1.54]],
        constraints_angle=[[0, 1, 2, 109.5]],
        constraints_dihedral=[[0, 1, 2, 3, 180.0]]
    )
    os.chdir(w_dir_main)

    file = os.path.normpath(f'{cmin_methods_dir}/CMIN/{sdf.split(".")[0]}.sdf')
    file2 = os.path.normpath(f'{cmin_methods_dir}/CMIN/All_confs/{sdf.split(".")[0]}_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)

@pytest.mark.parametrize(
    "test",
    [
        ("empty_values"),
    ],
)
def test_empty_values(
    test
):
    # runs the program with the different tests
    os.chdir(cmin_empty_dir)
    cmin(program='xtb', files=f'{cmin_empty_dir}/*.sdf')
    os.chdir(w_dir_main)

    file = os.path.normpath(f'{cmin_empty_dir}/CMIN/a.sdf')
    file2 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/a_all_confs.sdf')
    file3 = os.path.normpath(f'{cmin_empty_dir}/CMIN/b_fail.sdf')
    file4 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/b_fail_all_confs.sdf')
    file5 = os.path.normpath(f'{cmin_empty_dir}/CMIN/c.sdf')
    file6 = os.path.normpath(f'{cmin_empty_dir}/CMIN/All_confs/c_all_confs.sdf')

    assert os.path.exists(file)
    assert os.path.exists(file2)
    assert os.path.exists(file3)
    assert os.path.getsize(file3) == 0
    assert os.path.exists(file4)
    assert os.path.getsize(file4) == 0
    assert os.path.exists(file5)
    assert os.path.exists(file6)


@pytest.fixture(autouse=True, scope="module")
def cleanup_generated_files():
    yield

    for folder in [cmin_methods_dir, cmin_empty_dir, cmin_ts_dir]:
        cmin_out = os.path.join(folder, "CMIN")
        if os.path.exists(cmin_out):
            shutil.rmtree(cmin_out, ignore_errors=True)

        for dat_file in glob.glob(os.path.join(folder, "*.dat")):
            try:
                os.remove(dat_file)
            except OSError:
                pass

    os.chdir(w_dir_main)
