#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                   QDESCP module                    #
######################################################.

import os
import subprocess
import pytest
import pandas as pd
import numpy as np
import glob
import math
import shutil
from pathlib import Path
from aqme.qdescp import qdescp
from aqme.utils import load_sdf
from aqme.qdescp_utils import (
    read_json,
    get_descriptors,
    get_sdf_property,
    extract_smiles_from_file,
    validate_atom_mapping_consistency,
)
from aqme.csearch.utils import smiles_metadata_for_csearch
from rdkit.Chem import AllChem as Chem
from types import SimpleNamespace
from rdkit.Chem import AllChem as Chem
from types import SimpleNamespace

# saves the working directory
w_dir_main = os.getcwd()
qdescp_input_dir = Path(w_dir_main).joinpath("tests/qdescp_inputs")
qdescp_empty_dir = Path(w_dir_main).joinpath("tests/qdescp_empty")
qdescp_csv_dir = Path(w_dir_main).joinpath("tests/qdescp_csv")
qdescp_au_dir = Path(w_dir_main).joinpath("tests/qdescp_au")

GAS_CONSTANT = 8.3144621  # J / K / mol
J_TO_AU = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION
T = 298.15

# import descriptors
denovo_descriptors = get_descriptors('denovo')
interpret_descriptors = get_descriptors('interpret')
full_descriptors = get_descriptors('full')


class CaptureLog:
    def __init__(self):
        self.messages = []

    def write(self, message):
        self.messages.append(message)

    def finalize(self):
        pass


def test_qdescp_rejects_repeated_atom_map_number(tmp_path):
    sdf_file = tmp_path / "repeated_map.sdf"
    sdf_file.write_text(
        "\n"
        ">  <SMILES>\n"
        "[CH3:1][CH3:1]\n\n"
        "$$$$\n"
    )
    log = CaptureLog()

    assert not validate_atom_mapping_consistency(
        [sdf_file],
        {1},
        extract_smiles_from_file,
        log
    )
    assert "appears multiple times" in "".join(log.messages)


@pytest.mark.parametrize(
    "header",
    [
        ("SMILES,SMILES,code_name"),  # exact duplicate column name
        ("SMILES,smiles,code_name"),  # same name, different case
    ],
)
def test_qdescp_rejects_duplicate_smiles_columns(tmp_path, header):
    """Two SMILES columns (exact duplicate or differing only in case) must
    stop the run with a clear message instead of silently corrupting the
    SMILES column used downstream.
    """
    csv_path = tmp_path / "duplicate_smiles.csv"
    csv_path.write_text(f"{header}\nC,CC,mol_1\n", encoding="utf-8")

    processor = qdescp.__new__(qdescp)
    log = CaptureLog()
    processor.args = SimpleNamespace(csv_name=str(csv_path), log=log)

    with pytest.raises(SystemExit):
        processor._read_qdescp_csv()

    messages = "".join(log.messages)
    assert "Se han detectado dos columnas SMILES" in messages


def test_qdescp_mapped_atoms_keep_partial_charge_order():
    test_dir = qdescp_empty_dir / "mapped_charge_order"
    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir()
    input_csv = test_dir / "mapped_charge_order.csv"
    input_csv.write_text(
        "SMILES,code_name\n"
        "[H][C:2](=[O:3])[H:1],formaldehyde\n",
        encoding="utf-8",
    )
    output_files = [
        Path(w_dir_main) / f"AQME-ROBERT_{level}_mapped_charge_order.csv"
        for level in ("denovo", "interpret", "full")
    ]

    try:
        qdescp(
            input=str(input_csv),
            destination=str(test_dir / "QDESCP"),
            qdescp_atoms=[1, 2, 3],
            sample=1,
            nprocs=1,
        )

        for output_file in output_files:
            descriptors = pd.read_csv(output_file)
            charge_carbon = descriptors.loc[0, "Atom_2_C_Partial charge"]
            charge_hydrogen = descriptors.loc[0, "Atom_1_H_Partial charge"]
            charge_oxygen = descriptors.loc[0, "Atom_3_O_Partial charge"]

            assert charge_carbon > charge_hydrogen > charge_oxygen
    finally:
        for output_file in output_files:
            output_file.unlink(missing_ok=True)
        shutil.rmtree(test_dir, ignore_errors=True)


def test_qdescp_mapped_atom_charge_is_consistent_between_smiles():
    test_dir = qdescp_empty_dir / "mapped_charge_consistency"
    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir()
    smiles_inputs = [
        ("mapped_formaldehyde.csv", "[H][C:2](=[O:3])[H:1]", "mapped_formaldehyde"),
        ("mapped_carbonyl.csv", "[C:2]=O", "mapped_carbonyl"),
    ]
    output_files = [
        Path(w_dir_main) / f"AQME-ROBERT_{level}_{csv_name}"
        for csv_name, _, _ in smiles_inputs
        for level in ("denovo", "interpret", "full")
    ]

    try:
        partial_charges = []
        for csv_name, smiles, code_name in smiles_inputs:
            input_csv = test_dir / csv_name
            input_csv.write_text(
                f"SMILES,code_name\n{smiles},{code_name}\n",
                encoding="utf-8",
            )
            qdescp(
                input=str(input_csv),
                destination=str(test_dir / Path(csv_name).stem / "QDESCP"),
                qdescp_atoms=[2],
                sample=1,
                nprocs=1,
            )
            descriptors = pd.read_csv(
                Path(w_dir_main) / f"AQME-ROBERT_interpret_{csv_name}"
            )
            partial_charges.append(descriptors.loc[0, "Atom_2_C_Partial charge"])

        assert partial_charges[0] == pytest.approx(partial_charges[1], abs=5e-4)
    finally:
        for output_file in output_files:
            output_file.unlink(missing_ok=True)
        shutil.rmtree(test_dir, ignore_errors=True)


# tests for QDESCP-xTB
def _run_qdescp_xtb(input_csv, via_cli=False, stale_csearch=False, **qdescp_kwargs):
    """Reset the QDESCP/CSEARCH folders, run QDESCP-xTB on tests/qdescp_inputs/<input_csv>
    and return the folders and AQME-ROBERT databases generated."""
    folder_qdescp = qdescp_input_dir / "QDESCP"
    folder_csearch = qdescp_input_dir / "CSEARCH"
    for folder in (folder_qdescp, folder_csearch):
        if folder.exists():
            shutil.rmtree(folder)
    if stale_csearch:
        folder_csearch.mkdir(parents=True)
        folder_csearch.joinpath("stale_rdkit.sdf").write_text("stale")

    outputs = {
        level: Path(w_dir_main) / f"AQME-ROBERT_{level}_{input_csv}"
        for level in ("denovo", "interpret", "full")
    }
    for path in outputs.values():
        path.unlink(missing_ok=True)

    input_path = qdescp_input_dir / input_csv
    if via_cli:
        cmd_aqme = ["python", "-m", "aqme", "--qdescp",
                    "--input", str(input_path), "--destination", str(folder_qdescp)]
        for key, value in qdescp_kwargs.items():
            cmd_aqme += [f"--{key}", str(value)]
        subprocess.run(cmd_aqme)
    else:
        qdescp(input=str(input_path), destination=str(folder_qdescp), **qdescp_kwargs)

    for path in outputs.values():
        assert path.exists(), f"{path.name} was not generated"

    return SimpleNamespace(
        qdescp=folder_qdescp,
        csearch=folder_csearch,
        boltz=folder_qdescp / "boltz",
        input_cols=len(pd.read_csv(input_path).columns),
        **outputs,
    )


def _expected_descriptors(atom_prefix="P_"):
    """Molecular and atomic (with SMARTS prefix) descriptors expected at each level."""
    mol = {"denovo": denovo_descriptors["mol"]}
    mol["interpret"] = mol["denovo"] + interpret_descriptors["mol"]
    mol["full"] = mol["interpret"] + full_descriptors["mol"]

    atoms = {"denovo": denovo_descriptors["atoms"]}
    atoms["interpret"] = atoms["denovo"] + interpret_descriptors["atoms"]
    atoms["full"] = atoms["interpret"] + full_descriptors["atoms"]
    atoms = {level: [f"{atom_prefix}{d}" for d in descps] for level, descps in atoms.items()}
    return mol, atoms


def _check_levels(run, expected, desc_type, require_values=False, first_row_nan=()):
    """Each level must contain its descriptors and none of the higher-level ones.

    require_values: every value of the expected descriptors must be filled.
    first_row_nan: descriptors that must be NaN for the first molecule only.
    """
    for level in ("denovo", "interpret", "full"):
        df = pd.read_csv(getattr(run, level))
        for descp in expected[level]:
            assert descp in df.columns, f"{desc_type.capitalize()} descriptor {descp} is missing from {level}!"
            for i, val in enumerate(df[descp]):
                if require_values or (descp in first_row_nan and i > 0):
                    assert not pd.isna(val), f"{descp} is empty in row {i} of {level}"
                elif descp in first_row_nan:
                    assert pd.isna(val), f"{descp} should be empty in row {i} of {level}"
        for descp in expected["full"]:
            if descp not in expected[level]:
                assert descp not in df.columns, f"{desc_type.capitalize()} descriptor {descp} should not be in {level}!"


def _check_methane_without_targets(run, target_columns):
    """test_atom.csv/test_group.csv contain 4 molecules and mol_1 (methane) has none of
    the targeted atoms/groups, so only its atomic descriptors must be empty."""
    assert (run.qdescp / "mol_1_rdkit_conf_1.json").exists()
    assert (run.boltz / "mol_1_boltz.json").exists()

    df = pd.read_csv(run.interpret)
    assert len(df["HOMO"]) == 4
    assert df["HOMO"].isna().sum() == 0
    for col in target_columns:
        assert len(df[col]) == 4
        assert df[col].isna().sum() == 1, f"{col} should only be empty for methane"
    return df


def _check_atom_index(run):
    df = pd.read_csv(run.interpret)
    assert 'Atom_1_C_Partial charge' in df
    assert round(df['Atom_1_C_Partial charge'][1], 1) == -0.1


def _check_robert_outputs(run):
    """AQME-ROBERT databases start with code_name/SMILES and the raw QDESCP files are removed."""
    for level in ("denovo", "interpret", "full"):
        df = pd.read_csv(getattr(run, level))
        assert sorted(df.columns[:2].tolist(), key=str.lower) == ['code_name', 'SMILES']
        assert not (run.qdescp / "raw_data" / f"QDESCP_{level}_descriptors.csv").exists()


def test_qdescp_xtb_standard():
    run = _run_qdescp_xtb("test.csv")

    # 1) xTB parameters are stored correctly in the JSON files
    # 2 confs of mol_1 and 4 confs of mol_2
    assert len(glob.glob(f'{run.qdescp}/*.json')) == 6
    assert len(glob.glob(f'{run.qdescp}/*.xyz')) == 6

    energies_json, fermi_lvls_json = [], []
    for conf in ("mol_1_rdkit_conf_1", "mol_1_rdkit_conf_2"):
        assert (run.qdescp / f"{conf}.xyz").exists()
        json_data = read_json(run.qdescp / f"{conf}.json")
        energies_json.append(json_data["total energy"])
        fermi_lvls_json.append(json_data["Fermi-level"])

    energies_target = [-13.66512757, -13.66416741]
    fermi_lvls_target = [-4.5152, -4.3637]
    for i, _ in enumerate(energies_target):
        assert round(energies_target[i], 3) == round(energies_json[i], 3), \
            f"Energy mismatch at index {i}: Target = {round(energies_target[i],3)}, JSON = {round(energies_json[i],3)}"
        assert round(fermi_lvls_target[i], 1) == round(fermi_lvls_json[i], 1), \
            f"Fermi level mismatch at index {i}: Target = {round(fermi_lvls_target[i],1)}, JSON = {round(fermi_lvls_json[i],1)}"

    # 2) Boltzmann-averaged Fermi level matches the boltz JSON and the CSV
    boltz_terms = [
        math.exp(-(e - min(energies_json)) * J_TO_AU / GAS_CONSTANT / T)
        for e in energies_json
    ]
    fermi_lvl_boltz_calc = sum(
        fermi * term / sum(boltz_terms) for fermi, term in zip(fermi_lvls_json, boltz_terms)
    )
    fermi_lvl_boltz_file = read_json(run.boltz / "mol_1_boltz.json")["Fermi-level"]
    assert round(fermi_lvl_boltz_calc, 1) == round(fermi_lvl_boltz_file, 1), \
        f"Fermi level mismatch: calculated {fermi_lvl_boltz_calc} vs file {fermi_lvl_boltz_file}"

    # 3) CSV databases
    pd_boltz_interpret = pd.read_csv(run.interpret)
    assert round(fermi_lvl_boltz_calc, 1) == round(pd_boltz_interpret["Fermi-level"][0], 1)
    assert round(pd_boltz_interpret["HOMO"][0], 1) == -11.4
    assert round(pd_boltz_interpret["HOMO"][1], 1) == -11.3
    assert 'code_name' in pd_boltz_interpret.columns

    pd_boltz_full = pd.read_csv(run.full)
    assert pd_boltz_full["NumRotatableBonds"][0] == 3
    assert pd_boltz_full["NumRotatableBonds"][1] == 4

    mol, _ = _expected_descriptors()
    _check_levels(run, mol, 'mol', require_values=True)
    assert len(pd.read_csv(run.denovo).columns) == 11 == len(mol["denovo"]) + run.input_cols
    assert len(pd_boltz_interpret.columns) == 23 == len(mol["interpret"]) + run.input_cols
    assert len(pd_boltz_full.columns) == 240 # this might change in future RDKit versions

    # 4) raw databases (with atomic lists) are stored in the raw_data folder
    for level in ("denovo", "interpret", "full"):
        raw_csv = run.qdescp / "raw_data" / getattr(run, level).name
        assert raw_csv.exists()
        raw_df = pd.read_csv(raw_csv)
        assert 'code_name' in raw_df.columns
        raw_atom_vals = ["Partial charge", "Electrophil.", "Normaliz. nucleophil.", "Fukui+",
                         "Atom SASA", "Buried volume", "H bond H2O"]
        if level != "denovo":
            raw_atom_vals += ["Atom FOD", "Coord. numbers", "Atom Polarizability",
                              "Atom dispersion", "Pyramidalization", "Pyramidaliz. volume"]
        for raw_atom_val in raw_atom_vals:
            assert str(raw_df[raw_atom_val][0])[0] == '['


def test_qdescp_xtb_atom():
    run = _run_qdescp_xtb("test_atom.csv", qdescp_atoms=["P"])

    df = _check_methane_without_targets(run, ["P_Atom FOD", "P_Buried volume"])
    assert 'P_Electrophil.' in df
    assert round(df["HOMO"][0], 1) == -13.1


def test_qdescp_xtb_atom_index():
    run = _run_qdescp_xtb("test_idx.csv", qdescp_atoms=[1])
    _check_atom_index(run)


def test_qdescp_xtb_atom_index_command_line():
    run = _run_qdescp_xtb("test_idx.csv", via_cli=True, qdescp_atoms=[1])
    _check_atom_index(run)


def test_qdescp_xtb_group():
    run = _run_qdescp_xtb("test_group.csv", qdescp_atoms=["C=O"])

    df = _check_methane_without_targets(run, ["C=O_C_Atom FOD"])
    assert 'C=O_C_Partial charge' in df
    assert 'C=O_O_Partial charge' in df


def test_qdescp_xtb_multigroup():
    # Pd is included to check the try/except in the SMARTS pattern match,
    # and to check if the code works even if there are atoms that aren't used
    run = _run_qdescp_xtb("test_atom.csv", qdescp_atoms=["P", "CC", "Pd"])

    df = _check_methane_without_targets(run, ["P_Atom FOD", "P_Buried volume"])
    assert 'P_Electrophil.' in df


def test_qdescp_xtb_mapped_duplicates():
    # duplicate canonical SMILES with different mapped atoms
    run = _run_qdescp_xtb(
        "test_mapped_duplicates.csv", stale_csearch=True,
        qdescp_atoms=[1], sample=1, nprocs=1,
    )

    source_smiles = "[CH3:1]CC"
    alias_smiles = "C[CH2:1]C"
    source_metadata = smiles_metadata_for_csearch(source_smiles)
    alias_metadata = smiles_metadata_for_csearch(alias_smiles)
    source_file = run.csearch / "mol_a_rdkit.sdf"
    alias_file = run.csearch / "mol_b_rdkit.sdf"

    assert not (run.csearch / "stale_rdkit.sdf").exists()
    assert source_file.exists()
    assert alias_file.exists()
    assert get_sdf_property(source_file, "SMILES_INPUT") == source_smiles
    assert get_sdf_property(source_file, "AQME_ATOM_MAP") == source_metadata["atom_map"]
    assert get_sdf_property(alias_file, "SMILES_INPUT") == alias_smiles
    assert get_sdf_property(alias_file, "AQME_ATOM_MAP") == alias_metadata["atom_map"]
    assert get_sdf_property(alias_file, "AQME_ATOM_MAP") != source_metadata["atom_map"]

    pd_boltz_interpret = pd.read_csv(run.interpret)
    mapped_df = pd_boltz_interpret.set_index("code_name")
    assert len(pd_boltz_interpret["code_name"]) == 2
    assert not pd.isna(mapped_df.loc["mol_a", "Atom_1_C_Partial charge"])
    assert not pd.isna(mapped_df.loc["mol_b", "Atom_1_C_Partial charge"])


def test_qdescp_xtb_robert_atom():
    # AQME-ROBERT workflow with atomic descriptors
    run = _run_qdescp_xtb(
        "test_atom.csv", qdescp_atoms=["P", "CC", "Pd"],
        csv_name=str(qdescp_input_dir / "test_atom.csv"),
    )

    df = _check_methane_without_targets(run, ["P_Atom FOD", "P_Buried volume"])
    assert 'P_Electrophil.' in df
    assert 'Name' not in df
    assert 'SMILES' in df.columns

    mol, atoms = _expected_descriptors()
    _check_levels(run, mol, 'mol')
    _check_levels(run, atoms, 'atoms',
                  first_row_nan=('P_Partial charge', 'P_Buried volume', 'P_H bond H2O'))
    _check_robert_outputs(run)

    assert len(pd.read_csv(run.denovo).columns) == 20 == len(mol["denovo"]) + len(atoms["denovo"]) + run.input_cols
    assert len(df.columns) == 41 == len(mol["interpret"]) + len(atoms["interpret"]) + run.input_cols
    assert len(pd.read_csv(run.full).columns) == 258 # bunch of RDKit descps


def test_qdescp_xtb_robert_mol():
    # AQME-ROBERT workflow with NO atomic descriptors
    run = _run_qdescp_xtb("test_atom.csv", csv_name=str(qdescp_input_dir / "test_atom.csv"))

    assert 'SMILES' in pd.read_csv(run.interpret).columns

    mol, _ = _expected_descriptors()
    _check_levels(run, mol, 'mol', require_values=True)
    _check_robert_outputs(run)

    assert len(pd.read_csv(run.denovo).columns) == 11 == len(mol["denovo"]) + run.input_cols
    assert len(pd.read_csv(run.interpret).columns) == 23 == len(mol["interpret"]) + run.input_cols
    assert len(pd.read_csv(run.full).columns) == 240 # this might change in future RDKit versions


def test_qdescp_extra_column():
    """
    Extra columns in the input CSV besides code_name/SMILES (e.g. a target
    column used for AQME-ROBERT workflows) must be kept in the full, denovo
    and interpret output files, not just in full.
    """

    file = 'test_extra_column.csv'
    folder_qdescp = f'{qdescp_input_dir}/QDESCP_extra_column'
    if os.path.exists(folder_qdescp):
        shutil.rmtree(folder_qdescp)

    file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_{file}'
    file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_{file}'
    file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_{file}'
    for file_out in [file_descriptors_denovo, file_descriptors_interpret, file_descriptors_full]:
        if os.path.exists(file_out):
            os.remove(file_out)

    # QDESCP-xTB workflow with a methane input CSV that has an extra "my_target" column
    qdescp(
        input=f'{qdescp_input_dir}/{file}',
        destination=f'{folder_qdescp}',
    )

    input_cols = len(pd.read_csv(f'{qdescp_input_dir}/{file}').columns)

    for path in [file_descriptors_full, file_descriptors_denovo, file_descriptors_interpret]:
        df = pd.read_csv(path)
        assert 'code_name' in df.columns
        assert 'SMILES' in df.columns
        assert 'my_target' in df.columns, f"my_target column missing from {os.path.basename(path)}"
        assert df['my_target'][0] == 1.5

    # denovo/interpret must keep exactly descriptors + input CSV columns (code_name, SMILES, my_target)
    df_denovo = pd.read_csv(file_descriptors_denovo)
    df_interpret = pd.read_csv(file_descriptors_interpret)
    descp_denovo_mol = denovo_descriptors['mol']
    descp_interpret_mol = descp_denovo_mol + interpret_descriptors['mol']
    assert len(df_denovo.columns) == len(descp_denovo_mol) + input_cols
    assert len(df_interpret.columns) == len(descp_interpret_mol) + input_cols


@pytest.mark.parametrize(
    "test",
    [
        # tests for ignoring not working optimizations (i.e. qdescp keeps going even if an xTB calc fails)
        ("empty_values"),
    ],
)

def test_qdescp_missing(
    test
):

    # reset folder and files
    folder_qdescp = f'{qdescp_empty_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_AQME_run.csv'
    file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_AQME_run.csv'
    file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_AQME_run.csv'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    qdescp(
        files=f'{qdescp_empty_dir}/*.sdf',
        destination=f'{folder_qdescp}',
    )

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 3
    assert 'a' in df_interpret['code_name'][0]
    assert 'b_fail' in df_interpret['code_name'][1]
    assert 'c' in df_interpret['code_name'][2]
    assert round(df_interpret['HOMO'][0],1) == -11.8
    assert str(df_interpret['HOMO'][1]) == 'nan'

@pytest.mark.parametrize(
    "test",
    [
        # tests for using SDF as inputs and automated detection of common atom patterns
        ("sdf_input_n_auto"),
        ("sdf_input_[As]"),
    ],
)

def test_qdescp_sdf(
    test
):
    qdescp_sdf_dir = Path(w_dir_main).joinpath("tests/qdescp_sdf")

    if test == "sdf_input_[As]":
        qdescp_sdf_dir = Path(w_dir_main).joinpath("tests/qdescp_two_sdf")

    # reset folder and files
    folder_qdescp = f'{qdescp_sdf_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_AQME_run.csv'
    file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_AQME_run.csv'
    file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_AQME_run.csv'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    if test == "sdf_input_[As]":
        qdescp(
            files=f'{qdescp_sdf_dir}/*.sdf',
            destination=f'{folder_qdescp}',
            qdescp_atoms=['As'],
            charge=-1,
            mult=1,
        )
        name_1 = 'conf_72'
        name_2 = 'conf_73'
        atom = 'As'
        charge_1 = 0.35
        charge_2 = 0.31

    else:
        qdescp(
            files=f'{qdescp_sdf_dir}/*.sdf',
            destination=f'{folder_qdescp}',
        )
        name_1 = 'mol1'
        name_2 = 'mol_2'
        atom = 'P'
        charge_1 = 0.34
        charge_2 = 0.33

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 2
    assert name_1 == df_interpret['code_name'][0]
    assert name_2 == df_interpret['code_name'][1]
    assert len(df_interpret.columns) == 41

    # check if the automated detection of common pattern works
    assert charge_1 == round(df_interpret[f'{atom}_Partial charge'][0],2)
    assert charge_2 == round(df_interpret[f'{atom}_Partial charge'][1],2)


def test_qdescp_xyz_auto_charge_mult(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    xyz_structures = {
        'methane': '5\nmethane\nC 0.000 0.000 0.000\nH 0.629 0.629 0.629\nH -0.629 -0.629 0.629\nH -0.629 0.629 -0.629\nH 0.629 -0.629 -0.629\n',
        'ethane': '8\nethane\nC -0.770 0.000 0.000\nC 0.770 0.000 0.000\nH -1.157 0.513 0.889\nH -1.157 0.513 -0.889\nH -1.157 -1.026 0.000\nH 1.157 -0.513 -0.889\nH 1.157 -0.513 0.889\nH 1.157 1.026 0.000\n',
        'propane': '11\npropane\nC -1.270 0.000 0.000\nC 0.000 0.000 0.000\nC 1.270 0.000 0.000\nH -1.657 0.513 0.889\nH -1.657 0.513 -0.889\nH -1.657 -1.026 0.000\nH 0.000 0.000 1.089\nH 0.000 1.026 -0.363\nH 0.000 -1.026 -0.363\nH 1.657 0.513 0.889\nH 1.657 -0.513 0.889\n',
    }
    xyz_files = []
    for name, structure in xyz_structures.items():
        xyz_file = tmp_path / f'{name}.xyz'
        xyz_file.write_text(structure)
        xyz_files.append(str(xyz_file))

    qdescp(
        files=xyz_files,
        destination=str(tmp_path / 'QDESCP'),
    )

    for name in xyz_structures:
        sdf_file = tmp_path / 'CMIN' / f'{name}.sdf'
        assert get_sdf_property(sdf_file, 'Real charge') == '0'
        assert get_sdf_property(sdf_file, 'Mult') == '1'

@pytest.mark.parametrize(
    "file",
    [
        # tests for using SDF as inputs and automated detection of common atom patterns
        ("smiles_workflow.csv"),
    ],
)

def test_qdescp_csv(
    file
):

    # reset folder and files
    folder_qdescp = f'{qdescp_csv_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_{file}'
    file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_{file}'
    file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_{file}'
    if os.path.exists(file_descriptors_denovo):
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret):
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full):
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    qdescp(
        input=f'{qdescp_csv_dir}/{file}',
        destination=f'{folder_qdescp}',
    )

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 2
    assert 'mol_1' == df_interpret['code_name'][0]
    assert 'mol_2' == df_interpret['code_name'][1]
    assert len(df_interpret.columns) == 23

    csearch_dir = f'{os.path.dirname(folder_qdescp)}/CSEARCH'
    n_conformers_mol1 = len(load_sdf(f'{csearch_dir}/mol_1_rdkit.sdf'))
    n_conformers_mol2 = len(load_sdf(f'{csearch_dir}/mol_2_rdkit.sdf'))
    assert 2 <= n_conformers_mol1 <= 5
    assert 4 <= n_conformers_mol2 <= 5

    # check that the xTB version is printed
    f = open(f'{w_dir_main}/QDESCP_data.dat', "r")
    data = f.readlines()
    f.close()

    version_print = False
    for line in data:
        if 'xTB version used: 6.7.1' in line:
            version_print = True
    assert version_print

@pytest.mark.parametrize(
    "file, run_test",
    [
        # tests for using CSV inputs with charges/mult and automated detection of common atom patterns
        ("Au_test.csv",1),
        # overwrite charge/mult through the command line
        ("Au_test.csv",2),
    ],
)

def test_au_csv(
    file,run_test
):

    # reset folder and files
    folder_qdescp = f'{qdescp_au_dir}/QDESCP'
    folder_boltz = f'{folder_qdescp}/boltz'
    for folder in [folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    file_descriptors_interpret = f'{w_dir_main}/AQME-ROBERT_interpret_{file}'
    file_descriptors_full = f'{w_dir_main}/AQME-ROBERT_full_{file}'
    file_descriptors_denovo = f'{w_dir_main}/AQME-ROBERT_denovo_{file}'
    if os.path.exists(file_descriptors_denovo): 
        os.remove(file_descriptors_denovo)
    if os.path.exists(file_descriptors_interpret): 
        os.remove(file_descriptors_interpret)
    if os.path.exists(file_descriptors_full): 
        os.remove(file_descriptors_full)

    # QDESCP-xTB workflow
    qdescp_kwargs = {
        "input": f'{qdescp_au_dir}/{file}',
        "destination": f'{folder_qdescp}',
    }
    
    if run_test == 2:
        qdescp_kwargs.update({
            "charge": 0,
            "mult": 1
        })

    qdescp(**qdescp_kwargs)

    expected_charge_mult = {
        '200': ('0', '1'),
        '201': ('2', '3'),
    } if run_test == 1 else {
        '200': ('0', '1'),
        '201': ('0', '1'),
    }
    for code_name, (charge, mult) in expected_charge_mult.items():
        cmin_file = f'{qdescp_au_dir}/CMIN/{code_name}_rdkit.sdf'
        assert get_sdf_property(cmin_file, 'Real charge') == charge
        assert get_sdf_property(cmin_file, 'Mult') == mult

    # Checking molecular descriptors
    descp_denovo_mol = denovo_descriptors['mol'] 
    descp_denovo_atoms = denovo_descriptors['atoms']

    descp_interpret_mol = descp_denovo_mol + interpret_descriptors['mol']
    descp_interpret_atoms = descp_denovo_atoms + interpret_descriptors['atoms']

    descp_full_mol = descp_interpret_mol + full_descriptors['mol']
    descp_full_atoms = descp_interpret_atoms + full_descriptors['atoms']

    # add Au_ prefix for SMARTS pattern
    descp_denovo_atoms = [f'Au_{descp}' for descp in descp_denovo_atoms]
    descp_interpret_atoms = [f'Au_{descp}' for descp in descp_interpret_atoms]
    descp_full_atoms = [f'Au_{descp}' for descp in descp_full_atoms]

    # checking csv file
    df_interpret = pd.read_csv(file_descriptors_interpret)
    assert len(df_interpret['code_name']) == 2
    assert '200' == str(df_interpret['code_name'][0])
    assert '201' == str(df_interpret['code_name'][1])
    assert round(df_interpret['HOMO'][0],1) == -9.6
    if run_test == 1:
        assert round(df_interpret['HOMO'][1],1) == -22.6
    elif run_test == 2:
        assert round(df_interpret['HOMO'][1],1) == -8.1

    input_cols = len(pd.read_csv(f'{qdescp_au_dir}/{file}').columns)
    assert len(df_interpret.columns) == 41 == len(descp_interpret_mol)+len(descp_interpret_atoms)+input_cols # input CSV columns: SMILES and code_name

    # Checking molecular and atomic descriptors
    def check_descriptors_Au(pd_boltz, descriptors, excluded_descriptors, desc_type):
        """
        Function to check the presence and absence of descriptors in the DataFrame.
        pd_boltz: Pandas DataFrame with the calculated descriptors.
        descriptors: List of descriptors that should be present.
        excluded_descriptors: List of descriptors that should not be present.
        desc_type: Type of descriptors ('mol' or 'atoms') for printing in messages.
        """
        # Check for the presence of descriptors
        for descp in descriptors:
            for _,val in enumerate(pd_boltz[descp]):
                if descp in ['Au_Partial charge','Au_Buried volume','Au_H bond H2O']:
                    assert str(val).lower() != 'nan'
                assert descp in pd_boltz.columns, f"{desc_type.capitalize()} descriptor {descp} is missing from columns!"

        # Check for the absence of descriptors that should not be present
        for descp in excluded_descriptors:
            assert descp not in pd_boltz.columns, f"{desc_type.capitalize()} descriptor {descp} should not be present in columns!"

    check_descriptors_Au(df_interpret, descp_interpret_mol, [d for d in descp_full_mol if d not in descp_interpret_mol], 'mol')
    check_descriptors_Au(df_interpret, descp_interpret_atoms, [d for d in descp_full_atoms if d not in descp_interpret_atoms], 'atoms')
    
    # only two files should remain inside the xtb_data folders
    xtb_data_path = f'{folder_qdescp}'
    format_check = ['.json','.xyz']
    for file_Au in ["200_rdkit_conf_1","201_rdkit_conf_1"]:
        for fmt in format_check:
            assert os.path.exists(f'{xtb_data_path}/{file_Au}{fmt}')
    assert len(glob.glob(f'{xtb_data_path}/*.json')) == 2
    assert len(glob.glob(f'{xtb_data_path}/*.xyz')) == 2

# tests for QDESCP-NMR
@pytest.mark.parametrize(
    "json_files",
    [
        ("*.json")
    ]
)

def test_qdescp_nmr(json_files):

    # reset folder
    folder_csearch = f'{qdescp_input_dir}/CSEARCH'
    folder_qdescp = f'{qdescp_input_dir}/QDESCP'
    folder_boltz = f'{qdescp_input_dir}/boltz'
    for folder in [folder_csearch,folder_qdescp,folder_boltz]:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    json_files = f'{qdescp_input_dir}/test_conf_*.json'
    nmr_atoms = [6, 1]  # [C,H]
    nmr_slope=[-1.0537, -1.0784]
    nmr_intercept=[181.7815,31.8723]

    # QDESCP-NMR workflow
    qdescp(
        files=json_files,
        destination=qdescp_input_dir,
        program="nmr",
        nmr_slope=nmr_slope,
        nmr_intercept=nmr_intercept,
        nmr_experim=f'{qdescp_input_dir}/Experimental_NMR_shifts.csv',
    )

    # read json files and calculate Boltzmann shifts
    energies_json,nmr_json = [],[]
    for json_file in glob.glob(json_files):
        json_data = read_json(json_file)
        energies_json.append(json_data["optimization"]["scf"]["scf energies"][-1])
        # retrieves and scales NMR shifts from json files
        atoms = json_data["atoms"]["elements"]["number"]
        tensor = json_data["properties"]["NMR"]["NMR isotopic tensors"]
        shifts = {}
        i = 0
        for atom, ten in zip(atoms, tensor):
            if atom in nmr_atoms:
                # assigning values from arrays
                index = nmr_atoms.index(atom)
                slope_nuc = nmr_slope[index]
                intercept_nuc = nmr_intercept[index]
                scaled_nmr = (intercept_nuc - ten) / (-slope_nuc)
                shifts[i] = scaled_nmr
            else:
                pass
            i += 1
        nmr_json.append(shifts)

    # calculate Boltzmann averaged values
    energ = [number - min(energies_json) for number in energies_json]
    boltz_sum = 0.0
    for e in energ:
        boltz_sum += math.exp(-e * J_TO_AU / GAS_CONSTANT / T)
    weights = []
    for e in energ:
        weight = math.exp(-e * J_TO_AU / GAS_CONSTANT / T) / boltz_sum
        weights.append(weight)

    boltz_avg = []
    for i, p in enumerate(nmr_json):
        boltz_avg.append([number * weights[i] for number in p.values()])
    nmr_json_calc = np.sum(boltz_avg, 0)

    # check that the averaged NMR shifts are the same as the values included in the json file
    folder_boltz = f'{qdescp_input_dir}/boltz'
    json_data_boltz = read_json(f'{folder_boltz}/test_boltz.json')
    nmr_json_file = list(json_data_boltz["NMR Chemical Shifts"].values())
    for i,_ in enumerate(nmr_json_calc):
        assert round(nmr_json_calc[i],2) == round(nmr_json_file[i],2)

    # check that the averaged NMR shifts are the same as the values included in the csv file
    pd_boltz = pd.read_csv(f'{qdescp_input_dir}/Experimental_NMR_shifts_predicted.csv')
    nmr_boltz_csv = pd_boltz["boltz_avg"]
    nmr_json_calc_1H = nmr_json_calc[21:]
    for i,_ in enumerate(nmr_json_calc_1H):
        assert round(nmr_json_calc_1H[i],2) == round(nmr_boltz_csv[i],2)

    # check that the shift errors in the csv file are correct
    nmr_experim_csv = pd_boltz["experimental_ppm"]
    error_calc = abs(nmr_boltz_csv - nmr_experim_csv)
    error_csv = pd_boltz["error_boltz"]
    for i,_ in enumerate(error_calc):
        if str(error_calc[i]) not in ['nan']:
            assert round(error_calc[i],2) == round(error_csv[i],2)


@pytest.mark.parametrize('program', ['rdkit', 'crest'])
def test_qdescp_conformer_xyz_files_follow_sdf_order(tmp_path, program):
    """Each conformer XYZ must be paired with the charge/mult of that conformer.

    OpenBabel writes mol_conf_1.xyz ... mol_conf_12.xyz and glob() returns them in
    an arbitrary (usually lexicographic) order, so _conf_10 used to be paired with
    the charge and multiplicity of _conf_2. Checked for the SDF files produced by
    both CSEARCH programs (rdkit and crest).
    """
    n_confs = 12
    name = f'mol_{program}'
    sdf_file = tmp_path / f'{name}.sdf'

    # a different charge/mult per conformer, so any mismatch is detectable
    charges = list(range(n_confs))
    mults = [i + 1 for i in range(n_confs)]

    with Chem.SDWriter(str(sdf_file)) as writer:
        for i in range(n_confs):
            mol = Chem.AddHs(Chem.MolFromSmiles('C'))
            Chem.EmbedMolecule(mol, randomSeed=i + 1)
            mol.SetProp('_Name', f'{name} {i + 1}')
            mol.SetProp('Real charge', str(charges[i]))
            mol.SetProp('Mult', str(mults[i]))
            writer.write(mol)

    # stand in for the XYZ files that OpenBabel writes with the -m option
    for i in range(n_confs):
        (tmp_path / f'{name}_conf_{i + 1}.xyz').write_text('1\n\nC 0.0 0.0 0.0\n')

    processor = qdescp.__new__(qdescp)
    processor.args = SimpleNamespace(charge=None, mult=None)

    xyz_files, xyz_charges, xyz_mults = processor._process_other_conformers(
        str(sdf_file), name
    )

    assert len(xyz_files) == n_confs
    # conformer order, not lexicographic order (_conf_10 after _conf_9)
    assert [Path(f).name for f in xyz_files] == [
        f'{name}_conf_{i + 1}.xyz' for i in range(n_confs)
    ]
    assert [int(charge) for charge in xyz_charges] == charges
    assert [int(mult) for mult in xyz_mults] == mults


def test_qdescp_conformer_xyz_paths_are_not_duplicated(tmp_path):
    """XYZ inputs must return usable paths.

    The directory of the input file used to be prepended to paths that glob()
    already returned as absolute, giving unusable paths such as
    C:/dir/C:/dir/mol_conf_1.xyz.
    """
    n_confs = 3
    name = 'mol_xyz'
    xyz_file = tmp_path / f'{name}.xyz'
    xyz_file.write_text('1\n\nC 0.0 0.0 0.0\n')

    for i in range(n_confs):
        (tmp_path / f'{name}_conf_{i + 1}.xyz').write_text(
            f'1\ncharge={i} mult=1\nC 0.0 0.0 0.0\n'
        )

    processor = qdescp.__new__(qdescp)
    processor.args = SimpleNamespace(charge=None, mult=None)

    xyz_files, xyz_charges, xyz_mults = processor._process_xyz_conformers(
        str(xyz_file), name
    )

    assert len(xyz_files) == n_confs
    for conf_file in xyz_files:
        assert os.path.isfile(conf_file), f'{conf_file} is not a usable path'
    assert xyz_charges == list(range(n_confs))
    assert xyz_mults == [1] * n_confs


def _write_qdescp_sdf(path, smiles):
    """Write a single-molecule SDF, optionally carrying a <SMILES> property."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles if smiles is not None else 'C'))
    Chem.EmbedMolecule(mol, randomSeed=1)
    mol.SetProp('_Name', Path(path).stem)
    if smiles is not None:
        mol.SetProp('SMILES', smiles)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)


def test_qdescp_keeps_files_without_smiles(tmp_path):
    """Files with no SMILES property (xyz/json inputs, or sdf files not coming
    from a QDESCP CSV run) must be kept as-is: SMILES-based duplicate
    detection simply does not apply to them, and the run must not stop.
    """
    with_smiles = tmp_path / 'mol_1_rdkit.sdf'
    without_smiles = tmp_path / 'mol_2_rdkit.sdf'
    _write_qdescp_sdf(with_smiles, 'CC')
    _write_qdescp_sdf(without_smiles, None)

    processor = qdescp.__new__(qdescp)
    log = CaptureLog()
    processor.args = SimpleNamespace(
        files=[str(with_smiles), str(without_smiles)], log=log
    )

    unique_files = processor.get_unique_files()

    assert unique_files == [str(with_smiles), str(without_smiles)]
    assert ''.join(log.messages) == ''


def test_qdescp_keeps_unique_smiles_and_warns_about_duplicates(tmp_path):
    """Inputs that all carry a SMILES are deduplicated and the duplicate reported."""
    first = tmp_path / 'mol_1_rdkit.sdf'
    duplicate = tmp_path / 'mol_2_rdkit.sdf'
    other = tmp_path / 'mol_3_rdkit.sdf'
    _write_qdescp_sdf(first, 'CC')
    _write_qdescp_sdf(duplicate, 'CC')
    _write_qdescp_sdf(other, 'CCO')

    processor = qdescp.__new__(qdescp)
    log = CaptureLog()
    processor.args = SimpleNamespace(
        files=[str(first), str(duplicate), str(other)], log=log
    )

    unique_files = processor.get_unique_files()

    assert unique_files == [str(first), str(other)]
    messages = ''.join(log.messages)
    assert 'mol_2_rdkit.sdf' in messages
    assert 'same SMILES' in messages


def test_qdescp_csv_stops_when_smiles_column_has_missing_cells(tmp_path):
    """The run must stop right after reading the CSV (before conformer
    generation starts) if the SMILES column has fewer filled cells than
    another column, i.e. some rows are missing their SMILES while other
    columns of that same row are filled in.
    """
    csv_path = tmp_path / 'incomplete_smiles.csv'
    pd.DataFrame({
        'code_name': ['mol_1', 'mol_2', 'mol_3'],
        'SMILES': ['CC', '', 'CCO'],
    }).to_csv(csv_path, index=False)

    processor = qdescp.__new__(qdescp)
    log = CaptureLog()
    processor.args = SimpleNamespace(csv_name=str(csv_path), log=log)

    with pytest.raises(SystemExit):
        processor._read_qdescp_csv()

    messages = ''.join(log.messages)
    assert 'mol_2' in messages
    assert 'SMILES column are filled' in messages


def test_qdescp_csv_allows_fully_blank_rows(tmp_path):
    """A row that is entirely blank (e.g. a trailing empty row from Excel)
    must not be treated as a missing-SMILES error.
    """
    csv_path = tmp_path / 'blank_row.csv'
    pd.DataFrame({
        'code_name': ['mol_1', 'mol_2', ''],
        'SMILES': ['CC', 'CCO', ''],
    }).to_csv(csv_path, index=False)

    processor = qdescp.__new__(qdescp)
    log = CaptureLog()
    processor.args = SimpleNamespace(csv_name=str(csv_path), log=log)

    df_qdescp = processor._read_qdescp_csv()

    assert len(df_qdescp) == 3
    assert ''.join(log.messages) == ''


@pytest.mark.parametrize(
    "solvent_options,reported_option",
    [
        ({'qdescp_solvent': 'h2o'}, '--qdescp_solvent'),
        ({'xtb_keywords': '--alpb h2o'}, '--alpb'),
        ({'xtb_keywords': '--gbsa h2o'}, '--gbsa'),
        ({'qdescp_solvent': 'h2o', 'xtb_keywords': '--alpb h2o'}, '--qdescp_solvent'),
    ],
)
def test_qdescp_stops_when_a_solvent_is_requested(
    tmp_path, monkeypatch, capsys, solvent_options, reported_option
):
    """QDESCP descriptors are calculated with PTB, which has no implicit solvation
    model, so asking for a solvent with --qdescp_solvent or with xTB solvation
    keywords must stop the run instead of being silently ignored.
    """
    monkeypatch.chdir(tmp_path)
    xyz_file = tmp_path / 'methane.xyz'
    xyz_file.write_text(
        '5\nmethane\nC 0.000 0.000 0.000\nH 0.629 0.629 0.629\n'
        'H -0.629 -0.629 0.629\nH -0.629 0.629 -0.629\nH 0.629 -0.629 -0.629\n'
    )

    with pytest.raises(SystemExit):
        qdescp(
            files=[str(xyz_file)],
            destination=str(tmp_path / 'QDESCP'),
            **solvent_options,
        )

    output = capsys.readouterr().out
    assert 'does not support solvation' in output
    assert reported_option in output

