#!/usr/bin/env python

import pytest
from types import SimpleNamespace

from aqme.csearch.base import csearch
from aqme.csearch.utils import (
    resolve_racerts_charge,
    resolve_racerts_reacting_atoms,
)
from aqme.filter import conformer_filters
from rdkit.Chem import AllChem as Chem


class FakeLog:
    def __init__(self):
        self.messages = []
        self.finalized = False

    def write(self, message):
        self.messages.append(message)

    def finalize(self):
        self.finalized = True


def test_resolve_racerts_reacting_atoms_sdf_uses_atom_maps():
    mol = Chem.AddHs(Chem.MolFromSmiles("[CH3:1][CH2:2]O"))
    args = SimpleNamespace(freeze=[1, 2])
    log = FakeLog()

    reacting_atoms = resolve_racerts_reacting_atoms(args, mol, "example.sdf", log)

    assert reacting_atoms == [1, 2]


def test_resolve_racerts_reacting_atoms_xyz_uses_file_order():
    class FakeMol:
        def GetNumAtoms(self):
            return 5

    args = SimpleNamespace(freeze=[1, 3, 5])
    log = FakeLog()

    reacting_atoms = resolve_racerts_reacting_atoms(args, FakeMol(), "example.xyz", log)

    assert reacting_atoms == [1, 3, 5]


def test_resolve_racerts_charge_defaults_to_zero_for_xyz():
    class FakeMol:
        def HasProp(self, _):
            return False

    args = SimpleNamespace(charge=None)
    log = FakeLog()

    charge = resolve_racerts_charge(args, FakeMol(), "example.xyz", log)

    assert charge == 0
    assert any("recommended to set --charge explicitly" in msg for msg in log.messages)


def test_resolve_racerts_reacting_atoms_sdf_without_mapping_exits():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    args = SimpleNamespace(freeze=[1])
    log = FakeLog()

    with pytest.raises(SystemExit):
        resolve_racerts_reacting_atoms(args, mol, "example.sdf", log)


def test_resolve_racerts_reacting_atoms_sdf_invalid_map_reports_available_maps():
    mol = Chem.AddHs(Chem.MolFromSmiles("[CH3:1][CH2:2]O"))
    args = SimpleNamespace(freeze=[3])
    log = FakeLog()

    with pytest.raises(SystemExit):
        resolve_racerts_reacting_atoms(args, mol, "example.sdf", log)

    assert any("does not appear in the input molecule/SDF" in msg for msg in log.messages)
    assert any("Available map numbers are" in msg for msg in log.messages)


def test_validate_freeze_mode_blocks_crest():
    app = csearch.__new__(csearch)
    app.args = SimpleNamespace(
        freeze=[1],
        program="crest",
        log=FakeLog(),
    )

    with pytest.raises(SystemExit):
        app._validate_freeze_mode()


def test_log_racerts_citation_writes_additional_reference():
    app = csearch.__new__(csearch)
    log = FakeLog()
    app.args = SimpleNamespace(log=log, freeze=[1])

    app._log_racerts_citation()

    assert any("Rapid Generation of Transition-State Conformer Ensembles" in msg for msg in log.messages)


def test_select_conformers_force_filters_applies_filters_even_with_no_ff(monkeypatch):
    app = csearch.__new__(csearch)
    app.args = SimpleNamespace(log=FakeLog())
    outmols = [object(), object()]
    cenergy = [1.0, 0.0]
    called = {}

    def fake_conformer_filters(self, sorted_all_cids, cenergy_in, outmols_in, force_full_filters=False):
        called["value"] = (sorted_all_cids, cenergy_in, outmols_in, force_full_filters)
        return [1, 0]

    monkeypatch.setattr("aqme.csearch.base.conformer_filters", fake_conformer_filters)

    selected = app._select_conformers(
        outmols, cenergy, "example", "NO FF", force_filters=True
    )

    assert selected == [1, 0]
    assert called["value"][0] == [1, 0]
    assert called["value"][3] is True


def test_copy_racerts_atom_maps_preserves_original_mapping():
    source = Chem.AddHs(Chem.MolFromSmiles("[CH3:1][CH2:2]O"))
    target = Chem.Mol(source)
    for atom in target.GetAtoms():
        atom.SetAtomMapNum(0)

    app = csearch.__new__(csearch)
    copied = app._copy_racerts_atom_maps(source, target)

    assert [atom.GetAtomMapNum() for atom in copied.GetAtoms()] == [
        atom.GetAtomMapNum() for atom in source.GetAtoms()
    ]


def test_conformer_filters_freeze_bypasses_dot_smiles_rmsd_only(monkeypatch):
    app = csearch.__new__(csearch)
    app.args = SimpleNamespace(
        csearch=True,
        log=FakeLog(),
        ewin_cmin=5.0,
        initial_energy_threshold=0.0001,
        energy_threshold=0.25,
        rms_threshold=0.25,
        heavyonly=True,
        max_matches_rmsd=1000,
    )
    mol = Chem.MolFromSmiles("CC")
    mol.SetProp("SMILES", "A.B")
    outmols = [mol, mol]
    cenergy = [0.0, 1.0]

    def fail_if_called(*args, **kwargs):
        raise AssertionError("RMSD-only filter should not be used in freeze mode")

    monkeypatch.setattr("aqme.filter.apply_rmsd_only_filter", fail_if_called)
    monkeypatch.setattr("aqme.filter.apply_energy_window_filter", lambda *args, **kwargs: [0, 1])
    monkeypatch.setattr("aqme.filter.apply_pre_energy_filter", lambda *args, **kwargs: [0, 1])
    monkeypatch.setattr("aqme.filter.apply_rmsd_and_energy_filter", lambda *args, **kwargs: [0, 1])

    selected = conformer_filters(app, [0, 1], cenergy, outmols, force_full_filters=True)

    assert selected == [0, 1]
