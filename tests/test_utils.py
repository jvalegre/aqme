from types import SimpleNamespace
import os
import sys

import numpy as np
import pytest
from rdkit import Chem

from aqme.cmin import normalize_cmin_backend
from aqme.qdescp import normalize_qdescp_runtime_options, qdescp
from aqme.utils import check_run, command_line_args


class FakePath:
    def __init__(self, name) -> None:
        self.name = name

    def __truediv__(self, p):
        return self

    def as_posix(self):
        return self.name

    def joinpath(self, p):
        return self.name + p


@pytest.mark.parametrize(
    "w_dir, path_exists, expected_folder_count, expected_resume_qcorr",
    [
        pytest.param(
            FakePath("/an/ordinary/working/directory"),
            True,
            1,
            False,
            id="no failed directory",
        ),
        pytest.param(
            FakePath("/a/working/directory/with/dir/named/failed_jobs/directory"),
            True,
            1,
            False,
            id='Shall not trigger exception: "folder_count referenced before assignment"',
        ),
        pytest.param(
            FakePath("/a/working/directory/with/failed/run_1/directory"),
            True,
            2,
            True,
            id="failed directory",
        ),
        pytest.param(
            FakePath("/a/working/directory/with/failed/run_1/directory"),
            False,
            2,
            True,
            id="failed run 1",
        ),
        pytest.param(
            FakePath("/a/working/directory/with/failed/run_13/directory"),
            True,
            14,
            True,
            id="failed run 13",
        ),
        pytest.param(
            FakePath(
                "/a/working/directory/with/my_testrun_12/directory/failed_results/run_1/"
            ),
            True,
            1,
            False,
            id="match the last run_xx",
        ),
    ],
)
def test_check_run(
    mocker, w_dir, path_exists, expected_folder_count, expected_resume_qcorr
):
    mocker.patch("pathlib.Path.as_posix", return_value=w_dir)
    mocker.patch("os.listdir", return_value=w_dir.name.split("/"))
    mocker.patch("os.path.exists", return_value=path_exists)

    try:
        assert (expected_folder_count, expected_resume_qcorr) == check_run(
            w_dir=w_dir
        )
    except UnboundLocalError as e:
        pytest.fail(f":: {e}")


def test_command_line_args_geom_opt_default(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["aqme", "--qdescp"])

    args = command_line_args()

    assert args.program is None
    assert args.geom_opt is True
    assert not hasattr(args, "xtb_opt")


def test_command_line_args_geom_opt_disable(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["aqme", "--qdescp", "--geom_opt", "False"])

    args = command_line_args()

    assert args.geom_opt is False

    input_sdf = "tests/cmin_methods/pentane_rdkit.sdf"
    dest_dir = "tests/qdescp_inputs/QDESCP_NO_OPT"

    # Read the initial coordinates before running QDESCP.
    with Chem.SDMolSupplier(input_sdf, removeHs=False) as supplier:
        mol_initial = next(mol for mol in supplier if mol is not None)
        coords_initial = mol_initial.GetConformer().GetPositions()

    # Run the QDESCP workflow with geometry optimization disabled.
    qdescp(files=[input_sdf], geom_opt=args.geom_opt, destination=dest_dir)

    output_sdf = os.path.join(dest_dir, "pentane_rdkit.sdf")
    with Chem.SDMolSupplier(output_sdf, removeHs=False) as supplier:
        mol_final = next(mol for mol in supplier if mol is not None)
        coords_final = mol_final.GetConformer().GetPositions()

    # Disabling geometry optimization must preserve the input coordinates.
    assert np.allclose(coords_initial, coords_final, atol=1e-5)


def test_command_line_args_geom_opt_enable_changes_coordinates(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["aqme", "--qdescp"])
    args = command_line_args()

    assert args.geom_opt is True

    input_sdf = "tests/cmin_methods/pentane_rdkit.sdf"
    dest_dir = "tests/qdescp_inputs/QDESCP_GEOM_OPT"

    # Read the initial coordinates before running QDESCP.
    with Chem.SDMolSupplier(input_sdf, removeHs=False) as supplier:
        mol_initial = next(mol for mol in supplier if mol is not None)
        coords_initial = mol_initial.GetConformer().GetPositions()

    # geom_opt defaults to True, so QDESCP must optimize the input geometry.
    qdescp(files=[input_sdf], geom_opt=args.geom_opt, destination=dest_dir)

    output_sdf = os.path.join(dest_dir, "pentane_rdkit.sdf")
    with Chem.SDMolSupplier(output_sdf, removeHs=False) as supplier:
        mol_final = next(mol for mol in supplier if mol is not None)
        coords_final = mol_final.GetConformer().GetPositions()

    assert not np.allclose(coords_initial, coords_final, atol=1e-5)


def test_normalize_qdescp_runtime_options_defaults_to_xtb():
    args = SimpleNamespace(program=None, geom_opt=True)

    normalized = normalize_qdescp_runtime_options(args)

    assert normalized.program == "xtb"


def test_normalize_cmin_backend_uses_tblite_for_xtb():
    args = SimpleNamespace(program="xtb")

    normalized = normalize_cmin_backend(args)

    assert normalized.program == "tblite"
