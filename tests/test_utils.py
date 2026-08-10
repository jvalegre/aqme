import sys
from types import SimpleNamespace

import pytest

from aqme.cmin import normalize_cmin_backend
from aqme.qdescp import normalize_qdescp_runtime_options
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
    assert args.method == "xtb"
    assert args.geom_opt is True
    assert not hasattr(args, "xtb_opt")


def test_command_line_args_geom_opt_disable(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["aqme", "--qdescp", "--geom_opt", "False"])

    args = command_line_args()

    assert args.geom_opt is False


def test_normalize_qdescp_runtime_options_rejects_method_without_geom_opt():
    args = SimpleNamespace(method="mace", program="xtb", geom_opt=False)

    with pytest.raises(ValueError, match="geom_opt=False cannot be combined"):
        normalize_qdescp_runtime_options(args, method_was_provided=True)


def test_normalize_cmin_backend_uses_tblite_for_xtb():
    args = SimpleNamespace(method=None, model=None, program=None)

    normalized = normalize_cmin_backend(args)

    assert normalized.method == "tblite"
    assert normalized.model == "tblite"
    assert normalized.program == "tblite"
