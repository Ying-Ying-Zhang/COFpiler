import shutil

import numpy as np
import pytest
from ase.io import read
from conftest import CASES, DATA, EXAMPLE, REFERENCE

from cofpiler.cli import intercalate_main, stack_main


def test_stack_command_reproduces_reference(tmp_path):
    case = CASES["stacking"][1]
    shutil.copy(DATA / case["table"], tmp_path)
    shutil.copy(DATA / "hexagonal_monolayer.extxyz", tmp_path)
    stack_main(
        [
            "-data",
            case["table"],
            "-i",
            "hexagonal_monolayer.extxyz",
            "-path",
            str(tmp_path),
            "-T",
            str(case["temperature"]),
            "-s",
            case["symmetry"],
            "-L",
            str(case["layers"]),
            "-M",
            "2",
            "-o",
            "extxyz",
            "--seed",
            str(case["seed"]),
        ]
    )
    out = tmp_path / "statistical_structures"
    expected = read(REFERENCE / f"{case['name']}.extxyz")
    np.testing.assert_allclose(read(out / "8_0.extxyz").positions, expected.positions, atol=1e-6)
    assert (out / "8_1.extxyz").exists()
    record = (out / "record_8.txt").read_text()
    assert "model 0" in record and "model 1" in record


def test_intercalate_command_reproduces_reference(tmp_path):
    case = CASES["intercalate"][0]
    for f in EXAMPLE.iterdir():
        shutil.copy(f, tmp_path)
    intercalate_main(
        [
            "--data",
            "example_intercalate.xlsx",
            "--path",
            str(tmp_path),
            "-T",
            str(case["temperature"]),
            "-L",
            str(case["period"]),
            "--set",
            str(case["sets"]),
            "-o",
            "extxyz",
            "--seed",
            str(case["seed"]),
        ]
    )
    expected = read(REFERENCE / f"{case['name']}.extxyz")
    out = tmp_path / f"{expected.info['final_mode']}_1_{case['period']}l.extxyz"
    np.testing.assert_allclose(read(out).positions, expected.positions, atol=1e-6)
    assert (tmp_path / f"record_{case['period']}.txt").exists()


def test_input_errors_are_reported_without_traceback(tmp_path, capsys):
    shutil.copy(EXAMPLE / "example.xlsx", tmp_path)
    shutil.copy(DATA / "hexagonal_monolayer.extxyz", tmp_path)
    with pytest.raises(SystemExit) as exit_info:
        stack_main(
            [
                "-data",
                "example.xlsx",
                "-i",
                "hexagonal_monolayer.extxyz",
                "-path",
                str(tmp_path),
                "-L",
                "3",
            ]
        )
    assert exit_info.value.code == 2
    error = capsys.readouterr().err
    assert "cofpiler: error: the table has AB or ABC stacking types" in error
    assert "--ab-shift" in error


def test_missing_input_file_is_reported(tmp_path, capsys):
    with pytest.raises(SystemExit):
        intercalate_main(["-data", "missing.xlsx", "-path", str(tmp_path), "-L", "3", "--set", "1"])
    assert "missing.xlsx" in capsys.readouterr().err
