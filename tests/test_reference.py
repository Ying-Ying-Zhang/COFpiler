"""The package reproduces the 2022 scripts (see tests/data/make_reference.py).

For stacking, only models of AA layers: AB and ABC layers differ from 2022 on purpose.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.io import read
from conftest import CASES, DATA, EXAMPLE, REFERENCE

from cofpiler import StackingTable, build_intercalated_model, build_stacked_model


def assert_same_structure(actual: Atoms, expected: Atoms) -> None:
    assert actual.get_chemical_symbols() == expected.get_chemical_symbols()
    np.testing.assert_allclose(actual.positions, expected.positions, atol=1e-6)
    np.testing.assert_allclose(actual.cell.array, expected.cell.array, atol=1e-6)


@pytest.mark.parametrize("case", CASES["stacking"], ids=lambda c: c["name"])
def test_aa_stacking_matches_2022_script(case):
    structure, _ = build_stacked_model(
        read(DATA / "hexagonal_monolayer.extxyz"),
        StackingTable.read(DATA / case["table"]),
        n_layers=case["layers"],
        temperature=case["temperature"],
        symmetry=int(case["symmetry"][1:]),
        mirror=case["mirror"],
        ab_shift="slip",
        rng=np.random.RandomState(case["seed"]),
    )
    assert_same_structure(structure, read(REFERENCE / f"{case['name']}.extxyz"))


@pytest.mark.parametrize("case", CASES["intercalate"], ids=lambda c: c["name"])
def test_intercalation_matches_2022_script(case):
    structure, mode, _ = build_intercalated_model(
        EXAMPLE,
        StackingTable.read(EXAMPLE / "example_intercalate.xlsx"),
        period=case["period"],
        n_sets=case["sets"],
        temperature=case["temperature"],
        rng=np.random.RandomState(case["seed"]),
    )
    expected = read(REFERENCE / f"{case['name']}.extxyz")
    assert mode == expected.info["final_mode"]
    assert_same_structure(structure, expected)
