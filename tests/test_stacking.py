import numpy as np
import pytest
from ase.io import read
from conftest import EXAMPLE

from cofpiler import StackingMode, StackingTable, build_stacked_model, ureg


@pytest.fixture
def table(aa_ab_table_path):
    return StackingTable.read(aa_ab_table_path)


@pytest.fixture
def monolayer(monolayer_path):
    return read(monolayer_path)


def test_builds_one_monolayer_per_layer(monolayer, table):
    structure, records = build_stacked_model(monolayer, table, n_layers=5, temperature=293)
    assert len(structure) == 5 * len(monolayer)
    assert [r.layer for r in records] == [2, 3, 4, 5]
    assert set(r.mode for r in records) <= set(table.names)


def test_layers_move_up_by_the_interlayer_distance(monolayer, table):
    structure, records = build_stacked_model(monolayer, table, n_layers=4, temperature=293)
    z = structure.positions[:, 2].reshape(4, len(monolayer))[:, 0]
    np.testing.assert_allclose(np.diff(z), 3.4)
    np.testing.assert_allclose(records[-1].total_shift[2], 3 * 3.4)


def test_same_seed_gives_same_model(monolayer, table):
    first, _ = build_stacked_model(monolayer, table, 6, 5000, rng=np.random.RandomState(7))
    second, _ = build_stacked_model(monolayer, table, 6, 5000, rng=np.random.RandomState(7))
    np.testing.assert_array_equal(first.positions, second.positions)


def test_abc_types_are_not_supported(monolayer):
    with pytest.raises(NotImplementedError, match="ABC_ecl"):
        build_stacked_model(monolayer, StackingTable.read(EXAMPLE / "example.xlsx"), 4, 293)


def test_unknown_slip_kind_is_rejected(monolayer):
    table = StackingTable(
        modes=(StackingMode(name="AA_x", energy="0 kJ/mol", shift=ureg.Quantity([0, 0, 3.4], "Å")),)
    )
    with pytest.raises(ValueError, match="ecl, s1, s2"):
        build_stacked_model(monolayer, table, 3, 293)


@pytest.mark.parametrize(
    "kwargs, error",
    [
        ({"n_layers": 1}, "at least 2"),
        ({"symmetry": 5}, "symmetry must be one of"),
        ({"symmetry": 4}, "needs ab_slip"),
        ({"symmetry": 4, "ab_slip": "axis"}, "expects a stacking type 'AB_s2_diag'"),
    ],
)
def test_invalid_arguments(monolayer, table, kwargs, error):
    arguments = {"n_layers": 3, "temperature": 293, **kwargs}
    with pytest.raises(ValueError, match=error):
        build_stacked_model(monolayer, table, **arguments)


def test_c4_sets_the_other_ab_slip_energy_to_zero(monolayer):
    def m(name, energy):
        shift = ureg.Quantity([0.9, 0.7, 3.4], "angstrom")
        return StackingMode(name=name, energy=f"{energy} kJ/mol", shift=shift)

    table = StackingTable(modes=(m("AA_ecl", 40.0), m("AB_s2_axis", 40.0), m("AB_s2_diag", 40.0)))
    _, records = build_stacked_model(
        monolayer, table, 40, 293, symmetry=4, ab_slip="axis", rng=np.random.RandomState(0)
    )
    assert {r.mode for r in records} == {"AB_s2_diag"}
    assert table.energies()[2] == 40.0  # the input table is unchanged
