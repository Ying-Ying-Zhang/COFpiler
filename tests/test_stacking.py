import numpy as np
import pytest
from ase.io import read

from cofpiler import StackingMode, StackingTable, build_stacked_model, ureg
from cofpiler.geometry import rotate_xy
from cofpiler.stacking import ab_lattice_offset


@pytest.fixture
def table(aa_ab_table_path):
    return StackingTable.read(aa_ab_table_path)


@pytest.fixture
def monolayer(monolayer_path):
    return read(monolayer_path)


@pytest.fixture
def offset(monolayer):
    return ab_lattice_offset(monolayer.cell, 6)[:2]


def mode(name, energy, xy, z=3.4):
    shift = ureg.Quantity([*xy, z], "angstrom")
    return StackingMode(name=name, energy=f"{energy} kJ/mol", shift=shift)


def images(vector):
    """The vector turned by every multiple of 60 degrees."""
    return [rotate_xy(vector, i * np.pi / 3) for i in range(6)]


def test_builds_one_monolayer_per_layer(monolayer, table):
    structure, records = build_stacked_model(monolayer, table, 5, 293, ab_shift="slip")
    assert len(structure) == 5 * len(monolayer)
    assert [r.layer for r in records] == [2, 3, 4, 5]
    assert set(r.mode for r in records) <= set(table.names)


def test_layers_move_up_by_the_interlayer_distance(monolayer, table):
    structure, records = build_stacked_model(monolayer, table, 4, 293, ab_shift="slip")
    z = structure.positions[:, 2].reshape(4, len(monolayer))[:, 0]
    np.testing.assert_allclose(np.diff(z), 3.4)
    np.testing.assert_allclose(records[-1].total_shift[2], 3 * 3.4)


def test_same_seed_gives_same_model(monolayer, table):
    first, _ = build_stacked_model(
        monolayer, table, 6, 5000, ab_shift="slip", rng=np.random.RandomState(7)
    )
    second, _ = build_stacked_model(
        monolayer, table, 6, 5000, ab_shift="slip", rng=np.random.RandomState(7)
    )
    np.testing.assert_array_equal(first.positions, second.positions)


def test_abc_layers_shift_to_an_ab_position(monolayer, offset):
    table = StackingTable(modes=(mode("ABC_ecl", 0.0, (0.0, 0.0), z=2.8),))
    _, records = build_stacked_model(
        monolayer, table, 8, 293, ab_shift="slip", rng=np.random.RandomState(0)
    )
    for r in records:
        assert any(np.allclose(r.shift[:2], o) for o in images(offset))
        assert r.shift[2] == 2.8


def test_ab_layers_keep_the_slip_from_the_table(monolayer, offset):
    # In the 2022 script the offset replaced the slip; now each shift is slip + offset.
    slip = np.array([0.0, 1.7])
    table = StackingTable(modes=(mode("AB_s1", 0.0, slip),))
    _, records = build_stacked_model(
        monolayer, table, 8, 293, ab_shift="slip", rng=np.random.RandomState(0)
    )
    for r in records:
        assert any(np.allclose(r.shift[:2], s + o) for s in images(slip) for o in images(offset))


def test_slip_and_full_tables_give_the_same_model(monolayer, offset):
    slip = np.array([0.9, 0.7])
    aa = mode("AA_s1_inc", 1.0, (0.0, 1.7))
    as_slip = StackingTable(modes=(aa, mode("AB_s2", 0.0, slip), mode("ABC_s1", 0.5, slip)))
    as_full = StackingTable(
        modes=(aa, mode("AB_s2", 0.0, slip + offset), mode("ABC_s1", 0.5, slip + offset))
    )
    first, _ = build_stacked_model(
        monolayer, as_slip, 8, 2000, ab_shift="slip", rng=np.random.RandomState(3)
    )
    second, _ = build_stacked_model(
        monolayer, as_full, 8, 2000, ab_shift="full", rng=np.random.RandomState(3)
    )
    np.testing.assert_allclose(first.positions, second.positions, atol=1e-12)


def test_ab_and_abc_rows_need_ab_shift(monolayer, table):
    with pytest.raises(ValueError, match="ab_shift"):
        build_stacked_model(monolayer, table, 3, 293)


def test_abc_is_not_implemented_for_c4(monolayer):
    table = StackingTable(modes=(mode("ABC_ecl", 0.0, (0.0, 0.0)),))
    with pytest.raises(NotImplementedError, match="C4"):
        build_stacked_model(monolayer, table, 3, 293, symmetry=4, ab_slip="axis", ab_shift="slip")


@pytest.mark.parametrize("name, error", [("XY_ecl", "AA, AB_ or ABC"), ("AA_x", "ecl, s1, s2")])
def test_unknown_stacking_types_are_rejected(monolayer, name, error):
    table = StackingTable(modes=(mode(name, 0.0, (0.0, 0.0)),))
    with pytest.raises(ValueError, match=error):
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
    arguments = {"n_layers": 3, "temperature": 293, "ab_shift": "slip", **kwargs}
    with pytest.raises(ValueError, match=error):
        build_stacked_model(monolayer, table, **arguments)


def test_c4_sets_the_other_ab_slip_energy_to_zero(monolayer):
    table = StackingTable(
        modes=(
            mode("AA_ecl", 40.0, (0.9, 0.7)),
            mode("AB_s2_axis", 40.0, (0.9, 0.7)),
            mode("AB_s2_diag", 40.0, (0.9, 0.7)),
        )
    )
    _, records = build_stacked_model(
        monolayer,
        table,
        40,
        293,
        symmetry=4,
        ab_slip="axis",
        ab_shift="slip",
        rng=np.random.RandomState(0),
    )
    assert {r.mode for r in records} == {"AB_s2_diag"}
    assert table.energies()[2] == 40.0  # the input table is unchanged
