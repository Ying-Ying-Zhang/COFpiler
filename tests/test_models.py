import numpy as np
import pandas as pd
import pytest
from conftest import EXAMPLE
from pydantic import ValidationError

from cofpiler import StackingMode, StackingTable, ureg


def mode(name="AA_ecl", energy="1 kJ/mol", shift=None):
    shift = ureg.Quantity([0.0, 0.0, 3.4], "angstrom") if shift is None else shift
    return StackingMode(name=name, energy=energy, shift=shift)


def test_reads_example_table():
    table = StackingTable.read(EXAMPLE / "example.xlsx")
    assert len(table.modes) == 11
    assert table.names[0] == "AA_ecl"
    assert table.energies()[0] == pytest.approx(25.280678)
    np.testing.assert_allclose(table.shifts()[-1], [10.6, 1.2, 2.8])


def test_csv_and_excel_give_the_same_table(tmp_path):
    frame = pd.read_excel(EXAMPLE / "example.xlsx", sheet_name="data")
    frame.to_csv(tmp_path / "example.csv", index=False)
    from_csv = StackingTable.read(tmp_path / "example.csv")
    from_excel = StackingTable.read(EXAMPLE / "example.xlsx")
    assert from_csv.names == from_excel.names
    np.testing.assert_allclose(from_csv.energies(), from_excel.energies(), rtol=1e-12)
    np.testing.assert_allclose(from_csv.shifts(), from_excel.shifts(), rtol=1e-12)


def test_energy_per_unit_is_converted_to_per_mole():
    assert mode(energy="1 eV").energy.to("kJ/mol").magnitude == pytest.approx(96.485, rel=1e-4)


def test_shift_in_nanometre_is_converted_to_angstrom():
    m = mode(shift=ureg.Quantity([0.1, 0.0, 0.34], "nm"))
    np.testing.assert_allclose(m.shift.magnitude, [1.0, 0.0, 3.4])


def test_table_unit_options_convert_columns(tmp_path):
    pd.DataFrame(
        {"stacking_type": ["AA_ecl"], "Erel": [0.01], "x": [0.0], "y": [0.0], "z": [0.34]}
    ).to_csv(tmp_path / "t.csv", index=False)
    table = StackingTable.read(tmp_path / "t.csv", energy_unit="eV", length_unit="nm")
    assert table.energies()[0] == pytest.approx(0.96485, rel=1e-4)
    np.testing.assert_allclose(table.shifts()[0], [0.0, 0.0, 3.4])


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"energy": "1 angstrom"}, "expected an energy"),
        ({"energy": 1.0}, "quantity with units"),
        ({"energy": "nan kJ/mol"}, "finite"),
        ({"shift": ureg.Quantity([0.0, 3.4], "angstrom")}, "x, y and z"),
        ({"shift": ureg.Quantity([0.0, 0.0, 3.4], "kJ/mol")}, "expected a length"),
        ({"shift": ureg.Quantity([0.0, 0.0, -3.4], "angstrom")}, "must be positive"),
        ({"name": ""}, "at least 1 character"),
    ],
)
def test_invalid_modes_are_rejected(kwargs, message):
    with pytest.raises(ValidationError, match=message):
        mode(**kwargs)


def test_duplicate_names_are_rejected():
    with pytest.raises(ValidationError, match="duplicate stacking types: AA_ecl"):
        StackingTable(modes=(mode(), mode()))


def test_empty_table_is_rejected():
    with pytest.raises(ValidationError):
        StackingTable(modes=())


def test_missing_columns_are_reported():
    with pytest.raises(ValueError, match="missing columns: z"):
        StackingTable.from_frame(pd.DataFrame(columns=["stacking_type", "Erel", "x", "y"]))


def test_temperature_with_units():
    table = StackingTable(modes=(mode("AA_ecl", "0 kJ/mol"), mode("AA_s1", "1 kJ/mol")))
    np.testing.assert_allclose(
        table.probabilities(ureg.Quantity(20.0, "degC")), table.probabilities(293.15)
    )
