"""Typed, unit-aware input model for stacking-mode tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pint
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from cofpiler.boltzmann import boltzmann_probabilities

ureg = pint.get_application_registry()

ENERGY_UNIT = "kJ/mol"
LENGTH_UNIT = "angstrom"
COLUMNS = ("stacking_type", "Erel", "x", "y", "z")


def _as_quantity(value: Any) -> pint.Quantity:
    if isinstance(value, pint.Quantity):
        return value
    if isinstance(value, str):
        return ureg.Quantity(value)
    raise ValueError(f"expected a quantity with units, got {value!r}")


class StackingMode(BaseModel):
    """One stacking mode: its name, relative energy, and the shift between two layers."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str = Field(min_length=1)
    energy: pint.Quantity
    shift: pint.Quantity

    @field_validator("energy", mode="before")
    @classmethod
    def _energy_per_mole(cls, value: Any) -> pint.Quantity:
        q = _as_quantity(value)
        if q.check("[energy]"):  # per stacking unit, e.g. eV
            q = q * ureg.avogadro_constant
        if not q.check("[energy] / [substance]"):
            raise ValueError(f"expected an energy, got units of {q.units}")
        q = q.to(ENERGY_UNIT)
        if not np.isfinite(q.magnitude):
            raise ValueError("energy must be a finite number")
        return q

    @field_validator("shift", mode="before")
    @classmethod
    def _length_vector(cls, value: Any) -> pint.Quantity:
        q = _as_quantity(value)
        if not q.check("[length]"):
            raise ValueError(f"expected a length, got units of {q.units}")
        xyz = np.asarray(q.to(LENGTH_UNIT).magnitude, dtype=float)
        if xyz.shape != (3,):
            raise ValueError(f"shift needs x, y and z components, got shape {xyz.shape}")
        if not np.all(np.isfinite(xyz)):
            raise ValueError("shift components must be finite numbers")
        if xyz[2] <= 0:
            raise ValueError(
                f"the z component (interlayer distance) must be positive, got {xyz[2]}"
            )
        return ureg.Quantity(xyz, LENGTH_UNIT)


class StackingTable(BaseModel):
    """The stacking modes a statistical model samples from."""

    model_config = ConfigDict(frozen=True)

    modes: tuple[StackingMode, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def _unique_names(self) -> StackingTable:
        names = self.names
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"duplicate stacking types: {', '.join(duplicates)}")
        return self

    @classmethod
    def from_frame(
        cls,
        frame: pd.DataFrame,
        energy_unit: str = ENERGY_UNIT,
        length_unit: str = LENGTH_UNIT,
    ) -> StackingTable:
        """Build a table from columns stacking_type, Erel, x, y, z."""
        missing = [c for c in COLUMNS if c not in frame.columns]
        if missing:
            raise ValueError(f"missing columns: {', '.join(missing)}")
        return cls(
            modes=tuple(
                StackingMode(
                    name=str(row.stacking_type),
                    energy=ureg.Quantity(float(row.Erel), energy_unit),
                    shift=ureg.Quantity([float(row.x), float(row.y), float(row.z)], length_unit),
                )
                for row in frame.itertuples(index=False)
            )
        )

    @classmethod
    def read(
        cls,
        path: str | Path,
        energy_unit: str = ENERGY_UNIT,
        length_unit: str = LENGTH_UNIT,
        sheet: str = "data",
    ) -> StackingTable:
        """Read a .csv file, or the `sheet` of an Excel file."""
        path = Path(path)
        if path.suffix.lower() == ".csv":
            frame = pd.read_csv(path)
        else:
            frame = pd.read_excel(path, sheet_name=sheet)
        return cls.from_frame(frame, energy_unit, length_unit)

    @property
    def names(self) -> list[str]:
        return [m.name for m in self.modes]

    def energies(self) -> NDArray[np.float64]:
        """Relative energies in kJ/mol."""
        return np.array([m.energy.to(ENERGY_UNIT).magnitude for m in self.modes])

    def shifts(self) -> NDArray[np.float64]:
        """Shift vectors in angstrom, one row per mode."""
        return np.array([m.shift.to(LENGTH_UNIT).magnitude for m in self.modes])

    def probabilities(self, temperature: float | pint.Quantity) -> NDArray[np.float64]:
        """Boltzmann probability of each mode; a plain number is read as kelvin."""
        if isinstance(temperature, pint.Quantity):
            temperature = float(temperature.to("K").magnitude)
        return boltzmann_probabilities(self.energies(), temperature)
