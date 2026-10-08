"""Statistical stacking models of layered materials (port of the 2022 COFpiler.py)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from ase import Atoms
from numpy.typing import ArrayLike, NDArray

from cofpiler.geometry import reflect_xy, rotate_xy, wrap_c_vector
from cofpiler.models import ENERGY_UNIT, StackingTable, ureg

Rng = np.random.Generator | np.random.RandomState

SUPPORTED_SYMMETRIES = (3, 4, 6)
SLIP_KINDS = ("ecl", "s1", "s2")


@dataclass(frozen=True)
class LayerRecord:
    """The stacking type drawn for one layer, its shift, and the accumulated shift."""

    layer: int
    mode: str
    shift: NDArray[np.float64]
    total_shift: NDArray[np.float64]


def ab_lattice_offset(
    cell: ArrayLike, symmetry: int, ab_slip: Literal["axis", "diag"] | None = None
) -> NDArray[np.float64]:
    """In-plane offset between the A and B layers of an AB stacking."""
    a = np.asarray(cell, dtype=float)[0, 0]
    if symmetry == 4:
        if ab_slip == "axis":
            return np.array([0, 1 / 2 * a, 0])
        if ab_slip == "diag":
            return np.array([1 / 2 * a, 1 / 2 * a, 0])
        raise ValueError("symmetry C4 needs ab_slip='axis' or ab_slip='diag'")
    return np.array([0, np.sqrt(3) / 3 * a, 0])


def _check_supported(names: Sequence[str]) -> None:
    for name in names:
        if "AA" not in name and "AB_" not in name:
            raise NotImplementedError(
                f"stacking type {name!r}: only AA and AB_ types are supported. "
                "ABC types have no code path in the 2022 algorithm (see README, Known issues)."
            )
        if not any(kind in name for kind in SLIP_KINDS):
            raise ValueError(f"stacking type {name!r} must contain one of {', '.join(SLIP_KINDS)}")


def _c4_table(table: StackingTable, ab_slip: Literal["axis", "diag"]) -> StackingTable:
    # As in 2022: the AB_s2 mode of the other slip direction gets a relative energy of zero.
    other = {"axis": "AB_s2_diag", "diag": "AB_s2_axis"}[ab_slip]
    if other not in table.names:
        raise ValueError(f"symmetry C4 with ab_slip={ab_slip!r} expects a stacking type {other!r}")
    zero = ureg.Quantity(0.0, ENERGY_UNIT)
    return StackingTable(
        modes=tuple(
            m.model_copy(update={"energy": zero}) if m.name == other else m for m in table.modes
        )
    )


def build_stacked_model(
    monolayer: Atoms,
    table: StackingTable,
    n_layers: int,
    temperature: float,
    symmetry: int = 6,
    mirror: bool = True,
    mirror_plane: Sequence[float] = (0.001, 1.0),
    ab_slip: Literal["axis", "diag"] | None = None,
    rng: Rng | None = None,
) -> tuple[Atoms, list[LayerRecord]]:
    """Stack `n_layers` copies of `monolayer`, drawing each interlayer shift from `table`.

    Each layer draws a stacking type with its Boltzmann probability at `temperature` (K)
    and a random symmetry-equivalent direction for the shift. Results match the 2022
    script, including its known issues (see README).
    """
    if n_layers < 2:
        raise ValueError(f"n_layers must be at least 2, got {n_layers}")
    if symmetry not in SUPPORTED_SYMMETRIES:
        raise ValueError(f"symmetry must be one of {SUPPORTED_SYMMETRIES}, got {symmetry}")
    _check_supported(table.names)
    rng = np.random.default_rng() if rng is None else rng

    cell = monolayer.cell.array
    ab_lattice = ab_lattice_offset(cell, symmetry, ab_slip)
    if symmetry == 4:
        assert ab_slip is not None  # checked by ab_lattice_offset
        table = _c4_table(table, ab_slip)
    names = table.names
    shifts = table.shifts()
    probabilities = table.probabilities(temperature)
    uniform = np.full(symmetry, 1 / symmetry)

    def random_rotation(vector: NDArray[np.float64], base_angle: float) -> NDArray[np.float64]:
        # The step is pi/3 for every symmetry, as in 2022 (see README).
        step = rng.choice(symmetry, p=uniform)
        return rotate_xy(vector, base_angle + step * np.pi / 3)

    structure = monolayer.copy()
    angle = 0.0
    angle_lattice = 0.0
    total = np.zeros(3)
    records: list[LayerRecord] = []

    for layer in range(2, n_layers + 1):
        k = rng.choice(len(names), p=probabilities)
        mode = names[k]
        angle += (rng.choice(symmetry, p=uniform) + 1) * 2 * np.pi / symmetry
        angle_lattice += (rng.choice(symmetry, p=uniform) + 1) * 2 * np.pi / symmetry

        if "AA" in mode:
            shift = shifts[k].copy()
        else:
            shift = shifts[k] - ab_lattice
            sign = 1 if layer % 2 == 0 else -1
            # As in 2022: the rotated lattice offset replaces the in-plane slip (see README).
            shift[:2] = random_rotation(sign * ab_lattice[:2], angle_lattice)

        if "ecl" in mode:
            pass
        elif "s1" in mode:
            shift[:2] = random_rotation(shift[:2], angle)
        elif "s2" in mode:
            if mirror:
                shift[:2] = reflect_xy(shift[:2], mirror_plane)
            shift[:2] = random_rotation(shift[:2], angle)

        if "AB_" in mode:
            shift[:2] += ab_lattice[:2]
        total += shift

        layer_atoms = monolayer.copy()
        layer_atoms.positions += total
        structure.extend(layer_atoms)
        records.append(LayerRecord(layer, mode, shift.copy(), total.copy()))

    if "AB_" in mode:
        c = total - shift
        c[2] = total[2] + shift[2]
    else:
        c = total + shift
    structure.set_cell([cell[0], cell[1], wrap_c_vector(c, cell[0], cell[1])], scale_atoms=False)
    return structure, records
