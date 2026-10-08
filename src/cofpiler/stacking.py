"""Statistical stacking models of layered materials (based on the 2022 COFpiler.py)."""

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
Family = Literal["AA", "AB", "ABC"]


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


def stacking_family(name: str) -> Family:
    """Whether a stacking type is AA, AB or ABC, from its name."""
    if "AA" in name:
        return "AA"
    if "ABC" in name:
        return "ABC"
    if "AB_" in name:
        return "AB"
    raise ValueError(f"stacking type {name!r} must contain AA, AB_ or ABC")


def _check_supported(
    names: Sequence[str], symmetry: int, ab_shift: Literal["slip", "full"] | None
) -> list[Family]:
    families = [stacking_family(name) for name in names]
    for name in names:
        if not any(kind in name for kind in SLIP_KINDS):
            raise ValueError(f"stacking type {name!r} must contain one of {', '.join(SLIP_KINDS)}")
    if symmetry == 4 and "ABC" in families:
        raise NotImplementedError("ABC stacking is not implemented for C4 lattices")
    if ab_shift is None and set(families) != {"AA"}:
        raise ValueError(
            "the table has AB or ABC stacking types: set ab_shift (command line: --ab-shift) "
            "to 'slip' if their x, y are the slip on top of the AB position, "
            "or to 'full' if they include the offset to the AB position"
        )
    return families


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
    ab_shift: Literal["slip", "full"] | None = None,
    rng: Rng | None = None,
) -> tuple[Atoms, list[LayerRecord]]:
    """Stack `n_layers` copies of `monolayer`, drawing each interlayer shift from `table`.

    Each layer draws a stacking type with its Boltzmann probability at `temperature` (K)
    and a random symmetry-equivalent direction for the shift. AB and ABC layers also
    shift by the offset to the AB position, turned to a random symmetry-equivalent
    direction. AB flips the offset every other layer and ABC keeps it; since the
    direction is random, both give the same set of interlayer shifts and differ in
    energy, as in AB(C)stat of the 2022 paper.

    `ab_shift` says whether the x, y of AB and ABC rows are the slip on top of the
    offset ("slip") or include it ("full"); it is required when the table has such rows.

    Models of AA layers match the 2022 script, including its known issues (see README).
    """
    if n_layers < 2:
        raise ValueError(f"n_layers must be at least 2, got {n_layers}")
    if symmetry not in SUPPORTED_SYMMETRIES:
        raise ValueError(f"symmetry must be one of {SUPPORTED_SYMMETRIES}, got {symmetry}")
    families = _check_supported(table.names, symmetry, ab_shift)
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
        mode, family = names[k], families[k]
        angle += (rng.choice(symmetry, p=uniform) + 1) * 2 * np.pi / symmetry
        angle_lattice += (rng.choice(symmetry, p=uniform) + 1) * 2 * np.pi / symmetry

        shift = shifts[k].copy()
        offset = np.zeros(2)
        if family != "AA":
            if ab_shift == "full":
                shift[:2] -= ab_lattice[:2]
            # AB goes back and forth between two positions, ABC keeps moving the same way.
            sign = -1 if family == "AB" and layer % 2 == 1 else 1
            offset = random_rotation(sign * ab_lattice[:2], angle_lattice)

        if "ecl" in mode:
            pass
        elif "s1" in mode:
            shift[:2] = random_rotation(shift[:2], angle)
        elif "s2" in mode:
            if mirror:
                shift[:2] = reflect_xy(shift[:2], mirror_plane)
            shift[:2] = random_rotation(shift[:2], angle)

        shift[:2] += offset
        total += shift

        layer_atoms = monolayer.copy()
        layer_atoms.positions += total
        structure.extend(layer_atoms)
        records.append(LayerRecord(layer, mode, shift.copy(), total.copy()))

    if family == "AB":
        c = total - shift
        c[2] = total[2] + shift[2]
    else:
        c = total + shift
    structure.set_cell([cell[0], cell[1], wrap_c_vector(c, cell[0], cell[1])], scale_atoms=False)
    return structure, records
