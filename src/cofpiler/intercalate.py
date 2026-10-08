"""Statistical models with intercalated monomers (port of the 2022 COFpiler_inter.py)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import cast

import numpy as np
from ase import Atoms
from ase.io import read

from cofpiler.geometry import rotate_xy, wrap_c_vector
from cofpiler.models import StackingTable
from cofpiler.stacking import LayerRecord, Rng

N_ROT = 4


def mirror_diagonal(atoms: Atoms) -> Atoms:
    """Swap x and y about the centre of `atoms`, i.e. mirror across the diagonal."""
    atoms = atoms.copy()
    centre = atoms.positions.mean(axis=0)
    x = atoms.positions[:, 1] - centre[1] + centre[0]
    y = atoms.positions[:, 0] - centre[0] + centre[1]
    atoms.positions[:, 0] = x
    atoms.positions[:, 1] = y
    return atoms


def build_intercalated_model(
    directory: str | Path,
    table: StackingTable,
    period: int,
    n_sets: int,
    temperature: float,
    rng: Rng | None = None,
) -> tuple[Atoms, str, list[LayerRecord]]:
    """Stack layers with an intercalated monomer, `n_sets` times.

    For every stacking type `t` in `table`, `directory` holds `t1.xyz` (layers below the
    monomer), `t2.xyz` (the monomer) and `t3.xyz` (layers above). Each set has
    `period + 1` layers, one of them the monomer. Results match the 2022 script,
    including its known issues (see README).

    Returns the structure, the stacking type of the last monomer, and one record per monomer.
    """
    if period < 1 or n_sets < 1:
        raise ValueError(f"period and n_sets must be at least 1, got {period} and {n_sets}")
    rng = np.random.default_rng() if rng is None else rng
    directory = Path(directory)
    cache: dict[str, Atoms] = {}

    def load(mode: str, part: int) -> Atoms:
        name = f"{mode}{part}.xyz"
        if name not in cache:
            cache[name] = cast(Atoms, read(directory / name))
        return cache[name].copy()

    names = table.names
    shifts = table.shifts()
    probabilities = table.probabilities(temperature)
    uniform = np.full(N_ROT, 1 / N_ROT)
    cycle = period + 1

    def random_angle(base: float) -> float:
        return float(base + rng.choice(N_ROT, p=uniform) * 2 * np.pi / N_ROT)

    k = rng.choice(len(names), p=probabilities)
    mode = names[k]
    structure = load(mode, 1)
    cell = structure.cell.array
    shift = shifts[k].copy()
    total = np.zeros(3)
    angle_sum = 0.0
    m = 0
    records: list[LayerRecord] = []

    for layer in range(2, n_sets * cycle + 1):
        is_monomer = (layer + 1) % cycle == 0
        if is_monomer:
            k = rng.choice(len(names), p=probabilities)
            mode = names[k]
            shift = shifts[k].copy()
            part = load(mode, 2)
            part.positions[:, :2] -= shift[:2]
            if rng.choice(2, p=(1 / 2, 1 / 2)) == 0:
                # As in 2022: the monomer is mirrored but its shift is not (see README).
                part = mirror_diagonal(part)
        elif layer % cycle == 0:
            part = load(mode, 3)
        else:
            part = load(mode, 1)
            if (layer + period) % cycle != 0:
                m += 1

        # The input structures are bilayers; move them up for the layers below.
        part.positions[:, 2] += 1 / 3 * cell[2, 2] * m

        if is_monomer:
            # As in 2022: two random rotations, the first angle is not kept (see README).
            shift[:2] = rotate_xy(shift[:2], random_angle(angle_sum))
            angle_sum = random_angle(angle_sum)
            shift[:2] = rotate_xy(shift[:2], angle_sum)
            part.positions[:, :2] += shift[:2]
            records.append(LayerRecord(layer, mode, shift.copy(), total.copy()))
        part.positions[:, 2] += shift[2] * math.floor((layer - 1) / cycle)
        structure.extend(part)

        if layer % cycle == 0:
            total += shift + 1 / 3 * (period - 2) * shift

    c = wrap_c_vector(total + np.array([shift[0], shift[1], 0.0]), cell[0], cell[1])
    structure.set_cell([cell[0], cell[1], c], scale_atoms=False)
    return structure, mode, records
