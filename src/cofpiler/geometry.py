"""In-plane vector operations and the stacking-direction lattice vector."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

Vector = NDArray[np.float64]


def rotate_xy(vector: ArrayLike, angle: float) -> Vector:
    """Rotate a 2D vector clockwise by `angle` (radians)."""
    c, s = math.cos(angle), math.sin(angle)
    return np.dot(np.array([[c, s], [-s, c]]), np.asarray(vector, dtype=float))


def reflect_xy(vector: ArrayLike, mirror_plane: Sequence[float]) -> Vector:
    """Reflect a 2D vector across the line through the origin along `mirror_plane`."""
    phi = np.arctan2(mirror_plane[1], mirror_plane[0])
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    return np.array([[c, s], [s, -c]]) @ np.asarray(vector, dtype=float)


def wrap_c_vector(c: ArrayLike, a: ArrayLike, b: ArrayLike) -> Vector:
    """Bring the stacking vector `c` close to the cell axis using the in-plane vectors `a` and `b`.

    Shifts by whole lattice vectors until the in-plane components are within 0.8 of the
    cell length, then tilts toward positive x and y when `c` is flatter than 50 degrees.
    """
    c = np.array(c, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a[0] == 0 or b[1] == 0:
        raise ValueError("the cell needs a[0] != 0 and b[1] != 0")

    while np.abs(c[1]) > 0.8 * np.abs(b[1]):
        c -= np.sign(c[1] * b[1]) * b
    while np.abs(c[0]) > 0.8 * np.abs(a[0]):
        c -= np.sign(c[0] * a[0]) * a
    if c[1] < 0 and np.arctan(np.abs(c[2] / c[1])) < np.radians(50):
        c += b
    if c[0] < 0 and np.arctan(np.abs(c[2] / c[0])) < np.radians(50):
        c += a
    return c
