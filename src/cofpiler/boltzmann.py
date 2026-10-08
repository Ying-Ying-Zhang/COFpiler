"""Boltzmann weights of stacking modes."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

# 1/R in K mol/kJ as used since 2022. CODATA gives 120.27; kept for reproducibility.
INVERSE_GAS_CONSTANT = 120.37


def boltzmann_probabilities(energies: ArrayLike, temperature: float) -> NDArray[np.float64]:
    """Probability of each stacking mode.

    Args:
        energies: relative energies in kJ/mol.
        temperature: temperature in K.
    """
    if temperature <= 0:
        raise ValueError(f"temperature must be positive, got {temperature} K")
    weights = np.exp(-np.asarray(energies, dtype=float) * INVERSE_GAS_CONSTANT / temperature)
    probabilities = weights / weights.sum()
    # Random choice needs the probabilities to sum to exactly one.
    probabilities[-1] += 1 - probabilities.sum()
    return probabilities
