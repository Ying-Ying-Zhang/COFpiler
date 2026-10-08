import numpy as np
import pytest

from cofpiler.boltzmann import INVERSE_GAS_CONSTANT, boltzmann_probabilities


def test_probabilities_sum_to_exactly_one():
    p = boltzmann_probabilities([25.28, 21.54, 21.60, 40.31, 0.0, 1.43], 293)
    assert p.sum() == 1.0


def test_ratio_follows_boltzmann_factor():
    p = boltzmann_probabilities([0.0, 2.0], 300)
    assert p[1] / p[0] == pytest.approx(np.exp(-2.0 * INVERSE_GAS_CONSTANT / 300))


def test_equal_energies_are_equally_likely():
    np.testing.assert_allclose(boltzmann_probabilities([5.0, 5.0, 5.0, 5.0], 293), 0.25)


def test_lowest_energy_is_most_likely():
    assert np.argmax(boltzmann_probabilities([3.0, 0.5, 1.0], 293)) == 1


@pytest.mark.parametrize("temperature", [0, -10])
def test_rejects_non_positive_temperature(temperature):
    with pytest.raises(ValueError, match="temperature"):
        boltzmann_probabilities([0.0, 1.0], temperature)
