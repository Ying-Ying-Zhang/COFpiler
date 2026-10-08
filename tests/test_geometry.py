import numpy as np
import pytest

from cofpiler.geometry import reflect_xy, rotate_xy, wrap_c_vector


def test_rotation_is_clockwise():
    np.testing.assert_allclose(rotate_xy([1.0, 0.0], np.pi / 2), [0.0, -1.0], atol=1e-12)


def test_full_turn_is_identity():
    np.testing.assert_allclose(rotate_xy([0.3, 1.7], 2 * np.pi), [0.3, 1.7], atol=1e-12)


def test_reflection_keeps_vectors_on_the_line():
    np.testing.assert_allclose(reflect_xy([2.0, 2.0], (1.0, 1.0)), [2.0, 2.0])


def test_reflection_flips_the_normal_component():
    np.testing.assert_allclose(reflect_xy([1.0, -1.0], (1.0, 1.0)), [-1.0, 1.0])


def test_reflection_solves_the_2022_equations():
    # The 2022 script found the mirror image (x, y) of (x0, y0) by solving these with fsolve.
    x0, y0 = 0.9, 0.7
    k = 1 / 0.001
    x, y = reflect_xy([x0, y0], (0.001, 1.0))
    assert y0 + y - k * (x0 + x) == pytest.approx(0, abs=1e-9)
    assert k * (y - y0) + x - x0 == pytest.approx(0, abs=1e-9)


A = np.array([17.3, 0.0, 0.0])
B = np.array([-8.65, 14.98, 0.0])


def test_short_c_vector_is_unchanged():
    np.testing.assert_allclose(wrap_c_vector([1.0, 2.0, 20.0], A, B), [1.0, 2.0, 20.0])


def test_long_in_plane_components_are_wrapped():
    # Two steps of -b bring y to 1.0 and x to 17.3, one step of -a brings x to 0.
    c = wrap_c_vector([0.0, 2 * 14.98 + 1.0, 20.0], A, B)
    np.testing.assert_allclose(c, [0.0, 1.0, 20.0], atol=1e-12)


def test_flat_vector_is_tilted_to_positive_x_and_y():
    # +b makes y positive but x negative, then +a makes x positive.
    np.testing.assert_allclose(wrap_c_vector([0.0, -10.0, 3.4], A, B), [8.65, 4.98, 3.4])


def test_degenerate_cell_is_rejected():
    with pytest.raises(ValueError, match="cell"):
        wrap_c_vector([0.0, 1.0, 3.4], A, [1.0, 0.0, 0.0])
