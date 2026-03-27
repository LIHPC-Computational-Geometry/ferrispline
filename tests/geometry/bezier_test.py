import pytest
import numpy as np

from nurbs_math.geometry.bezier import (
    bernstein,
    rational_basis_bezier_function,
    eval_bezier_curve,
    compute_knot_insertion_matrix,
    bezier_curves,
)


def test_bernstein_extremes():
    """Test that Bernstein polynomials behave correctly at t=0 and t=1."""
    degree = 3
    # At t=0, only the first polynomial (v=0) is 1.0
    assert bernstein(0, degree, 0.0) == 1.0
    assert bernstein(1, degree, 0.0) == 0.0

    # At t=1, only the last polynomial (v=degree) is 1.0
    assert bernstein(degree, degree, 1.0) == 1.0
    assert bernstein(0, degree, 1.0) == 0.0


def test_bernstein_partition_of_unity():
    """Test that the sum of all Bernstein polynomials for a given degree is always 1.0."""
    degree = 3
    t_vals = np.linspace(0, 1, 10)

    for t in t_vals:
        total = sum(bernstein(v, degree, t) for v in range(degree + 1))
        np.testing.assert_allclose(total, 1.0)


def test_rational_basis_bezier_function_uniform_weights():
    """Test that uniform weights yield standard Bernstein polynomials."""
    degree = 2
    sample = 5
    weights = np.array([1.0, 1.0, 1.0])

    basis = rational_basis_bezier_function(weights, degree, sample)

    assert basis.shape == (degree + 1, sample)
    np.testing.assert_allclose(np.sum(basis, axis=0), np.ones(sample))


def test_eval_bezier_curve_endpoints():
    """Test that a Bezier curve starts and ends exactly on its first and last control points."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]])
    weights = np.array([1.0, 1.0, 1.0])
    degree = 2

    curve = eval_bezier_curve(control_points, weights, degree, sample=10)

    np.testing.assert_allclose(curve[0], control_points[0])
    np.testing.assert_allclose(curve[-1], control_points[-1])


def test_compute_knot_insertion_matrix_dimension():
    """Test that the extraction matrix has the correct (degree+1) x (degree+1) shape."""
    knots = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0]
    degree = 2
    segment_index = 3

    matrix = compute_knot_insertion_matrix(knots, degree, segment_index)
    assert matrix.shape == (degree + 1, degree + 1)


def test_compute_knot_insertion_matrix_invalid_segment_index():
    """Test that a ValueError is raised if the segment index is out of bounds."""
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1

    with pytest.raises(ValueError) as exc_info:
        compute_knot_insertion_matrix(knots, degree, segment_index=-1)
    error_msg = str(exc_info.value)
    assert (
        "segment_index" in error_msg
        and "is out of bounds for knots of length" in error_msg
    )


def test_bezier_curves_zero_weight_division():
    """Test that the function raises a ValueError if division by zero occurs during weight calculation."""
    knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]])
    weights = np.array([1.0, 0.0, 1.0])
    degree = 2

    with pytest.raises(ValueError, match="Zero weight encountered") as exc_info:
        bezier_curves(knots, control_points, weights, degree)
    assert "Cannot divide by zero." in str(exc_info.value)
