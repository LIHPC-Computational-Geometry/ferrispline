import pytest
import numpy as np

from nurbs_math.geometry.nurbs import evalNURBSCurve, evalBspline

def test_evalNURBSCurve_straight_line():
    """Test that a degree 1 NURBS curve evaluates to a straight line between control points."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1
    sample = 10

    curve = evalNURBSCurve(knots, control_points, weights, degree, sample)

    assert curve.shape == (10, 3)
    np.testing.assert_allclose(curve[0], control_points[0], err_msg="Curve should start at the first control point")
    np.testing.assert_allclose(curve[-1], control_points[-1], err_msg="Curve should end at the last control point")


def test_evalNURBSCurve_weight_invariance():
    """Test that scaling all weights by the same factor does not change the resulting curve."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]])
    knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    degree = 2
    sample = 10

    weights_base = np.array([1.0, 1.0, 1.0])
    curve_base = evalNURBSCurve(knots, control_points, weights_base, degree, sample)

    weights_scaled = np.array([2.0, 2.0, 2.0])
    curve_scaled = evalNURBSCurve(knots, control_points, weights_scaled, degree, sample)

    np.testing.assert_allclose(curve_base, curve_scaled, err_msg="Uniformly scaled weights should produce the exact same curve")


def test_evalNURBSCurve_inconsistent_sizes():
    """Test that a ValueError is raised if the arrays have inconsistent sizes."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1

    with pytest.raises(ValueError) as exc_info:
        evalNURBSCurve(knots, control_points, weights, degree)
    assert "Controle point and its weights are differents size" in str(exc_info.value)

def test_evalNURBSCurve_degree_out_of_bounds():
    """Test that a ValueError is raised if the degree is physically impossible for the knot vector."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 1.0]
    degree = 5

    with pytest.raises(ValueError) as exc_info:
        evalNURBSCurve(knots, control_points, weights, degree)
    assert "The degree is physically impossible" in str(exc_info.value)

def test_evalNURBSCurve_2d_points_zero_denominator():
    """Test that the function handles 2D control points properly, even if denominator is 0."""
    control_points = np.array([[0.0, 0.0], [1.0, 1.0]])
    weights = np.array([0.0, 0.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1
    sample = 2

    curve = evalNURBSCurve(knots, control_points, weights, degree, sample)
    
    assert curve.shape == (2, 2)
    np.testing.assert_array_equal(curve[0], [0.0, 0.0])

def test_evalNURBSCurve_invalid_sample():
    """Test that a ValueError is raised if sample size is zero or negative."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1

    with pytest.raises(ValueError) as exc_info:
        evalNURBSCurve(knots, control_points, weights, degree, sample=-5)
    assert "Sample size can not be zero or negative" in str(exc_info.value)