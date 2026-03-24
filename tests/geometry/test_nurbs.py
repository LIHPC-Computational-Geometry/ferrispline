import pytest
import numpy as np

from nurbs_math.geometry.nurbs import eval_nurbs_curve, eval_bspline

def test_eval_bspline_degree_zero():
    """Test that degree 0 basis function returns 1.0 inside the interval and 0.0 outside."""
    knots = [0.0, 1.0, 2.0, 3.0]
    
    # Inside the first interval [0.0, 1.0[
    assert eval_bspline(0, 0, knots, 0.5) == 1.0
    # Outside the first interval
    assert eval_bspline(0, 0, knots, 1.5) == 0.0

def test_eval_bspline_local_support():
    """Test that the basis function is 0 outside its local support domain."""
    knots = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0]
    degree = 2
    
    # For i=1, the support domain is strictly between knots[1] and knots[1 + 2 + 1] -> [0.0, 2.0[
    # Evaluating at u=2.5 should yield strictly 0.0
    assert eval_bspline(1, degree, knots, 2.5) == 0.0

def test_eval_bspline_partition_of_unity():
    """Test that the sum of all basis functions of degree p at a valid u equals 1.0."""
    knots = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0]
    degree = 2
    u = 1.5  # A valid parameter inside the domain
    
    # Number of control points = len(knots) - degree - 1 = 8 - 2 - 1 = 5
    num_ctrl_pts = 5 
    
    total_sum = sum(eval_bspline(i, degree, knots, u) for i in range(num_ctrl_pts))
    np.testing.assert_allclose(total_sum, 1.0)

@pytest.mark.parametrize(
    "i, degree, knots, expected_error",
    [
        (0, -1, [0.0, 1.0, 2.0, 3.0], "Degree cannot be negative"),
        (-1, 1, [0.0, 1.0, 2.0, 3.0], "is out of bounds for knot vector"),
        (2, 1, [0.0, 1.0, 2.0], "is out of bounds for knot vector"),
    ],
    ids=["negative_degree", "negative_index", "index_too_large"]
)
def test_eval_bspline_invalid_arguments(i, degree, knots, expected_error):
    """Test thata ValueError is raised if argument of eval_bspline are not valid."""
    
    u_constant = 0.5
    
    with pytest.raises(ValueError, match=expected_error):
        evalBspline(i, degree, knots, u_constant)

def test_eval_nurbs_curve_straight_line():
    """Test that a degree 1 NURBS curve evaluates to a straight line between control points."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1
    sample = 10

    curve = eval_nurbs_curve(knots, control_points, weights, degree, sample)

    assert curve.shape == (10, 3)
    np.testing.assert_allclose(curve[0], control_points[0], err_msg="Curve should start at the first control point")
    np.testing.assert_allclose(curve[-1], control_points[-1], err_msg="Curve should end at the last control point")

def test_eval_nurbs_curve_weight_invariance():
    """Test that scaling all weights by the same factor does not change the resulting curve."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]])
    knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    degree = 2
    sample = 10

    weights_base = np.array([1.0, 1.0, 1.0])
    curve_base = eval_nurbs_curve(knots, control_points, weights_base, degree, sample)

    weights_scaled = np.array([2.0, 2.0, 2.0])
    curve_scaled = eval_nurbs_curve(knots, control_points, weights_scaled, degree, sample)

    np.testing.assert_allclose(curve_base, curve_scaled, err_msg="Uniformly scaled weights should produce the exact same curve")

def test_eval_nurbs_curve_inconsistent_sizes():
    """Test that a ValueError is raised if the arrays have inconsistent sizes."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1

    with pytest.raises(ValueError) as exc_info:
        eval_nurbs_curve(knots, control_points, weights, degree)
    assert "Controle point and its weights are differents size" in str(exc_info.value)

def test_eval_nurbs_curve_degree_out_of_bounds():
    """Test that a ValueError is raised if the degree is physically impossible for the knot vector."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 1.0]
    degree = 5

    with pytest.raises(ValueError) as exc_info:
        eval_nurbs_curve(knots, control_points, weights, degree)
    assert "The degree is physically impossible" in str(exc_info.value)

def test_eval_nurbs_curve_2d_points_zero_denominator():
    """Test that the function handles 2D control points properly, even if denominator is 0."""
    control_points = np.array([[0.0, 0.0], [1.0, 1.0]])
    weights = np.array([0.0, 0.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1
    sample = 2

    curve = eval_nurbs_curve(knots, control_points, weights, degree, sample)
    
    assert curve.shape == (2, 2)
    np.testing.assert_array_equal(curve[0], [0.0, 0.0])

def test_eval_nurbs_curve_invalid_sample():
    """Test that a ValueError is raised if sample size is zero or negative."""
    control_points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    weights = np.array([1.0, 1.0])
    knots = [0.0, 0.0, 1.0, 1.0]
    degree = 1

    with pytest.raises(ValueError) as exc_info:
        eval_nurbs_curve(knots, control_points, weights, degree, sample=-5)
    assert "Sample size can not be zero or negative" in str(exc_info.value)