import pytest

import numpy as np
from unittest.mock import patch, MagicMock

from nurbs_math.load_nurbs import build_knot_vector, load_nurbs_from_vtk


def test_build_knot_vector():
    """Check the creation of knot VectorN with multiplicity"""
    definitions = [(0.0, 3), (0.5, 1), (1.0, 3)]

    result = build_knot_vector(definitions)

    expected = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    assert result == expected


@pytest.mark.parametrize(
    "knot_value, multiplicity",
    [
        (0.5, 0),
        (0.7, -2),
    ],
)
def test_build_knot_vector_multiplicity(knot_value, multiplicity):
    """Teste que les multiplicités de 0 ou négatives lèvent bien une ValueError."""

    invalid_definitions = [(knot_value, multiplicity)]

    with pytest.raises(ValueError) as exc_info:
        build_knot_vector(invalid_definitions)

    assert "must be > 0" in str(exc_info.value)


@patch("nurbs_math.load_nurbs.pv.read")
def test_load_nurbs_from_vtk_success(mock_pv_read):
    """Tests data extraction from a simulated VTK file."""

    mock_mesh = MagicMock()
    mock_mesh.points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
    mock_mesh.point_data = {"weights": np.array([1.0, 2.0, 1.0])}
    mock_mesh.field_data = {"knots": np.array([[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]])}
    mock_pv_read.return_value = mock_mesh

    ctrl_pts, weights, knots, degree = load_nurbs_from_vtk("fake_file.vtk")

    np.testing.assert_array_equal(ctrl_pts, mock_mesh.points)
    np.testing.assert_array_equal(weights, [1.0, 2.0, 1.0])
    np.testing.assert_array_equal(knots, [0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    assert degree == 2


@patch("nurbs_math.load_nurbs.pv.read")
def test_load_nurbs_from_vtk_fallback(mock_pv_read):
    """Test default behavior if weights and knots are absent of VTK file."""

    mock_mesh = MagicMock()
    mock_mesh.points = np.array(
        [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
    )  # 4 points
    mock_mesh.point_data = {}
    mock_mesh.field_data = {}
    mock_pv_read.return_value = mock_mesh

    ctrl_pts, weights, knots, degree = load_nurbs_from_vtk(
        "fake_file.vtk", default_degree=2
    )

    assert degree == 2
    np.testing.assert_array_equal(weights, [1.0, 1.0, 1.0, 1.0])

    expected_knots = [0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0]
    np.testing.assert_array_equal(knots, expected_knots)


@patch("nurbs_math.load_nurbs.pv.read")
def test_load_nurbs_from_vtk_degree_error(mock_pv_read):
    """Tests the secure shutdown of the program (sys.exit) if the calculated degree is < 1."""

    mock_mesh = MagicMock()
    mock_mesh.points = np.array([[0, 0, 0], [1, 0, 0]])
    mock_mesh.point_data = {"weights": np.array([1.0, 1.0])}
    # degree = len(knots) - len(points) - 1 => (2 - 2 - 1) = -1
    mock_mesh.field_data = {"knots": np.array([[0.0, 1.0]])}
    mock_pv_read.return_value = mock_mesh

    with pytest.raises(ValueError) as exc_info:
        load_nurbs_from_vtk("fake_break_file.vtk")

    assert "the number of control point and knot are invalid" in str(exc_info.value)


@patch("nurbs_math.load_nurbs.pv.read")
def test_load_nurbs_from_vtk_read_error(mock_pv_read):
    """Test secure shutdown if PyVista cannot read the file at all."""

    mock_pv_read.side_effect = Exception("Corrupt or non-existent file")

    with pytest.raises(ValueError) as exc_info:
        load_nurbs_from_vtk("corrupted_file.vtk")
    assert "Error while reading the file:" in str(exc_info.value)


@patch("nurbs_math.load_nurbs.pv.read")
def test_load_nurbs_from_vtk_no_points(mock_pv_read):
    """Tests safe shutdown if the file is readable but contains no dots."""

    mock_mesh = MagicMock()
    mock_mesh.points = np.array([])
    mock_pv_read.return_value = mock_mesh

    with pytest.raises(ValueError) as exc_info:
        load_nurbs_from_vtk("fichier_without_points.vtk")

    assert "The VTK file does not contain any control points." in str(exc_info.value)
