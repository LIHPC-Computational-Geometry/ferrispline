import sys
from enum import IntEnum

import numpy as np
import pyvista as pv

from .core_types import MatrixNx3, Vector

SIMPLE = 1
REDUSED = 2
FULL = 3
CLAMPED = 4

def load_nurbs_from_vtk(filepath: str, default_degree: int = 3) -> tuple[MatrixNx3, Vector, Vector, int]:
    """ Read a VTK file for extract all data. This function used a custom naming convention.

        -------------------------------------------
        ### Expected VTK ASCII format for a NURBS curve:
        #### vtk DataFile Version 3.0
        ```
        NURBS curve data
        ASCII
        DATASET POLYDATA

        POINTS <N> float
        X0 Y0 Z0
        X1 Y1 Z1
        ...
        XN YN ZN

        FIELD_DATA nurbs_data 1
        knots <K> 1 float
        t0 t1 t2 ... tK

        POINT_DATA <N>
        SCALARS weights float 1
        LOOKUP_TABLE default
        w0
        w1
        ...
        wN
        ```

        #### Important Formatting Rules:
        - `<N>` is the number of control points.
        - `<K>` is the number of knots (Must equal: `N + degree + 1`).
        - FIELD_DATA (`knots`) MUST be declared BEFORE POINT_DATA (`weights`) so the 
        VTK parser reads it as a global array, not as a point-specific attribute.

        --- 
    
        Args:
            filepath (str): path of the vtk file.
            default_degree (int): fallback degree if knots are missing in the file.

        Returns:
            tuple[MatrixNx3, Vector, Vector, int]: a tuple containing:
            - a matrix size N*3 of control points
            - a vector of size N for the weights of these control point
            - a vector of size control_points + degree + 1 with all knots
            - an intager for the curve's degree 
    """
    try:
        mesh = pv.read(filepath)
    except Exception as e:
        raise ValueError(f"Error while reading the file: {e} ")
        
    
    control_points: MatrixNx3 = np.array(mesh.points, dtype=np.float64)
    num_points: int = len(control_points)

    if num_points == 0:
        raise ValueError("Error: The VTK file does not contain any control points.")

    if "weights" not in mesh.point_data:
        ctrl_pt_weights: Vector = np.ones(num_points, dtype=np.float64)
    else:
        ctrl_pt_weights: Vector = np.array(mesh.point_data["weights"], dtype=np.float64)

    if "knots" not in mesh.field_data:
        degree = default_degree

        knots_list = (
            [0.0] * degree +
            list(range(num_points - degree + 1)) +
            [float(num_points - degree)] * degree
        )
        knots: Vector = np.array(knots_list, dtype=np.float64)
    else:
        knots: Vector = np.array(mesh.field_data["knots"], dtype=np.float64).flatten() # NOTE use flatten() because the return of field_data can be a 2D Matrix
        degree: int = len(knots) - len(control_points) - 1

    if degree < 1:
        raise ValueError("Error: the number of control point and knot are invalid: knots number = control point number + degree + 1")
        

    return control_points, ctrl_pt_weights, knots, degree

def build_knot_vector(knot_definitions: list[tuple[float, int]]) -> list:
    knot_vector = []
    for value, multiplicity in knot_definitions:
        if multiplicity <= 0:
            raise ValueError(f"The knot multiplicity {value} must be > 0. Actually multiplicity = {multiplicity}")
        knot_vector.extend([value] * multiplicity)
    return knot_vector

def default_value() -> tuple[MatrixNx3, Vector, Vector, int]:
    degree = 3

    knots_definitions: list[tuple[float, int]] = [
        (0.0, CLAMPED),
        (1/5, SIMPLE),
        (2/5, REDUSED),
        (3/5, FULL),
        (1, CLAMPED)
    ]

    knots: list = build_knot_vector(knots_definitions)

    control_points: MatrixNx3 = np.array(
        [
            [0, 6, 0],
            [1, 10, 0],
            [5, 12, 0],
            [9, 0, 0],
            [8, 3, 0],
            [5, 1.5, 0],
            [0, 0, 0],
            [2, -2, 0],
            [8, -2, 0],
            [10, 0, 0],
        ]
    )

    ctrl_pt_weights: Vector = np.array([1, 2, 2, 1, 0.5, 0.5, 1, 1, 2, 1])
    return control_points, ctrl_pt_weights, knots, degree

