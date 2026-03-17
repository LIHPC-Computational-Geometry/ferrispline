import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
import matplotlib.cm as cm
import sys
from enum import IntEnum

class Multiplicity3Degrees(IntEnum):
    SIMPLE = 1
    REDUSED = 2
    FULL = 3
    CLAMPED = 4

globals().update(Multiplicity3Degrees.__members__)

def buildKnotVector(knot_definitions):
    knot_vector = []
    for value, multiplicity in knot_definitions:
        knot_vector.extend([value] * multiplicity)
    return knot_vector

# NOTE segment_index is increment number `i` of the `figure` function with range (3, 10)
def computeKnotInsertionMatrix(knots, degree: int, segment_index: int):
    r"""Creation of a knot insertion matrix using Boehm's knot insertion algorithm.

        This conversion matrix allows the transformation of a local segment of a 
        B-Spline (or NURBS) curve into an independent Bezier curve.

        Mathematical Mechanics (Pyramid Construction):
        ----------------------------------------------
        The construction of this matrix is not done in a single step. It uses an 
        iterative and cumulative process. To obtain the final matrix for a curve 
        of degree p, the algorithm must sequentially calculate the matrices for the 
        underlying degrees (degree 1, then 2... up to p). Each new level of the 
        pyramid uses the results from the previous level.

        Detailed Calculations per Step:
        -------------------------------
        1. Interpolation ratios (alpha and beta): 
        The algorithm calculates the blending percentages based on the temporal 
        distance between the knots. `alpha` handles the influence for the iterative 
        insertion of the knot marking the start of the segment (`knot[segment_index]`), 
        while `beta` handles the insertion of the knot marking the end of the 
        segment (`knot[segment_index + 1]`).

        2. Temporary matrices (A and B): 
        These ratios populate two transition matrices:
        - Matrix A applies the insertion for the lower bound `knot[segment_index]`, 
            forcing the curve to anchor at the beginning of the segment.
        - Matrix B applies the insertion for the upper bound `knot[segment_index + 1]`, 
            forcing the curve to terminate cleanly at the end of the segment.

        3. Assembly and Accumulation: 
        The global extraction matrix grows at each iteration. It is multiplied 
        by matrix A (for the upper part), and its last row is multiplied by 
        matrix B (for the lower part). Both halves are then stacked vertically 
        (`np.vstack`) to form the calculation base for the next step.
        
        At the final pass (degree p), the resulting matrix will multiply the local 
        B-Spline control points to output the exact Bezier control points for this segment.
        
        Args:
            knot (array_like): Vector of knots.
            degree (int): Curve's degree.
            segment_index (int): Indicates where we want to extract the segment 
                (between knot[segment_index] and knot[segment_index + 1]).

        Returns:
            knot_insertion_matrix (array_like): Matrix of size (degree + 1) x (degree + 1).
    """

    extraction_matrix = np.eye(1) # NOTE 2D Identity Matrice creation 1*1

    for degree_step in range(1, degree + 1):
        start_idx = max(0, segment_index - degree_step)
        end_idx = min(len(knots), segment_index + degree_step + 2)
        local_knots = knots[start_idx:end_idx]

        tmp_matrix_A = np.zeros((degree_step, degree_step + 1))
        tmp_matrix_B = np.zeros((degree_step, degree_step + 1))

        for row in range(degree_step):
            knot_start = row + 1
            knot_end = knot_start + degree_step
            
            # NOTE: The total distance between our two anchor nodes
            distance = (local_knots[knot_end] - local_knots[knot_start])

            if distance != 0:
                alpha = (local_knots[degree_step] - local_knots[knot_start]) / distance
                beta = (local_knots[degree_step + 1] - local_knots[knot_start]) / distance
            else:
                alpha = 0.0
                beta = 0.0

            tmp_matrix_A[row, row] = 1 - alpha 
            tmp_matrix_A[row, row + 1] = alpha      
            
            tmp_matrix_B[row, row] = 1 - beta  
            tmp_matrix_B[row, row + 1] = beta

        upper_half = extraction_matrix @ tmp_matrix_A  
        lower_half = extraction_matrix[-1:] @ tmp_matrix_B

        extraction_matrix = np.vstack([upper_half, lower_half])

    return extraction_matrix


def bernstein(i, degree, t):
    return comb(degree, i) * pow(t, i) * pow((1 - t), (degree - i))


def rationalBasisFunction(weights, degree, sample=100):
    r""" Calcule the rational basis function

        This function corresponds to the first part of the mathematical formula, 
        calculating the weighted influence of each control point at a time t.

        Args:
            weights (array_like): weights vector of all control points 
            degree (array_like): degree of the curve
            sample (int, optional): number of points evaluate into the bezier curve. Defaults to 100.

        Returns:
            array_like: a matrix of (degree + 1) * sample size
    """

    t = np.linspace(0, 1, sample)
    weighted_strength = np.zeros((degree + 1, sample))
    for i in range(degree + 1):
        force = bernstein(i, degree, t)
        weighted_strength += weights[i] * force
    denominator = np.sum(weighted_strength, axis=0)
    return weighted_strength / denominator



def evalBezierCurve(control_points, weights, degree, sample=100):
    r"""
        Evaluates rational Bezier curve and returns it.

        The formula is:
        .. math::
            \text{curve}(t) = \frac{1}{\sum_{i=0}^{\text{degree}} \text{weights}[i] \start_idx{pmatrix} \text{degree} \\ i \end_idx{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)}} \sum_{i=0}^{\text{degree}} \text{weights}[i] \start_idx{pmatrix} \text{degree} \\ i \end_idx{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)} \text{control_points}[i]
 
        Args:
        control_points (array_like): Control points vector.
        weights (array_like): Weights vector.
        degree (int): Bezier basis degree.
        sample (int, optional): Render sample.

        Returns:
            curve (array_like): Rational Bezier curve.
    """

    rational_basis = rationalBasisFunction(weights, degree)
    transposed_rational_basis = rational_basis.T
    curve_points = transposed_rational_basis @ control_points
    return curve_points

def evalNURBSCurve(nodes, control_points, weights, degree, sample=300):
    u_min = nodes[degree]
    u_max = nodes[-degree - 1]
    u_vals = np.linspace(u_min, u_max, sample)
    curve = np.zeros((sample, control_points.shape[1]))

    for idx, u in enumerate(u_vals):
        numerator = np.zeros(control_points.shape[1])
        denominator = 0.0
        for i in range(len(control_points)):
            N = evalBspline(i, degree, nodes, u)
            numerator += weights[i] * N * control_points[i]
            denominator += weights[i] * N
        curve[idx] = numerator / denominator
    return curve

def bezierSegment():
    knot_insertion_matrix = computeKnotInsertionMatrix(nodes, degree, i)
    return bezier_segment

def evalNURBSSurface(
    nodes_u,
    nodes_v,
    control_points,
    weights,
    degree_u,
    degree_v,
    nb_points_u=50,
    nb_points_v=50,
):
    u_min = nodes_u[degree_u]
    u_max = nodes_u[-degree_u - 1]
    u_vals = np.linspace(u_min, u_max, nb_points_u)

    v_min = nodes_v[degree_v]
    v_max = nodes_v[-degree_v - 1]
    v_vals = np.linspace(v_min, v_max, nb_points_v)

    surface = np.zeros((nb_points_u, nb_points_v, 3))

    for iu, u in enumerate(u_vals):
        for iv, v in enumerate(v_vals):
            numerator = np.zeros(3)
            denominator = 0.0
            for i in range(control_points.shape[0]):
                Ni = evalBsplineSurface(i, degree_u, nodes_u, u)
                for j in range(control_points.shape[1]):
                    Mj = evalBsplineSurface(j, degree_v, nodes_v, v)
                    weights_ij = weights[i, j]
                    NMi_w = Ni * Mj * weights_ij
                    numerator += NMi_w * control_points[i, j]
                    denominator += NMi_w
            surface[iu, iv] = numerator / denominator if denominator != 0 else numerator
    return surface


def evalBspline(i, degree, nodes, u):
    n = len(nodes) - 1
    if degree == 0:
        if i >= n:
            return 0.0
        return 1.0 if nodes[i] <= u < nodes[i + 1] else 0.0
    first_part = 0.0
    second_part = 0.0
    if (i + degree) < n:
        denom1 = nodes[i + degree] - nodes[i]
        if denom1 != 0:
            first_part = (u - nodes[i]) / denom1 * evalBspline(i, degree - 1, nodes, u)
    if (i + degree + 1) < n:
        denom2 = nodes[i + degree + 1] - nodes[i + 1]
        if denom2 != 0:
            second_part = (
                (nodes[i + degree + 1] - u)
                / denom2
                * evalBspline(i + 1, degree - 1, nodes, u)
            )
    return first_part + second_part


# we add this specifically for surfaces
def evalBsplineSurface(i, degree, nodes, parameter):
    n = len(nodes) - 1
    if degree == 0:
        if i >= n:
            return 0.0
        return 1.0 if nodes[i] <= parameter < nodes[i + 1] else 0.0
    first_part = 0.0
    second_part = 0.0
    if (i + degree) < n:
        denom1 = nodes[i + degree] - nodes[i]
        if denom1 != 0:
            first_part = (
                (parameter - nodes[i])
                / denom1
                * evalBsplineSurface(i, degree - 1, nodes, parameter)
            )
    if (i + degree + 1) < n:
        denom2 = nodes[i + degree + 1] - nodes[i + 1]
        if denom2 != 0:
            second_part = (
                (nodes[i + degree + 1] - parameter)
                / denom2
                * evalBsplineSurface(i + 1, degree - 1, nodes, parameter)
            )
    return first_part + second_part


def figure(degree: int, knots, control_points, ctrl_pt_weights):
    bezier_segments = []

    for i in range(degree, len(knots) - degree - 1):
        if knots[i] == knots[i + 1]:
            continue
    
        ctrl_pt_start_idx = i - degree
        ctrl_pt_end_idx = i
    
        if first < 0 or last >= len(control_points):
            continue

        knot_insertion_matrix = computeKnotInsertionMatrix(knots, degree, i)
        
        local_points = control_points[ctrl_pt_start_idx : ctrl_pt_end_idx + 1]
        local_weights = ctrl_pt_weights[ctrl_pt_start_idx : ctrl_pt_end_idx + 1]

        weighted_points = local_weights[:, np.newaxis] * local_points
        
        bezier_weighted_points = knot_insertion_matrix @ weighted_points
        bezier_weights = knot_insertion_matrix @ local_weights

        bezier_points = bezier_weighted_points / bezier_weights[:, np.newaxis]

        curve = evalBezierCurve(bezier_points, bezier_weights, degree)
        bezier_segments.append(curve)
    nurbs_curve = evalNURBSCurve(knots, control_points, weights, degree)


    # DRAW
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")
    colors = cm.get_cmap("tab10", len(bezier_segments))

    for idx, segment in enumerate(bezier_segments):
        ax.plot(
            segment[:, 0],
            segment[:, 1],
            segment[:, 2],
            color=colors(idx),
            label=f"Edge {idx+1}",
        )

    # NURBS
    ax.plot(
        nurbs_curve[:, 0],
        nurbs_curve[:, 1],
        nurbs_curve[:, 2],
        "k--",
        linewidth=2,
        label="NURBS curve",
    )

    # Control points
    ax.plot(
        control_points[:, 0],
        control_points[:, 1],
        control_points[:, 2],
        "ks--",
        label="Control points",
    )

    # u_vals = np.linspace(knots[degree], knots[-degree - 1], 1000)
    x_vals = []
    y_vals = []
    z_vals = []
    for u in knots[degree:-degree]:
        numerator = np.zeros(control_points.shape[1])
        denominator = 0.0
        for i in range(len(control_points)):
            N = evalBspline(i, degree, knots, u)
            numerator += weights[i] * N * control_points[i]
            denominator += weights[i] * N
        point = numerator / denominator
        x_vals.append(point[0])
        y_vals.append(point[1])
        z_vals.append(point[2])

    ax.scatter(x_vals, y_vals, z_vals, color="red", s=50, label="Bezier points")

    ax.set_title("NURBS to Bezier")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    ax.grid(True)

    return fig


def article():
    degree = 3
    knots_definitions = [
        (0.0, CLAMPED),
        (1/5, SIMPLE),
        (2/5, REDUSED),
        (3/5, FULL),
        (1, CLAMPED)
    ]

    knots = buildKnotVector(knots_definitions)

    control_points = np.array(
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

    ctrl_pt_weights = np.array([1, 2, 2, 1, 0.5, 0.5, 1, 1, 2, 1])

    return figure(
        knots=knots, degree=degree, control_points=control_points, ctrl_pt_weights=ctrl_pt_weights
    )


if __name__ == "__main__":
    article()
    plt.show()
