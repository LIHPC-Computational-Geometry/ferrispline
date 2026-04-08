import numpy as np
from scipy.special import comb

from ..core_types import MatrixMxN, MatrixNx3, MatrixNxN, VectorN


# NOTE segment_index is increment number `i` of the `figure` function with range (3, 10)
def compute_knot_insertion_matrix(
    knots: list, degree: int, segment_index: int
) -> MatrixNxN:
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
        knot (array_like): VectorN of knots.
        degree (int): Curve's degree.
        segment_index (int): Indicates where we want to extract the segment
            (between knot[segment_index] and knot[segment_index + 1]).

    Returns:
        knot_insertion_matrix (array_like): Matrix of size (degree + 1) x (degree + 1).
    """
    if segment_index < 0 or segment_index >= len(knots) - 1:
        raise ValueError(
            f"segment_index ({segment_index}) is out of bounds for knots of length {len(knots)}"
        )

    extraction_matrix: MatrixNxN = np.eye(1)  # NOTE 2D Identity Matrice creation 1*1

    for degree_step in range(1, degree + 1):
        start_idx: int = max(0, segment_index - degree_step)
        end_idx: int = min(len(knots), segment_index + degree_step + 2)
        local_knots: VectorN = knots[start_idx:end_idx]

        tmp_matrix_A: MatrixNxN = np.zeros((degree_step, degree_step + 1))
        tmp_matrix_B: MatrixNxN = np.zeros((degree_step, degree_step + 1))

        for row in range(degree_step):
            knot_start: int = row + 1
            knot_end: int = knot_start + degree_step

            # NOTE: The total distance between our two anchor nodes
            distance: float = local_knots[knot_end] - local_knots[knot_start]

            if distance != 0:
                # NOTE: ratio d'interpolation d'insertion du noeud de début du segment
                alpha: float = (
                    local_knots[degree_step] - local_knots[knot_start]
                ) / distance
                # NOTE: ratio d'interpolation d'insertion du noeud de fin du segment
                beta: float = (
                    local_knots[degree_step + 1] - local_knots[knot_start]
                ) / distance
            else:
                alpha: float = 0.0
                beta: float = 0.0

            tmp_matrix_A[row, row] = 1 - alpha
            tmp_matrix_A[row, row + 1] = alpha

            tmp_matrix_B[row, row] = 1 - beta
            tmp_matrix_B[row, row + 1] = beta

        upper_half: MatrixNxN = extraction_matrix @ tmp_matrix_A
        lower_half: MatrixNxN = extraction_matrix[-1:] @ tmp_matrix_B

        extraction_matrix = np.vstack([upper_half, lower_half])

    return extraction_matrix


def bernstein(v: int, degree: int, t: VectorN) -> VectorN:
    return comb(degree, v) * pow(t, v) * pow((1 - t), (degree - v))


def rational_basis_bezier_function(
    weights: VectorN, degree: int, sample: int
) -> MatrixMxN:
    r"""Calcule the rational basis function

    This function corresponds to the first part of the mathematical formula,
    calculating the weighted influence of each control point at a time t.

    Args:
        weights (np.ndarray): weights VectorN of all control points
        degree (int): degree of the curve
        sample (int, optional): number of points evaluate into the bezier curve. Defaults to 100.

    Returns:
        np.ndarray: a matrix of (degree + 1) * sample size
    """

    if len(weights) != degree + 1:
        raise ValueError("Length of weights doesn't corresponds to degree + 1")
    t: VectorN = np.linspace(0, 1, sample)
    weighted_strength: MatrixMxN = np.zeros((degree + 1, sample))
    for i in range(degree + 1):
        force: VectorN = bernstein(i, degree, t)
        weighted_strength[i] = weights[i] * force
    denominator: VectorN = np.sum(weighted_strength, axis=0)
    if np.any(denominator == 0):
        raise ValueError("Weighted strength can not be divid by 0")
    return weighted_strength / denominator


def eval_bezier_curve(
    control_points: MatrixNx3, weights: VectorN, degree: int, sample: int = 100
) -> MatrixNx3:
    r"""
    Evaluates rational Bezier curve and returns it.

    The formula is:
    .. math::
        \text{curve}(t) = \frac{1}{\sum_{i=0}^{\text{degree}} \text{weights}[i] \start_idx{pmatrix} \text{degree} \\ i \end_idx{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)}} \sum_{i=0}^{\text{degree}} \text{weights}[i] \start_idx{pmatrix} \text{degree} \\ i \end_idx{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)} \text{control_points}[i]

    Args:
    control_points (array_like): Control points VectorN.
    weights (array_like): Weights VectorN.
    degree (int): Bezier basis degree.
    sample (int, optional): Render sample.

    Returns:
        curve (array_like): Rational Bezier curve.
    """
    if len(control_points) != degree + 1:
        raise ValueError("Length of controle points doesn't corresponds to degree + 1")
    rational_basis: MatrixMxN = rational_basis_bezier_function(weights, degree, sample)
    transposed_rational_basis: MatrixMxN = rational_basis.T
    curve_points: MatrixNx3 = transposed_rational_basis @ control_points
    return curve_points


def bezier_curves(
    knots: list, control_points: MatrixNx3, ctrl_pt_weights: VectorN, degree: int
) -> list:
    bezier_segments: list = []
    for i in range(degree, len(knots) - degree - 1):
        if knots[i] == knots[i + 1]:
            continue

        ctrl_pt_start_idx: int = i - degree
        ctrl_pt_end_idx: int = i

        if ctrl_pt_start_idx < 0 or ctrl_pt_end_idx >= len(control_points):
            continue

        knot_insertion_matrix: MatrixNxN = compute_knot_insertion_matrix(
            knots, degree, i
        )

        # NOTE: correspond aux points et poids de l'intervalle ou l'on souhaite insérer le nouveau point
        local_ctrl_pt: VectorN = control_points[ctrl_pt_start_idx : ctrl_pt_end_idx + 1]
        local_ctrl_pt_weights: VectorN = ctrl_pt_weights[
            ctrl_pt_start_idx : ctrl_pt_end_idx + 1
        ]

        # NOTE: Coordonnées homogènes des points de l'intervalle (matrice contenant n + 1 vecteur de coordonnées spatial multiplié par leur poids np.zeros(n+1, 3))
        weighted_points: MatrixNx3 = (
            local_ctrl_pt_weights[:, np.newaxis] * local_ctrl_pt
        )

        # NOTE: coordonnées homogènes dont l'influence est modifier par la matrice d'insertion (matrice contenant n + 1 vecteur de coordonnées spatial multiplié par leur poids)
        bezier_weighted_points: MatrixNx3 = knot_insertion_matrix @ weighted_points
        # NOTE: Nouveaux poids des points de contrôle de Bézier, calculés par la matrice d'insertion (vecteur 1D contenant n + 1 valeurs scalaires)
        bezier_weights: VectorN = knot_insertion_matrix @ local_ctrl_pt_weights

        if np.any(bezier_weights == 0):
            raise ValueError(
                "Zero weight encountered during Bezier extraction. Cannot divide by zero."
            )
        # NOTE: Supprime l'influence des poids ajoutés artificiellement dans weighted_points par une division, pour obtenir la matrice finale des n+1 vecteurs de coordonnées spatiales 3D (les points de contrôle de Bézier réels).
        bezier_points: MatrixNx3 = (
            bezier_weighted_points / bezier_weights[:, np.newaxis]
        )

        curve: MatrixNx3 = eval_bezier_curve(bezier_points, bezier_weights, degree)
        bezier_segments.append(curve)
    return bezier_segments
