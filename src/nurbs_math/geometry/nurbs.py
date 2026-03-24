import numpy as np

from ..core_types import MatrixMxN, MatrixNx3, MatrixNxN, Vector, Vector3

def evalNURBSCurve(knots: list, control_points: MatrixNx3, ctrl_pt_weights: Vector, degree: int, sample: int=300) -> MatrixNx3:
    if sample <= 0:
        raise ValueError("Sample size can not be zero or negative")
    if len(control_points) != len(ctrl_pt_weights):
        raise ValueError("Controle point and its weights are differents size")
    if degree != len(knots) - len(control_points) - 1:
        print({degree}, {len(knots)}, {len(control_points)}, {len(knots) - len(control_points) - 1})
        raise ValueError("The degree is physically impossible")

    u_min: int = knots[degree] # NOTE: Start of the valid parameter domain (ensures partition of unity)
    u_max: int = knots[-degree - 1] # NOTE: End of the valid parameter domain
    u_vals: Vector = np.linspace(u_min, u_max, sample)
    curve: MatrixNx3 = np.zeros((sample, control_points.shape[1]))

    # NOTE: Pour tous points 'u' entre 'u_min' et 'u_max',
    # trouve la position sur la courbe grace à la continuité des 'knots' et la force d'attraction des points de control
    for idx, u in enumerate(u_vals):
        numerator: Vector3 = np.zeros(control_points.shape[1])
        denominator: float = 0.0
        for i in range(len(control_points)):
            N: float = evalBspline(i, degree, knots, u)
            numerator += ctrl_pt_weights[i] * N * control_points[i]
            denominator += ctrl_pt_weights[i] * N
         
        curve[idx] = np.zeros(control_points.shape[1]) if denominator == 0 else numerator / denominator
    return curve


def evalBspline(i: int, degree: int, knots: list, u: float) -> float:
    n: int = len(knots) - 1

    if degree < 0:
        raise ValueError(f"Degree cannot be negative. Received: {degree}")
    if i < 0 or i >= n:
        raise ValueError(f"Index i ({i}) is out of bounds for knot vector of length {len(knots)}")
    if degree == 0:
        if i < n and knots[i] <= u < knots[i + 1]:
            return 1.0
        elif i < n and knots[i] <= u <= knots[i + 1] and u == knots[-1]:
            return 1.0
        else:
            return 0.0
    
    first_part: float = 0.0
    second_part: float  = 0.0

    if (i + degree) < n:
        denom1 = knots[i + degree] - knots[i]
        if denom1 != 0:
            first_part = (u - knots[i]) / denom1 * evalBspline(i, degree - 1, knots, u)

    if (i + degree + 1) < n:
        denom2 = knots[i + degree + 1] - knots[i + 1]
        if denom2 != 0:
            second_part = ((knots[i + degree + 1] - u) / denom2 * evalBspline(i + 1, degree - 1, knots, u))

    return first_part + second_part
