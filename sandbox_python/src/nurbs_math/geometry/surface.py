import numpy as np


def eval_nurbs_surface(
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
                Ni = eval_bspline_surface(i, degree_u, nodes_u, u)
                for j in range(control_points.shape[1]):
                    Mj = eval_bspline_surface(j, degree_v, nodes_v, v)
                    weights_ij = weights[i, j]
                    NMi_w = Ni * Mj * weights_ij
                    numerator += NMi_w * control_points[i, j]
                    denominator += NMi_w
            surface[iu, iv] = numerator / denominator if denominator != 0 else numerator
    return surface


# we add this specifically for surfaces
def eval_bspline_surface(i, degree, nodes, parameter):
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
                * eval_bspline_surface(i, degree - 1, nodes, parameter)
            )
    if (i + degree + 1) < n:
        denom2 = nodes[i + degree + 1] - nodes[i + 1]
        if denom2 != 0:
            second_part = (
                (nodes[i + degree + 1] - parameter)
                / denom2
                * eval_bspline_surface(i + 1, degree - 1, nodes, parameter)
            )
    return first_part + second_part
