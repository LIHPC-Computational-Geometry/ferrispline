import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
import matplotlib.cm as cm


def computeConversionMatrix(nodes, degree, interval_index):
    S = np.eye(1)
    for k in range(1, degree + 1):
        begin = max(0, interval_index - k)
        end = min(len(nodes), interval_index + k + 2)
        local_nodes = nodes[begin:end]

        A = np.zeros((k, k + 1))
        B = np.zeros((k, k + 1))

        for l in range(k):
            idx = l + 1
            denominator = local_nodes[idx + k] - local_nodes[idx]
            alpha = (
                (local_nodes[k] - local_nodes[idx]) / denominator
                if denominator != 0
                else 0.0
            )
            beta = (
                (local_nodes[k + 1] - local_nodes[idx]) / denominator
                if denominator != 0
                else 0.0
            )

            A[l, l] = 1 - alpha
            A[l, l + 1] = alpha
            B[l, l] = 1 - beta
            B[l, l + 1] = beta

        S_up = S @ A
        S_down = S[-1:] @ B
        S = np.vstack([S_up, S_down])

    return S


def bernstein(i, degree, t):
    return comb(degree, i) * (t**i) * (1 - t) ** (degree - i)


def evalBezierCurve(control_points, weights, degree, sample=100):
    r"""
    Evaluates rational Bezier curve and returns it.

    The formula is:

        # Math: \text{curve}(t) = \frac{1}{\sum_{i=0}^{\text{degree}} \text{weights}[i] \begin{pmatrix} \text{degree} \\ i \end{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)}} \sum_{i=0}^{\text{degree}} \text{weights}[i] \begin{pmatrix} \text{degree} \\ i \end{pmatrix} t^{i} (1 - t)^{(\text{degree} - i)} \text{control_points}[i]

    :param control_points: control points vector.
    :param weights: weights vector.
    :param degree: Bezier basis degree.
    :param sample: render sample.
    """

    t = np.linspace(0, 1, sample)
    nominator = np.zeros((sample, control_points.shape[1]))
    denominator = np.zeros(sample)

    for i in range(degree + 1):
        B = bernstein(i, degree, t)
        WB = weights[i] * B
        denominator += WB
        nominator += WB[:, np.newaxis] * control_points[i]

    return nominator / denominator[:, np.newaxis]


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
    u_min = nodes[degree_u]
    u_max = nodes[-degree_u - 1]
    u_vals = np.linspace(u_min, u_max, nb_points_u)

    v_min = nodes[degree_v]
    v_max = nodes[-degree_v - 1]
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


def figure(degree, nodes, control_points, weights):
    bezier_segments = []

    for i in range(degree, len(nodes) - degree - 1):
        if nodes[i] != nodes[i + 1]:
            first = i - degree
            last = i
            if first < 0 or last >= len(control_points):
                continue

            S = computeConversionMatrix(nodes, degree, i)
            print(f"Matrix S for interval {i} :")
            print(S)
            print()

            local_points = control_points[first : last + 1]
            local_weights = weights[first : last + 1]

            weighted_points = local_weights[:, np.newaxis] * local_points
            bezier_weighted_points = S @ weighted_points
            bezier_weights = S @ local_weights

            bezier_points = bezier_weighted_points / bezier_weights[:, np.newaxis]

            curve = evalBezierCurve(bezier_points, bezier_weights, degree)
            bezier_segments.append(curve)

    nurbs_curve = evalNURBSCurve(nodes, control_points, weights, degree)

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

    # u_vals = np.linspace(nodes[degree], nodes[-degree - 1], 1000)
    x_vals = []
    y_vals = []
    z_vals = []
    for u in nodes[degree:-degree]:
        numerator = np.zeros(control_points.shape[1])
        denominator = 0.0
        for i in range(len(control_points)):
            N = evalBspline(i, degree, nodes, u)
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
    plt.savefig("segment.png")
