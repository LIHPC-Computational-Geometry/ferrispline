import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from .core_types import MatrixNx3, VectorN
from .geometry.bezier import bezier_curves
from .geometry.nurbs import eval_nurbs_curve, cox_deboor


def figure(
    degree: int, knots: list, control_points: MatrixNx3, ctrl_pt_weights: VectorN
):
    bezier_segments: list = bezier_curves(
        knots, control_points, ctrl_pt_weights, degree
    )
    nurbs_curve: MatrixNx3 = eval_nurbs_curve(
        knots, control_points, ctrl_pt_weights, degree
    )

    # DRAW
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")
    colors = mpl.colormaps.get_cmap("tab10")
    colors = colors.resampled(len(bezier_segments))

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

    if bezier_segments:
        bezier_points = [segment[0] for segment in bezier_segments]
        bezier_points.append(bezier_segments[-1][-1])
        bezier_points = np.array(bezier_points)

        ax.scatter(
            bezier_points[:, 0],
            bezier_points[:, 1],
            bezier_points[:, 2],
            color="red",
            s=50,
            label="Bezier points",
        )

    ax.set_title("NURBS to Bezier")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    ax.grid(True)

    return fig
