import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from .core_types import MatrixNx3, VectorN
from .geometry.bezier import bezier_curves
from .geometry.nurbs import eval_nurbs_curve, cox_deboor

import nurbslib


def n_figure(
    degree: int, knots: list, control_points: MatrixNx3, ctrl_pt_weights: VectorN
):
    try:
        curve = nurbslib.PySplineCurve(degree, control_points, ctrl_pt_weights, knots)
        print("Super ! La courbe NURBS a été créée avec succès.")
    except ValueError as e:
        # PyO3 renverra proprement tes erreurs Rust (PyValueError) ici
        print(f"Erreur lors de la création : {e}")

    nurbs_curve = np.array(curve.eval_nurbs_curve(100))

    bezier_segments = curve.to_bezier()

    # DRAW
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")
    colors = mpl.colormaps.get_cmap("tab10")
    colors = colors.resampled(len(bezier_segments))

    evaluated_segments = []
    for idx, segment in enumerate(bezier_segments):
        seg_points = np.array(segment.evaluate(100, rational=True))
        evaluated_segments.append(seg_points)

        ax.plot(
            seg_points[:, 0],
            seg_points[:, 1],
            seg_points[:, 2],
            color=colors(idx),
            label=f"Bézier {idx+1}",
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

    if evaluated_segments:
        bezier_points = [seg[0] for seg in evaluated_segments]
        bezier_points.append(evaluated_segments[-1][-1])
        bezier_points = np.array(bezier_points)

        ax.scatter(
            bezier_points[:, 0],
            bezier_points[:, 1],
            bezier_points[:, 2],
            color="red",
            s=50,
            label="Bezier points",
        )

    ax.set_title("NURBS to Bezier Rust Implementation")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    ax.grid(True)

    return fig
