import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from .core_types import MatrixNx3, VectorN

import ferrispline


def n_figure(
    degree: int, knots: list, control_points: MatrixNx3, ctrl_pt_weights: VectorN
):
    model = ferrispline.PyModel()
    try:
        # Use the model to create a NURBS curve, which returns a curve ID.
        curve_id = model.create_nurbs(degree, control_points, knots, ctrl_pt_weights)
        print("Super ! La courbe NURBS a été créée avec succès.")
    except Exception as e:
        # Catch exceptions from the Rust binding
        print(f"Erreur lors de la création : {e}")
        return

    # Evaluate the NURBS curve using its ID.
    nurbs_curve = np.array(model.evaluate(curve_id, 100))

    # Convert the NURBS curve to Bézier segments.
    # NOTE: This assumes the `PyModel` wrapper exposes a `to_bezier` method
    # that returns a list of new curve IDs, one for each segment.
    bezier_segment_ids = model.convert([curve_id], "bezier")

    # DRAW
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")
    colors = mpl.colormaps.get_cmap("tab10")
    colors = colors.resampled(len(bezier_segment_ids))

    evaluated_segments = []
    # Evaluate each Bézier segment by its ID.
    for idx, segment_id in enumerate(bezier_segment_ids):
        seg_points = np.array(model.evaluate(segment_id, 100))
        evaluated_segments.append(seg_points)

        ax.plot(
            seg_points[:, 0],
            seg_points[:, 1],
            seg_points[:, 2],
            color=colors(idx),
            label=f"Bézier {idx + 1}",
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
