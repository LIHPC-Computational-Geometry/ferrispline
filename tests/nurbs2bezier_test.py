import numpy as np
from matplotlib.testing.compare import compare_images
from nurbs_math.visualisation import figure
from nurbs_math.load_nurbs import default_value


def test_segment():

    control_points, ctrl_pt_weights, knots, degree = default_value()

    fig = figure(
        degree=degree,
        knots=knots,
        control_points=control_points,
        ctrl_pt_weights=ctrl_pt_weights,
    )

    fig.savefig("segment.png")

    assert compare_images("./segment_ref.png", "segment.png", 0.0) == None
