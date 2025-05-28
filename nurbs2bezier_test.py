import numpy as np
from matplotlib.testing.compare import compare_images
from nurbs2bezier import figure


def test_segment():

    nodes = [0, 0, 0, 0, 1 / 5, 2 / 5, 2 / 5, 3 / 5, 3 / 5, 3 / 5, 1, 1, 1, 1]

    degree = 3

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

    weights = np.array([1, 2, 2, 1, 0.5, 0.5, 1, 1, 2, 1])

    figure(degree=degree, nodes=nodes, control_points=control_points, weights=weights)

    assert compare_images("segment_ref.png", "segment.png", 0.001) == None
