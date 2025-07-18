import numpy as np
from matplotlib.testing.compare import compare_images
from nurbs2bezier import article


def test_segment():

    fig = article()

    fig.savefig("segment.png")

    assert compare_images("segment_ref.png", "segment.png", 0.0) == None
