import argparse
import os

import matplotlib.pyplot as plt

from .visualisation import figure
from .load_nurbs import default_value, load_nurbs_from_vtk


def set_parser():
    parser = argparse.ArgumentParser(
        prog="nurbs2bezier.py", description="A script for convert NURBS to Bezier Curve"
    )
    parser.add_argument("-f", "--file", type=str, help="vtk file with the good format")
    return parser.parse_args()


def main():
    args = set_parser()

    if args.file:
        if not os.path.exists(args.file):
            raise FileNotFoundError(f"File {args.file} not found.")
        filepath: str = args.file

        print(f"Loading file: {filepath}...")
        control_points, ctrl_pt_weights, knots, degree = load_nurbs_from_vtk(filepath)
        print(f"Extraction success ! Curve of degree {degree} detected.")
    else:
        control_points, ctrl_pt_weights, knots, degree = default_value()

    fig = figure(
        degree=degree,
        knots=knots,
        control_points=control_points,
        ctrl_pt_weights=ctrl_pt_weights,
    )

    plt.show()


if __name__ == "__main__":
    main()
