import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

# Ajout du dossier 'src' au sys.path pour permettre les imports absolus
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from nurbs_math.visualisation import figure
from nurbs_math.nurbslib_visualisation import n_figure
from nurbs_math.load_nurbs import default_value, load_nurbs_from_vtk


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

    n_figure(
        degree=degree,
        knots=knots,
        control_points=control_points,
        ctrl_pt_weights=ctrl_pt_weights,
    )

    figure(
        degree=degree,
        knots=knots,
        control_points=control_points,
        ctrl_pt_weights=ctrl_pt_weights,
    )

    plt.show()


if __name__ == "__main__":
    main()
