import numpy as np
import numpy.typing as npt
from typing import Annotated

# array of floating numbers
FloatArray = npt.NDArray[np.float64]

# --- 1D (Vector) ---
Vector3 = Annotated[FloatArray, "3"]      # Fixed vector (eg: a spatial coordinate x, y, z)
Vector  = Annotated[FloatArray, "N"]      # Dynamic size vector (e.g. nodes, weights, u_vals)

# --- 2D (Matrix) ---
MatrixNx3 = Annotated[FloatArray, "N, 3"] # Matrix with 3 columns (eg: table of points)
MatrixNxN = Annotated[FloatArray, "N, N"] # Square matrix (e.g. insertion matrix)
MatrixMxN = Annotated[FloatArray, "M, N"] # Rectangular matrix (e.g. basis of evaluation)
