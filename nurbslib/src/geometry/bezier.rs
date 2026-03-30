use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::geometry::bezier::compute_knot_insertion_matrix as rs_compute_knot_insertion_matrix;

#[pyfunction]
#[pyo3(name = "compute_knot_insertion_matrix")]
pub fn compute_knot_insertion_matrix(
    knots: Vec<f64>,
    degree: usize,
    segment_index: usize,
) -> PyResult<Vec<Vec<f64>>> {
    let matrix = rs_compute_knot_insertion_matrix(&knots, degree, segment_index)
        .map_err(PyValueError::new_err)?;

    let mut py_matrix = Vec::with_capacity(matrix.nrows());
    for r in 0..matrix.nrows() {
        let mut row = Vec::with_capacity(matrix.ncols());
        for c in 0..matrix.ncols() {
            row.push(matrix[(r, c)]);
        }
        py_matrix.push(row);
    }
    Ok(py_matrix)
}
