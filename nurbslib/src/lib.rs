#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;

pub mod geometry;

#[pymodule]
fn nurbslib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(
        geometry::bezier::compute_knot_insertion_matrix,
        m
    )?)?;
    Ok(())
}
