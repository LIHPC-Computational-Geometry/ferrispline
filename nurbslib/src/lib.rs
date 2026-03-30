#![allow(unsafe_op_in_unsafe_fn)]

use core_rust::add as rust_add;
use pyo3::prelude::*;

pub mod geometry;

#[pyfunction]
#[pyo3(name = "add_numbers")]
fn add(left: u64, right: u64) -> PyResult<u64> {
    Ok(rust_add(left, right))
}

#[pymodule]
fn nurbslib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(
        geometry::bezier::compute_knot_insertion_matrix,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(add, m)?)?;
    Ok(())
}
