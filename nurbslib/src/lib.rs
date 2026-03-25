use pyo3::prelude::*;

use core_rust::add as rust_add;

#[pyfunction]
#[pyo3(name = "add_numbers")]
fn add(left: u64, right: u64) -> PyResult<u64> {
    Ok(rust_add(left, right))
}

#[pymodule]
fn api_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    Ok(())
}