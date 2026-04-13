use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::core::knot::KnotVector;

#[pyclass]
pub struct PyKnotVector {
    _inner: KnotVector,
}

#[pymethods]
impl PyKnotVector {
    #[new]
    pub fn new(knots: Vec<f64>) -> PyResult<Self> {
        let _inner = KnotVector::new(knots).map_err(PyValueError::new_err)?;
        Ok(Self { _inner })
    }
}
