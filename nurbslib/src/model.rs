use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::{ids::CurveId, model::Model};

#[pyclass]
pub struct PyModel {
    inner: Model,
}

impl Default for PyModel {
    fn default() -> Self {
        Self::new()
    }
}

#[pymethods]
impl PyModel {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: Model::new(),
        }
    }

    #[pyo3(signature = (degree, control_points, weights=None))]
    pub fn create_bezier(
        &mut self,
        degree: usize,
        control_points: PyReadonlyArray2<f64>,
        weights: Option<PyReadonlyArray1<f64>>,
    ) -> PyResult<String> {
        let ctrl = control_points.as_array().to_owned();
        let w = weights.map(|w| w.as_array().to_owned());
        let id = self
            .inner
            .create_bezier(degree, ctrl, w)
            .map_err(PyValueError::new_err)?;
        Ok(id.to_string())
    }

    pub fn delete_curve(&mut self, curve_id: &str) -> PyResult<bool> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        Ok(self.inner.delete_curve(&id))
    }

    pub fn evaluate<'py>(
        &self,
        py: Python<'py>,
        curve_id: &str,
        sample: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let pts = self
            .inner
            .evaluate(&id, sample)
            .map_err(PyValueError::new_err)?;
        Ok(pts.into_pyarray(py))
    }
}
