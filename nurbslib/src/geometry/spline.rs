use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use core_rust::core::knot::KnotVector;
use core_rust::geometry::spline::SplineCurve;

use crate::geometry::bezier::PyBezierCurve;

#[pyclass]
pub struct PySplineCurve {
    pub inner: SplineCurve,
}

#[pymethods]
impl PySplineCurve {
    #[new]
    pub fn new(
        degree: usize,
        points_py: PyReadonlyArray2<f64>,
        weight_py: PyReadonlyArray1<f64>,
        knots: Vec<f64>,
    ) -> PyResult<Self> {
        let control_points = points_py.as_array().to_owned();
        let weights = weight_py.as_array().to_owned();
        let knot_vector = KnotVector::new(knots).map_err(PyValueError::new_err)?;

        let inner = SplineCurve::builder()
            .degree(degree)
            .build_nurbs(control_points, weights, knot_vector)
            .map_err(PyValueError::new_err)?;

        Ok(Self { inner })
    }

    pub fn eval_nurbs_curve<'py>(
        &self,
        py: Python<'py>,
        sample: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let curve_points = self
            .inner
            .eval_nurbs_curve(sample)
            .map_err(PyValueError::new_err)?;

        Ok(curve_points.into_pyarray(py))
    }

    pub fn to_bezier(&self) -> PyResult<Vec<PyBezierCurve>> {
        let beziers = self.inner.to_bezier().map_err(PyValueError::new_err)?;
        Ok(beziers
            .into_iter()
            .map(|b| PyBezierCurve { inner: b })
            .collect())
    }
}
