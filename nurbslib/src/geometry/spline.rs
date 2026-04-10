use ndarray::{Array1, Array2};
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
        points_py: Vec<[f64; 3]>,
        weight_py: Vec<f64>,
        knots: Vec<f64>,
    ) -> PyResult<Self> {
        let mut control_points = Array2::<f64>::zeros((points_py.len(), 3));
        for (i, p) in points_py.iter().enumerate() {
            control_points[[i, 0]] = p[0];
            control_points[[i, 1]] = p[1];
            control_points[[i, 2]] = p[2];
        }

        let weights = Array1::from(weight_py);

        let knot_vector = KnotVector::new(knots).map_err(PyValueError::new_err)?;

        let inner = SplineCurve::builder()
            .degree(degree)
            .build_nurbs(control_points, weights, knot_vector)
            .map_err(PyValueError::new_err)?;

        Ok(Self { inner })
    }

    pub fn eval_nurbs_curve(&self, sample: usize) -> PyResult<Vec<[f64; 3]>> {
        let curve_points = self
            .inner
            .eval_nurbs_curve(sample)
            .map_err(PyValueError::new_err)?;

        let cols = curve_points.ncols();
        let mut py_points = Vec::with_capacity(cols);

        for i in 0..cols {
            py_points.push([
                curve_points[[0, i]],
                curve_points[[1, i]],
                curve_points[[2, i]],
            ]);
        }

        Ok(py_points)
    }

    pub fn to_bezier(&self) -> PyResult<Vec<PyBezierCurve>> {
        let beziers = self.inner.to_bezier().map_err(PyValueError::new_err)?;
        Ok(beziers
            .into_iter()
            .map(|b| PyBezierCurve { inner: b })
            .collect())
    }
}
