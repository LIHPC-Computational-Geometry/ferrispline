use nalgebra::Point3;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use core_rust::core::knot::KnotVector;
use core_rust::geometry::spline::SplineCurve;

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
        let controle_points: Vec<Point3<f64>> = points_py
            .into_iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect();

        let inner = SplineCurve::builder()
            .degree(degree)
            .build_nurbs(controle_points, weight_py, KnotVector(knots))
            .map_err(PyValueError::new_err)?;
        Ok(Self { inner })
    }

    pub fn eval_nurbs_curve(&self, sample: usize) -> PyResult<Vec<[f64; 3]>> {
        let curve_points = self
            .inner
            .eval_nurbs_curve(sample)
            .map_err(PyValueError::new_err)?;

        Ok(curve_points.iter().map(|p| [p.x, p.y, p.z]).collect())
    }
}
