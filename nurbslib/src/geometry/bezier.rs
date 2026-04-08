use ndarray::Array2;
use pyo3::prelude::*;

use core_rust::geometry::bezier::BezierCurve;

#[pyclass]
pub struct PyBezierCurve {
    pub inner: BezierCurve,
}

#[pymethods]
impl PyBezierCurve {
    // NOTE: Clone data points during the conversion maybe will be opti
    #[new]
    pub fn new(degree: usize, points: Vec<[f64; 3]>) -> PyResult<Self> {
        let mut controle_points = Array2::<f64>::zeros((points.len(), 3));
        for (i, p) in points.iter().enumerate() {
            controle_points[[i, 0]] = p[0];
            controle_points[[i, 1]] = p[1];
            controle_points[[i, 2]] = p[2];
        }

        let inner = BezierCurve::new(degree, controle_points);
        Ok(Self { inner })
    }

    pub fn evaluate(&self, sample: usize) -> PyResult<Vec<[f64; 3]>> {
        let curve_points = self.inner.evaluate(sample);

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
}
