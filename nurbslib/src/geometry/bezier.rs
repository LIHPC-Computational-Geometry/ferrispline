use ndarray::{Array1, Array2};
use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::geometry::bezier::BezierCurve;

#[pyclass]
pub struct PyBezierCurve {
    pub inner: BezierCurve,
}

#[pymethods]
impl PyBezierCurve {
    // NOTE: Clone data points during the conversion maybe will be opti
    #[new]
    pub fn new(degree: usize, points: Vec<[f64; 3]>, weights: Option<Vec<f64>>) -> PyResult<Self> {
        let mut control_points = Array2::<f64>::zeros((points.len(), 3));
        for (i, p) in points.iter().enumerate() {
            control_points[[i, 0]] = p[0];
            control_points[[i, 1]] = p[1];
            control_points[[i, 2]] = p[2];
        }

        let inner = match weights {
            Some(w) => {
                let weights_array = Array1::from(w);
                BezierCurve::new_with_weights(degree, control_points, weights_array)
                    .map_err(PyValueError::new_err)?
            }
            None => BezierCurve::new(degree, control_points).map_err(PyValueError::new_err)?,
        };
        Ok(Self { inner })
    }

    pub fn degree(&self) -> usize {
        self.inner.degree
    }

    // On définit la signature Python : sample est obligatoire, rational est optionnel (None par défaut)
    #[pyo3(signature = (sample, rational=None))]
    pub fn evaluate(&self, sample: usize, rational: Option<bool>) -> PyResult<Vec<[f64; 3]>> {
        let use_rational =
            rational.unwrap_or_else(|| self.inner.weights.iter().any(|&w| (w - 1.0).abs() > 1e-9));

        let curve_points = if use_rational {
            self.inner
                .evaluate_rational(sample)
                .map_err(PyValueError::new_err)?
        } else {
            self.inner.evaluate(sample)
        };

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
