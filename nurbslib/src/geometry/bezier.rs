use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::geometry::bezier::BezierCurve;

#[pyclass]
pub struct PyBezierCurve {
    pub inner: BezierCurve,
}

#[pymethods]
impl PyBezierCurve {
    #[new]
    pub fn new(
        degree: usize,
        points: PyReadonlyArray2<f64>,
        weights: Option<PyReadonlyArray1<f64>>,
    ) -> PyResult<Self> {
        // NOTE: 'as_array()' crée une vue (aucune copie)
        // NOTE: 'to_owned() fait une copie en bloc très rapide (optimiser en C)
        let control_points = points.as_array().to_owned();

        let inner = match weights {
            Some(w) => {
                let weights_array = w.as_array().to_owned();
                BezierCurve::new_with_weights(degree, control_points, weights_array)
                    .map_err(PyValueError::new_err)?
            }
            None => BezierCurve::new(degree, control_points).map_err(PyValueError::new_err)?,
        };
        Ok(Self { inner })
    }

    pub fn get_degree(&self) -> usize {
        self.inner.degree
    }

    // NOTE: struct BezierCurve doesn't impl copy
    pub fn get_control_points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.control_points.clone().into_pyarray(py)
    }

    pub fn get_weights<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.weights.clone().into_pyarray(py)
    }

    // On définit la signature Python : sample est obligatoire, rational est optionnel (None par défaut)
    #[pyo3(signature = (sample, rational=None))]
    pub fn evaluate<'py>(
        &self,
        py: Python<'py>,
        sample: usize,
        rational: Option<bool>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let use_rational =
            rational.unwrap_or_else(|| self.inner.weights.iter().any(|&w| (w - 1.0).abs() > 1e-9));

        let curve_points = if use_rational {
            self.inner
                .evaluate_rational(sample)
                .map_err(PyValueError::new_err)?
        } else {
            self.inner.evaluate(sample)
        };

        Ok(curve_points.into_pyarray(py))
    }
}
