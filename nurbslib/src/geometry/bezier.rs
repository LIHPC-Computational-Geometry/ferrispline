use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::{exceptions::PyValueError, prelude::*};

use core_rust::geometry::bezier::BezierCurve;

#[pyclass]
pub struct PyBezierCurve {
    pub inner: BezierCurve,
}

//
#[pymethods]
impl PyBezierCurve {
    #[new]
    pub fn new(
        degree: usize,
        points: PyReadonlyArray2<f64>,
        weights: Option<PyReadonlyArray1<f64>>,
    ) -> PyResult<Self> {
        // NOTE: `as_array()` creates a view (no copy is made).
        // NOTE: `to_owned()` performs a very fast block copy (optimized in C).
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

    // NOTE: `into_pyarray()` is a conversion function that consumes the variable.
    // NOTE: `to_pyarray()` creates a new NumPy array by copying the data from a Rust reference.
    // NOTE: We cannot use `into_pyarray()` here because `self` is borrowed, and dynamically allocated types like `Array1` and `Array2` do not implement the `Copy` trait.
    pub fn get_control_points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.inner.control_points.to_pyarray(py)
    }

    pub fn get_weights<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.weights.to_pyarray(py)
    }

    // NOTE: On définit la signature Python : sample est obligatoire, rational est optionnel (None par défaut)
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

    pub fn degree_elevation(&mut self, new_degree: usize) {
        self.inner.degree_elevation(new_degree);
    }
}
