use core_rust::{core::knot::KnotVector, ids::CurveId, model::Model};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{Python, exceptions::PyValueError, prelude::*};

#[pyclass(module = "ferrispline")]
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
    /// Creates a new Bezier curve in the model and returns its unique ID.
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

    #[pyo3(signature = (degree, control_points, knots, weights=None))]
    /// Creates a new NURBS curve in the model and assigns it a unique ID.
    pub fn create_nurbs(
        &mut self,
        degree: usize,
        control_points: PyReadonlyArray2<f64>,
        knots: Vec<f64>,
        weights: Option<PyReadonlyArray1<f64>>,
    ) -> PyResult<String> {
        let ctrl = control_points.as_array().to_owned();
        let k = KnotVector::new(knots).map_err(PyValueError::new_err)?;
        let w = weights.map(|w| w.as_array().to_owned());
        let id = self
            .inner
            .create_nurbs(degree, ctrl, k, w)
            .map_err(PyValueError::new_err)?;
        Ok(id.to_string())
    }

    /// Deletes the curve associated with the given ID.
    pub fn delete_curve(&mut self, curve_id: &str) -> PyResult<bool> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        Ok(self.inner.delete_curve(&id))
    }

    /// Evaluates the points of the curve associated with the given ID.
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

    /// Evaluates a preview of the curve based on provided geometric parameters.
    /// The `curve_data` tuple contains: (control_points, weights, knots).
    pub fn preview_evaluate<'py>(
        &self,
        py: Python<'py>,
        kind: String,
        degree: usize,
        curve_data: (
            PyReadonlyArray2<'py, f64>,
            Option<PyReadonlyArray1<'py, f64>>,
            Option<Vec<f64>>,
        ),
        sample: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let kind = match kind.as_str() {
            "bezier" => core_rust::model::CurveKind::Bezier,
            "nurbs" => core_rust::model::CurveKind::Nurbs,
            _ => return Err(PyValueError::new_err("Invalid curve kind")),
        };

        // Unpack the tuple parameters
        let (cp_py, cp_w_py, knots_vec) = curve_data;

        let cp = cp_py.as_array().to_owned();
        let cp_w = cp_w_py.map(|w| w.as_array().to_owned());

        // Map directly using the function reference and convert Result to Option
        let knots = knots_vec.map(KnotVector::new).and_then(Result::ok);

        let pts = core_rust::model::Model::preview_evaluate(kind, degree, cp, cp_w, knots, sample)
            .map_err(PyValueError::new_err)?;

        Ok(pts.into_pyarray(py))
    }

    /// Returns an array of control points for the specified curve.
    pub fn get_control_points<'py>(
        &self,
        py: Python<'py>,
        curve_id: &str,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let pts = self
            .inner
            .get_control_points(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(pts.into_pyarray(py))
    }

    /// Returns the degree of the curve.
    pub fn get_degree(&self, curve_id: &str) -> PyResult<usize> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        self.inner
            .get_degree(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
    }

    pub fn get_knots<'py>(
        &self,
        py: Python<'py>,
        curve_id: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let knots = self
            .inner
            .get_knots(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(knots.0.into_pyarray(py))
    }

    pub fn get_weights<'py>(
        &self,
        py: Python<'py>,
        curve_id: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let weights = self
            .inner
            .get_weights(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(weights.into_pyarray(py))
    }

    pub fn curve_kind(&self, curve_id: &str) -> PyResult<String> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let kind = self
            .inner
            .curve_kind(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        match kind {
            core_rust::model::CurveKind::Bezier => Ok("bezier".to_string()),
            core_rust::model::CurveKind::Nurbs => Ok("nurbs".to_string()),
        }
    }

    /// Converts a set of curves to a different curve kind (e.g., Bezier to NURBS).
    pub fn convert(&mut self, curve_ids: Vec<String>, new_kind: String) -> PyResult<Vec<String>> {
        let ids = curve_ids
            .into_iter()
            .map(|s| CurveId::try_from_str(&s).map_err(PyValueError::new_err))
            .collect::<Result<Vec<_>, _>>()?;
        let new_kind = match new_kind.as_str() {
            "bezier" => core_rust::model::CurveKind::Bezier,
            "nurbs" => core_rust::model::CurveKind::Nurbs,
            _ => return Err(PyValueError::new_err("Invalid curve kind")),
        };
        let out = self
            .inner
            .convert(&ndarray::Array1::from_vec(ids), new_kind)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(out.into_iter().map(|id| id.to_string()).collect())
    }

    pub fn move_control_point(
        &mut self,
        curve_id: &str,
        index: usize,
        new_pos: PyReadonlyArray1<f64>,
    ) -> PyResult<()> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        let pos = new_pos.as_array().to_owned();
        self.inner
            .move_control_point(&id, index, pos)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(())
    }

    pub fn set_control_point_weight(
        &mut self,
        curve_id: &str,
        index: usize,
        weight: f64,
    ) -> PyResult<()> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        self.inner
            .set_control_point_weight(&id, index, weight)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))?;
        Ok(())
    }

    pub fn is_dirty(&self, curve_id: &str) -> PyResult<bool> {
        let id = CurveId::try_from_str(curve_id).map_err(PyValueError::new_err)?;
        self.inner
            .is_dirty(&id)
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
    }
}
