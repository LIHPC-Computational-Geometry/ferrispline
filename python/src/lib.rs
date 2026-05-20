#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;

use crate::geometry::{bezier::PyBezierCurve, spline::PySplineCurve};
use crate::model::PyModel;

pub mod core;
pub mod geometry;
pub mod model;

#[pymodule]
fn ferrispline(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBezierCurve>()?;
    m.add_class::<PySplineCurve>()?;
    m.add_class::<PyModel>()?;
    Ok(())
}
