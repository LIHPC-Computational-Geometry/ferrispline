use super::BezierCurve;
use crate::geometry::spline::SplineCurve;

impl BezierCurve {
    pub fn to_nurbs(&self) -> Result<SplineCurve, String> {
        todo!("to_nurbs")
    }
}
