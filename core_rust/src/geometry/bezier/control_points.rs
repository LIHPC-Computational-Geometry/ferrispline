use super::BezierCurve;
use ndarray::Array1;

impl BezierCurve {
    pub fn move_control_point(
        &mut self,
        _index: usize,
        _new_pos: Array1<f64>,
    ) -> Result<(), String> {
        todo!("move_control_point")
    }

    pub fn set_control_point_weight(&mut self, _index: usize, _weight: f64) -> Result<(), String> {
        todo!("set_control_point_weight")
    }
}
