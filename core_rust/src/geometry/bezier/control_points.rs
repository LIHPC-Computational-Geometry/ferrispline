use super::BezierCurve;
use ndarray::Array1;

impl BezierCurve {
    pub fn move_control_point(
        &mut self,
        _index: usize,
        _new_pos: Array1<f64>,
    ) -> Result<(), String> {
        self.control_points.row_mut(_index).assign(&_new_pos);
        Ok(())
    }

    pub fn set_control_point_weight(&mut self, _index: usize, _weight: f64) -> Result<(), String> {
        self.weights[_index] = _weight;
        Ok(())
    }
}
