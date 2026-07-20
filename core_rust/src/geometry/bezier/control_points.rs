use super::BezierCurve;
use ndarray::Array1;

impl BezierCurve {
    pub fn move_control_point(
        &mut self,
        index: usize,
        new_pos: Array1<f64>,
    ) -> Result<(), String> {
        self.control_points.row_mut(index).assign(&new_pos);
        Ok(())
    }

    pub fn set_control_point_weight(&mut self, index: usize, weight: f64) -> Result<(), String> {
        self.weights[index] = weight;
        Ok(())
    }
}
