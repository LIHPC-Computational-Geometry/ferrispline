use ndarray::Array1;

use crate::geometry::spline::SplineCurve;

impl SplineCurve {
    /// Moves a specific control point of a curve to a new position.
    pub fn move_control_point(&mut self, _index: usize, _new_pos: Array1<f64>) {
        todo!("move_control_point")
    }

    /// Sets the weight of a specific control point of a curve.
    pub fn set_control_point_weight(&mut self, _index: usize, _weight: f64) {
        todo!("set_control_point_weight")
    }
}

#[cfg(test)]
mod tests {

    // ==========================================
    // 1. move_control_point Tests
    // ==========================================

    #[test]
    fn test_move_control_point_success() {
        todo!("test_move_control_point_success")
    }

    #[test]
    fn test_move_control_point_invalid_index() {
        todo!("test_move_point_invalid_index")
    }

    // ==========================================
    // 2. set_control_point_weight Tests
    // ==========================================

    #[test]
    fn test_set_control_point_weight_success() {
        todo!("test_set_control_point_weight_success")
    }

    #[test]
    fn test_set_control_point_weight_invalid_index() {
        todo!("test_set_control_point_weight_invalid_index")
    }
}
