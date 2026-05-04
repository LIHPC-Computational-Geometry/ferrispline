use crate::geometry::spline::SplineCurve;

impl SplineCurve {
    pub fn set_degree(&mut self, new_degree: usize) {
        if new_degree > self.degree {
            self.degree_elevation(new_degree);
        } else if new_degree < self.degree {
            self.degree_reduction(new_degree);
        } else {
            // Do nothing
        }
    }

    fn degree_elevation(&mut self, _new_degree: usize) {
        todo!("degree_elevation")
    }

    fn degree_reduction(&mut self, _new_degree: usize) {
        todo!("degree_reduction")
    }
}

#[cfg(test)]
mod tests {

    // ==========================================
    // 1. degree_elevation Tests
    // ==========================================

    #[test]
    fn test_degree_elevation_preserves_exact_geometric_shape() {
        todo!("test_degree_elevation1_preserves_exact_geometric_shape")
    }

    #[test]
    fn test_degree_elevation_increases_knot_multiplicity_by_one() {
        todo!("test_degree_elevation_increases_knot_multiplicity_by_one")
    }

    #[test]
    fn test_degree_elevation_adds_expected_control_points() {
        todo!("test_degree_elevation_adds_expected_control_points")
    }

    #[test]
    fn test_degree_elevation_preserves_joint_continuity() {
        todo!("test_degree_elevation_preserves_joint_continuity")
    }

    // ==========================================
    // 2. degree_reduction Tests
    // ==========================================

    #[test]
    fn test_degree_reduction_minimizes_distance_to_original_curve() {
        todo!("test_degree_reduction_minimizes_distance_to_original_curve")
    }

    #[test]
    fn test_degree_reduction_succeeds_if_error_within_tolerance() {
        todo!("test_degree_reduction_succeeds_if_error_within_tolerance")
    }

    #[test]
    fn test_degree_reduction_fail_if_error_exceeds_threshold() {
        todo!("test_degree_reduction_fail_if_error_exceeds_threshold")
    }

    #[test]
    fn test_degree_reduction_decrease_control_points_and_knots() {
        todo!("test_degree_reduction_decrease_control_points_and_knots")
    }

    // ==========================================
    // 3. set_degree Tests
    // ==========================================

    #[test]
    fn test_set_degree_elevation() {
        todo!("test_set_degree_elevation")
    }

    #[test]
    fn test_set_degree_reduction() {
        todo!("test_set_degree_reduction")
    }

    #[test]
    fn test_set_degree_no_change() {
        todo!("test_set_degree_no_change")
    }
}
