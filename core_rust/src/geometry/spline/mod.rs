use crate::core::knot::KnotVector;
use crate::traits::ParametricCurve;
use ndarray::{Array1, Array2};

mod control_points;
mod degree;
mod evaluation;
mod extraction;
mod knots;

//
// Main SplineCurve Struct
//
#[derive(Debug, Clone)]
pub struct SplineCurve {
    pub degree: usize,
    pub control_points: Array2<f64>,
    pub weights: Array1<f64>,
    pub knots: KnotVector,
}

//
// Main impl block for SplineCurve
//
impl SplineCurve {
    /// Creates a new NURBS curve with the given degree, control points, weights, and knots.
    /// Returns an error if the degree is zero, if the knot count is invalid, or if the number of weights does not match the number of control points.
    pub fn new_with_weights(
        degree: usize,
        control_points: Array2<f64>,
        weights: Array1<f64>,
        knots: KnotVector,
    ) -> Result<Self, String> {
        if degree == 0 {
            return Err("Degree must be at least 1.".to_string());
        }
        knots.lenght_check(control_points.nrows(), degree)?;
        if control_points.nrows() != weights.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                control_points.nrows(),
                weights.len()
            ));
        }
        Ok(Self {
            degree,
            weights,
            control_points,
            knots,
        })
    }

    /// Creates a new B-Spline curve with the given degree, control points, and knots.
    /// All weights are initialized to 1.0.
    /// Returns an error if the degree is zero or if the knot count is invalid based on the control points.
    pub fn new(
        degree: usize,
        control_points: Array2<f64>,
        knots: KnotVector,
    ) -> Result<Self, String> {
        if degree == 0 {
            return Err("Degree must be at least 1.".to_string());
        }
        knots.lenght_check(control_points.nrows(), degree)?;
        Ok(Self {
            degree,
            weights: Array1::ones(control_points.nrows()),
            control_points,
            knots,
        })
    }
}

//
// ParametricCurve Trait Implementation
//
impl ParametricCurve for SplineCurve {
    fn domain(&self) -> (f64, f64) {
        let n = self.knots.as_slice().len() - self.degree - 1;
        let p = self.degree;
        (self.knots.as_slice()[p], self.knots.as_slice()[n])
    }
}

#[cfg(test)]
mod tests {
    use crate::{core::knot::KnotVector, geometry::spline::SplineCurve, traits::ParametricCurve};
    use ndarray::{Array1, Array2};

    // ==========================================
    // 1. Construction & Validation Tests
    // ==========================================

    #[test]
    /// Tests that a B-Spline is successfully created and weights are defaulted to 1.0.
    fn test_new_bspline_success() {
        // 4 control points, degree 3. Required knots (m = n + p + 1): 3 + 3 + 1 = 7 (so 8 elements)
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));

        let spline = SplineCurve::new(3, ctrl_pts, knots).expect("Should build successfully");

        assert_eq!(
            spline.weights.len(),
            4,
            "Weights should be generated for each point"
        );
        assert_eq!(spline.weights[0], 1.0, "Default weights should be 1.0");
    }

    #[test]
    /// Tests that constructor fails if the number of weights doesn't match the control points.
    fn test_new_nurbs_weight_mismatch() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));
        let invalid_weights = Array1::ones(3); // 3 weights for 4 points -> Error!

        let result = SplineCurve::new_with_weights(3, ctrl_pts, invalid_weights, knots);

        assert!(result.is_err(), "Should fail due to weight mismatch");
        assert!(result.unwrap_err().contains("Weight count mismatch"));
    }

    #[test]
    /// Tests that a degree of 0 is rejected.
    fn test_new_degree_zero() {
        let knots = KnotVector::new(vec![0.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((1, 3));

        let result = SplineCurve::new(0, ctrl_pts, knots);

        assert!(result.is_err(), "Degree 0 should be rejected");
        assert!(result.unwrap_err().contains("Degree must be at least 1"));
    }

    // ==========================================
    // 5. ParametricCurve Trait Tests
    // ==========================================

    #[test]
    /// Tests that the domain correctly ignores clamped bounds according to the degree.
    fn test_domain_extraction() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));
        let spline = SplineCurve::new(2, ctrl_pts, knots).unwrap();

        let (u_min, u_max) = spline.domain();

        assert_eq!(u_min, 0.0, "Domain start should be 0.0");
        assert_eq!(u_max, 1.0, "Domain end should be 1.0");
    }
}
