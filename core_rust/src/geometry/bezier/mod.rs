use ndarray::{Array1, Array2};

mod control_points;
mod conversion;
mod degree;
mod evaluation;
mod subdivision;

#[derive(Debug, Clone)]
pub struct BezierCurve {
    pub degree: usize,
    pub control_points: Array2<f64>,
    pub weights: Array1<f64>,
}

impl BezierCurve {
    pub fn new_with_weights(
        degree: usize,
        control_points: Array2<f64>,
        weights: Array1<f64>,
    ) -> Result<Self, String> {
        if control_points.nrows() != weights.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                control_points.nrows(),
                weights.len()
            ));
        }
        if control_points.nrows() != degree + 1 {
            return Err(format!(
                "Degree count mismatch: {} control points vs {} degree",
                control_points.nrows(),
                degree
            ));
        }
        Ok(Self {
            degree,
            weights,
            control_points,
        })
    }

    pub fn new(degree: usize, control_points: Array2<f64>) -> Result<Self, String> {
        if control_points.nrows() != degree + 1 {
            return Err(format!(
                "Degree count mismatch: {} control points vs {} degree",
                control_points.nrows(),
                degree
            ));
        }
        Ok(Self {
            degree,
            weights: Array1::from(vec![1.0; control_points.nrows()]),
            control_points,
        })
    }
}

// ==========================================
// SECTION Unit Tests
// ==========================================

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, Array2};

    // ==========================================
    // SECTION 1. Constructor Validation Tests
    // ==========================================

    #[test]
    fn test_new_degree_mismatch() {
        let degree = 2;
        // Degree 2 needs 3 control points, we provide 4
        let control_points = Array2::zeros((4, 3));

        let result = BezierCurve::new(degree, control_points);
        assert!(
            result.is_err(),
            "Should fail because 4 points != degree 2 + 1"
        );
        assert!(result.unwrap_err().contains("Degree count mismatch"));
    }

    #[test]
    fn test_new_with_weights_mismatch() {
        let degree = 2;
        let control_points = Array2::zeros((3, 3));
        let weights = Array1::zeros(2);

        let result = BezierCurve::new_with_weights(degree, control_points, weights);
        assert!(
            result.is_err(),
            "Should fail because point count != weight count"
        );
        assert!(result.unwrap_err().contains("Weight count mismatch"));
    }
}
