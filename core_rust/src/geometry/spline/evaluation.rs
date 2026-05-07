use crate::geometry::spline::SplineCurve;
use crate::traits::ParametricCurve;
use nalgebra::Vector3;
use ndarray::{Array1, Array2};

impl SplineCurve {
    pub fn evaluate(&self, sample: usize) -> Result<Array2<f64>, String> {
        let domain = self.domain();
        let u_vals = Array1::linspace(domain.0, domain.1, sample);

        let mut points: Array2<f64> = Array2::zeros((3, sample));
        for (idx, u) in u_vals.iter().enumerate() {
            let mut numerator = Vector3::zeros();
            let mut denominator = 0.0;

            for i in 0..self.control_points.nrows() {
                let n = self.cox_de_boor(i, self.degree, *u)?;
                let weight_n = self.weights[i] * n;

                let cp_row = self.control_points.row(i);
                numerator += Vector3::new(cp_row[0], cp_row[1], cp_row[2]) * weight_n;
                denominator += weight_n;
            }

            let point = if denominator.abs() < 1e-9 {
                Vector3::zeros()
            } else {
                numerator / denominator
            };

            points[[0, idx]] = point.x;
            points[[1, idx]] = point.y;
            points[[2, idx]] = point.z;
        }
        Ok(points)
    }

    /// Implementation of the Cox-De Boor algorithm
    pub(crate) fn cox_de_boor(&self, i: usize, degree: usize, u: f64) -> Result<f64, String> {
        let n = self.knots.as_slice().len() - 1;

        if i >= n {
            return Err(format!(
                "Index i {} is out of bounds for knot vector of length {}",
                i,
                self.knots.as_slice().len()
            ));
        }
        if degree == 0 {
            if i < n
                && u >= self.knots.as_slice()[i]
                && (u < self.knots.as_slice()[i + 1]
                    || (u <= self.knots.as_slice()[i + 1] && u == self.knots.as_slice()[n - 1]))
            {
                return Ok(1.0);
            } else {
                return Ok(0.0);
            }
        }

        let mut first_part = 0.0;
        let mut second_part = 0.0;

        if i + self.degree < n {
            let denom1 = self.knots.as_slice()[i + self.degree] - self.knots.as_slice()[i];
            if denom1 != 0.0 {
                first_part =
                    (u - self.knots.as_slice()[i]) / denom1 * self.cox_de_boor(i, degree - 1, u)?;
            }
            if i + degree + 1 < n {
                let denom2 = self.knots.as_slice()[i + degree + 1] - self.knots.as_slice()[i + 1];
                if denom2 != 0.0 {
                    second_part = (self.knots.as_slice()[i + degree + 1] - u) / denom2
                        * self.cox_de_boor(i + 1, degree - 1, u)?;
                }
            }
        }
        Ok(first_part + second_part)
    }
}

#[cfg(test)]
mod tests {
    use crate::{core::knot::KnotVector, geometry::spline::SplineCurve};
    use ndarray::{Array2, array};

    // ==========================================
    // 2. Core Algorithm Tests (Cox-De Boor)
    // ==========================================

    #[test]
    /// Tests the basis step function (degree 0) of the Cox-De Boor algorithm.
    fn test_cox_de_boor_degree_0() {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((2, 3));
        let spline = SplineCurve::new(1, ctrl_pts, knots).unwrap();

        // At degree 0, the basis function should return 1.0 if u is within the knot span [0.0, 1.0), else 0.0
        let val_inside = spline.cox_de_boor(1, 0, 0.5).unwrap();
        assert_eq!(val_inside, 1.0, "Should be 1.0 inside the span");

        let val_outside = spline.cox_de_boor(0, 0, 0.5).unwrap();
        assert_eq!(val_outside, 0.0, "Should be 0.0 outside the span");
    }

    #[test]
    /// Tests the Partition of Unity: the sum of all basis functions for any u must equal 1.0.
    fn test_partition_of_unity() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((3, 3)); // degree 2, 3 points
        let spline = SplineCurve::new(2, ctrl_pts, knots).unwrap();

        let u_val = 0.5;
        let mut sum = 0.0;
        for i in 0..spline.control_points.nrows() {
            sum += spline.cox_de_boor(i, spline.degree, u_val).unwrap();
        }

        assert!(
            (sum - 1.0).abs() < 1e-9,
            "Sum of basis functions should be exactly 1.0, got {}",
            sum
        );
    }

    // ==========================================
    // 3. Evaluation Tests
    // ==========================================

    #[test]
    /// Tests that a clamped curve starts and ends exactly on its first and last control points.
    fn test_evaluate_clamped_properties() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(); // Clamped
        let ctrl_pts = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let spline = SplineCurve::new(2, ctrl_pts, knots).unwrap();

        let samples = 10;
        let evaluated_points = spline.evaluate(samples).unwrap();

        // Check start point (first column) matches [1.0, 2.0, 3.0]
        assert!((evaluated_points[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((evaluated_points[[1, 0]] - 2.0).abs() < 1e-6);
        assert!((evaluated_points[[2, 0]] - 3.0).abs() < 1e-6);

        // Check end point (last column) matches [7.0, 8.0, 9.0]
        let last_idx = samples - 1;
        assert!((evaluated_points[[0, last_idx]] - 7.0).abs() < 1e-6);
        assert!((evaluated_points[[1, last_idx]] - 8.0).abs() < 1e-6);
        assert!((evaluated_points[[2, last_idx]] - 9.0).abs() < 1e-6);
    }

    #[test]
    /// Tests that evaluation returns a matrix with 3 rows (X, Y, Z) and the correct number of samples.
    fn test_evaluate_dimensions() {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((2, 3));
        let spline = SplineCurve::new(1, ctrl_pts, knots).unwrap();

        let samples = 100;
        let points = spline.evaluate(samples).unwrap();

        assert_eq!(
            points.nrows(),
            3,
            "Matrix should always have 3 rows for X, Y, Z"
        );
        assert_eq!(
            points.ncols(),
            samples,
            "Matrix should have columns equal to requested sample count"
        );
    }
}
