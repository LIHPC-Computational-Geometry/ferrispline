use super::BezierCurve;
use ndarray::{Array1, Array2, Axis};
use num_integer::binomial;

impl BezierCurve {
    pub fn bernstein(&self, v: usize, t: &Array1<f64>) -> Array1<f64> {
        t.iter()
            .map(|t| {
                binomial(self.degree, v) as f64
                    * t.powf(v as f64)
                    * (1.0 - t).powf((self.degree - v) as f64)
            })
            .collect()
    }

    /// Evaluate Bezier curve for a number of points equal to `sample`
    pub fn evaluate(&self, sample: usize) -> Array2<f64> {
        let t: Array1<f64> = Array1::linspace(0.0, 1.0, sample);

        let mut points: Array2<f64> = Array2::zeros((sample, 3));

        for i in 0..=self.degree {
            let forces: Array1<f64> = self.bernstein(i, &t);

            for dir in 0..3 {
                let cp = self.control_points[[i, dir]];
                let mut col = points.column_mut(dir);
                col += &(&forces * cp);
            }
        }

        points
    }

    /// Evaluate Rationnal Bezier curve with weights for a number of points equal to `sample`
    pub fn evaluate_rational(&self, sample: usize) -> Result<Array2<f64>, String> {
        let basis: Array2<f64> = self.rational_basis(sample)?;
        let t_basis = basis.t();

        let curve_points = t_basis.dot(&self.control_points);
        Ok(curve_points)
    }

    fn rational_basis(&self, sample: usize) -> Result<Array2<f64>, String> {
        let t: Array1<f64> = Array1::linspace(0.0, 1.0, sample);
        let mut weighted_strength: Array2<f64> = Array2::zeros((self.degree + 1, sample));

        for i in 0..=self.degree {
            let forces: Array1<f64> = Array1::from(self.bernstein(i, &t));
            weighted_strength
                .row_mut(i)
                .assign(&(forces * self.weights[i]));
        }
        let denominator: Array1<f64> = weighted_strength.sum_axis(Axis(0));

        if denominator.iter().any(|&x| x.abs() < 1e-9) {
            return Err(
                "Division by zero in rational basis: denominator contains zero values".to_string(),
            );
        }
        Ok(weighted_strength / denominator)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, Array2, array};

    // ==========================================
    // #SECTION 2. Bernstein Basis Tests
    // ==========================================

    #[test]
    fn test_bernstein_extremes() {
        let degree = 3;
        let dummy_points = Array2::zeros((4, 3));
        let curve = BezierCurve::new(degree, dummy_points).unwrap();

        // At t=0, only the first polynomial (v=0) is 1.0
        assert_eq!(curve.bernstein(0, &array![0.0]), array![1.0]);
        assert_eq!(curve.bernstein(1, &array![0.0]), array![0.0]);

        // At t=1, only the last polynomial (v=degree) is 1.0
        assert_eq!(curve.bernstein(degree, &array![1.0]), array![1.0]);
        assert_eq!(curve.bernstein(0, &array![1.0]), array![0.0]);
    }

    #[test]
    fn test_bernstein_partition_of_unity() {
        let degree = 3;
        let dummy_points = Array2::zeros((4, 3));
        let curve = BezierCurve::new(degree, dummy_points).unwrap();

        let t_vals = Array1::linspace(0.0, 1.0, 10);
        let mut totals = vec![0.0; 10];

        for v in 0..=degree {
            let results = curve.bernstein(v, &t_vals);
            for (i, &res) in results.iter().enumerate() {
                totals[i] += res;
            }
        }

        // Sum of all basis polynomials for any t must be exactly 1.0
        for total in totals {
            assert!((total - 1.0).abs() < 1e-10);
        }
    }

    // ==========================================
    // #SECTION 3. Evaluation Tests
    // ==========================================

    #[test]
    fn test_bezier_evaluate_simple() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(degree, control_points).unwrap();

        // Evaluate with 3 samples: t=0.0, t=0.5, t=1.0
        let points = curve.evaluate(3);

        assert_eq!(
            points.nrows(),
            3,
            "Matrix should have rows equal to requested sample count"
        );
        assert_eq!(
            points.ncols(),
            3,
            "Matrix should always have 3 columns for X, Y, Z"
        );

        // At t=0.0, point should be p0
        assert!((points[[0, 0]] - 0.0).abs() < 1e-6);
        assert!((points[[0, 1]] - 0.0).abs() < 1e-6);

        // At t=0.5, point should be 0.25*p0 + 0.5*p1 + 0.25*p2 = (1.0, 1.0, 0.0)
        assert!((points[[1, 0]] - 1.0).abs() < 1e-6);
        assert!((points[[1, 1]] - 1.0).abs() < 1e-6);

        // At t=1.0, point should be p2
        assert!((points[[2, 0]] - 2.0).abs() < 1e-6);
        assert!((points[[2, 1]] - 0.0).abs() < 1e-6);
    }

    // ==========================================
    // SECTION 4. Rational Basis Tests
    // ==========================================

    #[test]
    fn test_rational_basis_division_by_zero() {
        let degree = 2;
        let control_points = Array2::zeros((3, 3));
        // All weights set to 0.0 to trigger division by zero
        let weights = Array1::zeros(3);

        let curve = BezierCurve::new_with_weights(degree, control_points, weights).unwrap();

        let result = curve.rational_basis(10);
        assert!(
            result.is_err(),
            "Should return an error when denominator is zero"
        );
        assert!(result.unwrap_err().contains("Division by zero"));
    }
}
