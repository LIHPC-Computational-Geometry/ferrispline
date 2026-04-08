use ndarray::{Array1, Array2, Axis};
use num_integer::binomial;

#[derive(Debug)]
pub struct BezierCurve {
    pub degree: usize,
    pub weights: Array1<f64>,
    pub controle_points: Array2<f64>,
}

impl BezierCurve {
    pub fn new_with_weights(
        degree: usize,
        controle_points: Array2<f64>,
        weights: Array1<f64>,
    ) -> Result<Self, String> {
        if controle_points.nrows() != weights.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                controle_points.nrows(),
                weights.len()
            ));
        }
        if controle_points.nrows() != degree + 1 {
            return Err(format!(
                "Degree count mismatch: {} control points vs {} degree",
                controle_points.nrows(),
                degree
            ));
        }
        Ok(Self {
            degree,
            weights,
            controle_points,
        })
    }

    pub fn new(degree: usize, controle_points: Array2<f64>) -> Result<Self, String> {
        if controle_points.nrows() != degree + 1 {
            return Err(format!(
                "Degree count mismatch: {} control points vs {} degree",
                controle_points.nrows(),
                degree
            ));
        }
        Ok(Self {
            degree,
            weights: Array1::from(vec![1.0; controle_points.nrows()]),
            controle_points,
        })
    }

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

        let mut points: Array2<f64> = Array2::zeros((3, sample));

        for i in 0..(self.degree + 1) {
            let forces: Array1<f64> = self.bernstein(i, &t);

            let cp_x = self.controle_points[[i, 0]];
            let cp_y = self.controle_points[[i, 1]];
            let cp_z = self.controle_points[[i, 2]];

            let mut row_x = points.row_mut(0);
            row_x += &(&forces * cp_x);

            let mut row_y = points.row_mut(1);
            row_y += &(&forces * cp_y);

            let mut row_z = points.row_mut(2);
            row_z += &(&forces * cp_z);
        }

        points
    }

    /// Evaluate Rationnal Bezier curve with weightsfor a number of points equal to `sample`
    pub fn evaluate_rational(&self, sample: usize) -> Result<Array2<f64>, String> {
        let curve: Array2<f64> = self.rational_basis(sample)?;
        curve.t();
        curve.dot(&self.controle_points);

        todo!()
    }

    fn rational_basis(&self, sample: usize) -> Result<Array2<f64>, String> {
        let t: Array1<f64> = Array1::linspace(0.0, 1.0, sample);
        let mut weighted_strength: Array2<f64> = Array2::zeros((self.degree + 1, sample));

        for i in 0..(self.degree + 1) {
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
    use ndarray::array;

    use super::*;

    #[test]
    fn test_bernstein_extremes() {
        let degree = 3;

        // At t=0, only the first polynomial (v=0) is 1.0 (using 1x0 matrix for dummy)
        let curve_v0: BezierCurve = BezierCurve::new(degree, Array2::zeros((1, 0))).unwrap();
        assert_eq!(curve_v0.bernstein(0, &array![0.0]), array![1.0]);
        let curve_v1: BezierCurve = BezierCurve::new(degree, Array2::zeros((1, 0))).unwrap();
        assert_eq!(curve_v1.bernstein(1, &array![0.0]), array![0.0]);

        // At t=1, only the last polynomial (v=degree) is 1.0
        let curve: BezierCurve = BezierCurve::new(degree, Array2::zeros((1, 0))).unwrap();
        assert_eq!(curve.bernstein(degree, &array![1.0]), array![1.0]);
        let curve: BezierCurve = BezierCurve::new(degree, Array2::zeros((1, 0))).unwrap();
        assert_eq!(curve.bernstein(0, &array![1.0]), array![0.0]);
    }

    #[test]
    fn test_bernstein_partition_of_unity() {
        let degree = 3;

        let mut t_vec = Vec::new();
        for i in 0..10 {
            t_vec.push(i as f64 / 9.0);
        }
        let t_vals = Array1::from(t_vec);
        let mut totals = vec![0.0; 10];

        for v in 0..=degree {
            let curve: BezierCurve = BezierCurve::new(degree, Array2::zeros((0, 0))).unwrap();
            let results = curve.bernstein(v, &t_vals);
            for (i, &res) in results.iter().enumerate() {
                totals[i] += res;
            }
        }

        for total in totals {
            assert!((total - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_bezier_evaluate_simple() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(degree, control_points).unwrap();

        // Evaluate with 3 samples: t=0.0, t=0.5, t=1.0
        let points = curve.evaluate(3);

        assert_eq!(points.nrows(), 3);

        // At t=0.0, point should be p0
        assert!((points[[0, 0]] - 0.0).abs() < 1e-6);
        assert!((points[[1, 0]] - 0.0).abs() < 1e-6);
        assert!((points[[2, 0]] - 0.0).abs() < 1e-6);

        // At t=0.5, point should be 0.25*p0 + 0.5*p1 + 0.25*p2 = (1.0, 1.0, 0.0)
        assert!((points[[0, 1]] - 1.0).abs() < 1e-6);
        assert!((points[[1, 1]] - 1.0).abs() < 1e-6);
        assert!((points[[2, 1]] - 0.0).abs() < 1e-6);

        // At t=1.0, point should be p2
        assert!((points[[0, 2]] - 2.0).abs() < 1e-6);
        assert!((points[[1, 2]] - 0.0).abs() < 1e-6);
        assert!((points[[2, 2]] - 0.0).abs() < 1e-6);
    }
}
