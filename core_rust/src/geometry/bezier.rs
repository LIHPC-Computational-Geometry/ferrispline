use nalgebra::min;
use ndarray::{Array1, Array2, Axis};
use num_integer::binomial;

#[derive(Debug)]
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

        for i in 0..=self.degree {
            let forces: Array1<f64> = self.bernstein(i, &t);

            for dir in 0..3 {
                let cp = self.control_points[[i, dir]];
                let mut row = points.row_mut(dir);
                row += &(&forces * cp);
            }
        }

        points
    }

    /// Evaluate Rationnal Bezier curve with weights for a number of points equal to `sample`
    pub fn evaluate_rational(&self, sample: usize) -> Result<Array2<f64>, String> {
        let basis: Array2<f64> = self.rational_basis(sample)?;
        let t_basis = basis.t();
        let curve_points = t_basis.dot(&self.control_points);
        Ok(curve_points.t().to_owned())
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

    /// Increasing the degree of a Bezier curve without changing its shape
    pub fn degree_elevation(&mut self, new_degree: usize) {
        if new_degree <= self.degree {
            return;
        }

        let t = new_degree - self.degree;

        let mut new_control_points: Array2<f64> = Array2::zeros((new_degree + 1, 3));
        let mut new_weights: Array1<f64> = Array1::zeros(new_degree + 1);

        for i in 0..=new_degree {
            let start = i.saturating_sub(t);
            let end = min(self.degree, i);

            for j in start..=end {
                let num = (binomial(self.degree, j) * binomial(t, i - j)) as f64;
                let den = binomial(new_degree, i) as f64;
                let coef = num / den;

                new_control_points[[i, 0]] += self.control_points[[j, 0]] * coef;
                new_control_points[[i, 1]] += self.control_points[[j, 1]] * coef;
                new_control_points[[i, 2]] += self.control_points[[j, 2]] * coef;

                new_weights[i] += self.weights[j] * coef;
            }
        }
        self.degree = new_degree;
        self.control_points = new_control_points;
        self.weights = new_weights;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, Array2, array};

    // ==========================================
    // 1. Constructor Validation Tests
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

    // ==========================================
    // 2. Bernstein Basis Tests
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
    // 3. Evaluation Tests
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
            "Matrix should always have 3 rows for X, Y, Z"
        );
        assert_eq!(
            points.ncols(),
            3,
            "Matrix should have columns equal to requested sample count"
        );

        // At t=0.0, point should be p0
        assert!((points[[0, 0]] - 0.0).abs() < 1e-6);
        assert!((points[[1, 0]] - 0.0).abs() < 1e-6);

        // At t=0.5, point should be 0.25*p0 + 0.5*p1 + 0.25*p2 = (1.0, 1.0, 0.0)
        assert!((points[[0, 1]] - 1.0).abs() < 1e-6);
        assert!((points[[1, 1]] - 1.0).abs() < 1e-6);

        // At t=1.0, point should be p2
        assert!((points[[0, 2]] - 2.0).abs() < 1e-6);
        assert!((points[[1, 2]] - 0.0).abs() < 1e-6);
    }

    // ==========================================
    // 4. Rational Basis Tests
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

    // ==========================================
    // 5. Degree Elevation Tests
    // ==========================================

    #[test]
    fn test_degree_elevation_point_count() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let mut curve = BezierCurve::new(degree, control_points).unwrap();

        // Élévation de n (2) à n+1 (3)
        curve.degree_elevation(3);

        assert_eq!(
            curve.degree, 3,
            "Le degré de la courbe doit être mis à jour à 3"
        );
        assert_eq!(
            curve.control_points.nrows(),
            4,
            "Passer au degré 3 nécessite exactement 4 points de contrôle (n + 1)"
        );
    }

    #[test]
    fn test_degree_elevation_endpoints() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let mut curve = BezierCurve::new(degree, control_points.clone()).unwrap();
        curve.degree_elevation(3);

        // Vérification du premier point (P_0 == Q_0)
        assert!((curve.control_points[[0, 0]] - control_points[[0, 0]]).abs() < 1e-9);
        assert!((curve.control_points[[0, 1]] - control_points[[0, 1]]).abs() < 1e-9);
        assert!((curve.control_points[[0, 2]] - control_points[[0, 2]]).abs() < 1e-9);

        // Vérification du dernier point (P_n == Q_{n+1})
        let last_idx_orig = 2;
        let last_idx_elev = 3;
        assert!(
            (curve.control_points[[last_idx_elev, 0]] - control_points[[last_idx_orig, 0]]).abs()
                < 1e-9
        );
        assert!(
            (curve.control_points[[last_idx_elev, 1]] - control_points[[last_idx_orig, 1]]).abs()
                < 1e-9
        );
        assert!(
            (curve.control_points[[last_idx_elev, 2]] - control_points[[last_idx_orig, 2]]).abs()
                < 1e-9
        );
    }

    #[test]
    fn test_degree_elevation_shape_invariance() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let mut curve = BezierCurve::new(degree, control_points).unwrap();

        let samples = 50; // Nombre de points à évaluer pour vérifier la forme

        // Utilisation de la méthode evaluate existante qui retourne une matrice (3, samples)
        let original_eval: Array2<f64> = curve.evaluate(samples);

        curve.degree_elevation(3);
        let elevated_eval = curve.evaluate(samples);

        // Comparaison coordonnée par coordonnée (X, Y, Z) pour chaque échantillon
        for axis in 0..3 {
            for t_step in 0..samples {
                let p_orig = original_eval[[axis, t_step]];
                let p_elev = elevated_eval[[axis, t_step]];
                let diff = (p_orig - p_elev).abs();

                assert!(
                    diff < 1e-9,
                    "Invariance de forme non respectée sur l'axe {} à l'échantillon {}. Différence: {}",
                    axis,
                    t_step,
                    diff
                );
            }
        }
    }

    #[test]
    fn test_degree_elevation_multiple_steps() {
        let degree = 3;
        // Une courbe cubique (degré 3) nécessite 4 points de contrôle
        // On utilise des coordonnées 3D distinctes pour bien tester chaque axe
        let control_points = array![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 1.0],
            [2.0, -1.0, 2.0],
            [3.0, 0.0, 0.0]
        ];
        let samples = 100; // Un échantillonnage plus fin pour une courbe de plus haut degré
        let mut curve = BezierCurve::new(degree, control_points.clone()).unwrap();
        let original_eval = curve.evaluate(samples);

        // Élévation directe de 3 à 6 (saut de 3 degrés)
        let target_degree = 6;
        curve.degree_elevation(target_degree);

        // 1. Vérification de la structure dimensionnelle
        assert_eq!(
            curve.degree, target_degree,
            "Le degré final doit être exactement de 6"
        );
        assert_eq!(
            curve.control_points.nrows(),
            target_degree + 1,
            "Une courbe de degré 6 doit posséder exactement 7 points de contrôle"
        );

        // 2. Vérification de l'ancrage strict des extrémités
        // L'extrémité de départ (t=0) doit rester identique
        assert!((curve.control_points[[0, 0]] - control_points[[0, 0]]).abs() < 1e-9);
        assert!((curve.control_points[[0, 1]] - control_points[[0, 1]]).abs() < 1e-9);
        assert!((curve.control_points[[0, 2]] - control_points[[0, 2]]).abs() < 1e-9);

        // L'extrémité d'arrivée (t=1) doit rester identique
        let last_idx_orig = degree; // Indice 3
        let last_idx_elev = target_degree; // Indice 6
        assert!(
            (curve.control_points[[last_idx_elev, 0]] - control_points[[last_idx_orig, 0]]).abs()
                < 1e-9
        );
        assert!(
            (curve.control_points[[last_idx_elev, 1]] - control_points[[last_idx_orig, 1]]).abs()
                < 1e-9
        );
        assert!(
            (curve.control_points[[last_idx_elev, 2]] - control_points[[last_idx_orig, 2]]).abs()
                < 1e-9
        );

        // 3. Validation de l'invariance totale de la forme
        let elevated_eval = curve.evaluate(samples);

        for axis in 0..3 {
            for t_step in 0..samples {
                let p_orig = original_eval[[axis, t_step]];
                let p_elev = elevated_eval[[axis, t_step]];
                let diff = (p_orig - p_elev).abs();

                assert!(
                    diff < 1e-9,
                    "Divergence géométrique détectée suite au saut de degré. Axe {}, Echantillon {}. Différence: {}",
                    axis,
                    t_step,
                    diff
                );
            }
        }
    }
}
