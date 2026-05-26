use super::BezierCurve;
use ndarray::Array2;
use num_integer::binomial;
use std::cmp::min;

impl BezierCurve {
    /// Calculates and returns the degree elevation matrix for a Bezier curve.
    ///
    /// # Arguments
    /// * `new_degree` - The desired degree after elevation.
    ///
    /// # Returns
    /// An `Array2<f64>` matrix containing the elevation coefficients.
    fn degree_elevation_matrix(current_degree: usize, new_degree: usize) -> Array2<f64> {
        let rows = new_degree + 1;
        let cols = current_degree + 1;
        let mut mat = Array2::zeros((rows, cols));

        let t = new_degree - current_degree;

        for i in 0..rows {
            let start = i.saturating_sub(t);
            let end = min(current_degree, i);

            for j in start..=end {
                let num = (binomial(current_degree, j) * binomial(t, i - j)) as f64;
                let den = binomial(new_degree, i) as f64;

                mat[[i, j]] = num / den;
            }
        }

        mat
    }

    /// Elevates the degree of the Bezier curve to the specified new degree.
    pub fn degree_elevation(&mut self, new_degree: usize) -> Result<(), String> {
        if new_degree <= self.degree {
            return Err(format!(
                "degree_elevation: new degree must be greater than current degree. {} <= {}",
                new_degree, self.degree
            ));
        }

        let elevation_mat: Array2<f64> =
            BezierCurve::degree_elevation_matrix(self.degree, new_degree);
        self.control_points = elevation_mat.dot(&self.control_points);
        self.weights = elevation_mat.dot(&self.weights);
        self.degree = new_degree;

        Ok(())
    }

    /// reduce the degree of a Bezier curve with the inverse formulas
    /// Pli = (Qi - i/n * Pi-1) / (1 - i/n)
    /// Pri = (Qi+1 - (1 - (i+1)/n) * Pi+1) / ((i+1)/n)
    // NOTE: This function will be not implemented yet because not used in the bot's project
    fn _degree_reduction_by_one(&mut self, _tolerance: f64) -> Result<Vec<BezierCurve>, String> {
        todo!("degree_reduction_by_one: not implemented")
    }

    /// Reduce degree locally by solving `E * Q ≈ P` in least squares sense, where
    /// `E` is the degree elevation matrix from `target_degree` to `self.degree`.
    ///
    /// For rational curves, the system is solved on homogeneous coordinates.
    pub fn _reduce_local_least_squares(
        &self,
        _target_degree: usize,
    ) -> Result<BezierCurve, String> {
        todo!("reduce_local_least_squares: not implemented")
    }

    pub fn set_degree(&mut self, degree: usize, _tolerance: Option<f64>) -> Result<(), String> {
        if degree > self.degree {
            self.degree_elevation(degree)
        } else if degree < self.degree {
            todo!("reduce_local_least_squares: not implemented")
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array2, array};

    // ==========================================
    // SECTION 5. Degree Elevation Tests
    // ==========================================

    #[test]
    fn test_degree_elevation_point_count() {
        let degree = 2;
        let control_points = array![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let mut curve = BezierCurve::new(degree, control_points).unwrap();

        // Élévation de n (2) à n+1 (3)
        let _ = curve.degree_elevation(3);

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
        let _ = curve.degree_elevation(3);

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

        let _ = curve.degree_elevation(3);
        let elevated_eval = curve.evaluate(samples);

        // Comparaison coordonnée par coordonnée (X, Y, Z) pour chaque échantillon
        for axis in 0..3 {
            for t_step in 0..samples {
                let p_orig = original_eval[[t_step, axis]];
                let p_elev = elevated_eval[[t_step, axis]];
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
        let _ = curve.degree_elevation(target_degree);

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

        for t_step in 0..samples {
            for axis in 0..3 {
                let p_orig = original_eval[[t_step, axis]];
                let p_elev = elevated_eval[[t_step, axis]];
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
