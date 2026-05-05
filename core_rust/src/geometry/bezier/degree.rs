use super::BezierCurve;
use nalgebra::min;
use ndarray::{Array1, Array2};
use num_integer::binomial;

impl BezierCurve {
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

                for k in 0..3 {
                    new_control_points[[i, k]] += self.control_points[[j, k]] * coef;
                }

                new_weights[i] += self.weights[j] * coef;
            }
        }
        self.degree = new_degree;
        self.control_points = new_control_points;
        self.weights = new_weights;
    }

    fn degree_reduction(&mut self, _new_degree: usize) {
        todo!("create a function degree_reduction")
    }

    pub fn set_degree(&mut self, degree: usize) -> Result<(), String> {
        if degree > self.degree {
            self.degree_elevation(degree);
        } else if degree < self.degree {
            self.degree_reduction(degree);
        }
        Ok(())
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
