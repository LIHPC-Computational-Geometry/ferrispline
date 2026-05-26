use crate::geometry::bezier::BezierCurve;
use crate::geometry::spline::SplineCurve;
use ndarray::{Array1, Array2, Axis, s};

impl SplineCurve {
    fn compute_knot_insertion_matrix(&self, segment_index: usize) -> Result<Array2<f64>, String> {
        let p = self.degree;
        let num_knots = self.knots.as_slice().len();
        let i = segment_index;

        if i >= num_knots.saturating_sub(1) {
            return Err(format!(
                "segment_index ({}) is out of bound for knots of length {}",
                i, num_knots
            ));
        }

        // Vecteur de nœuds local affectant le segment
        let start_idx = i.saturating_sub(p);
        let end_idx = (i + p + 1).min(num_knots - 1);
        let mut current_knots = self.knots.as_slice()[start_idx..=end_idx].to_vec();

        // Matrice identité représentant les p+1 points de contrôle locaux initiaux
        let mut c = Array2::<f64>::eye(p + 1);

        let u_left = self.knots.as_slice()[i];
        let u_right = self.knots.as_slice()[i + 1];

        // Évaluation des multiplicités existantes
        let mut mult_left = 0;
        let mut mult_right = 0;
        for &k_val in &current_knots {
            if (k_val - u_left).abs() < 1e-9 {
                mult_left += 1;
            }
            if (k_val - u_right).abs() < 1e-9 {
                mult_right += 1;
            }
        }

        // Algorithme de Boehm pour l'insertion d'un nœud `t`
        let mut insert_knot = |t: f64| {
            let mut k = 0;
            for (idx, &val) in current_knots.iter().enumerate() {
                if val <= t + 1e-9 {
                    k = idx;
                } else {
                    break;
                }
            }

            let num_pts = c.nrows();
            let mut new_c = Array2::<f64>::zeros((num_pts + 1, p + 1));

            for j in 0..=num_pts {
                let alpha = if j <= k.saturating_sub(p) {
                    1.0
                } else if j > k {
                    0.0
                } else {
                    let num = t - current_knots[j];
                    let den = current_knots[j + p] - current_knots[j];
                    if den == 0.0 { 0.0 } else { num / den }
                };

                for col in 0..=p {
                    let p_prev = if j == 0 { 0.0 } else { c[[j - 1, col]] };
                    let p_curr = if j == num_pts { 0.0 } else { c[[j, col]] };
                    new_c[[j, col]] = alpha * p_curr + (1.0 - alpha) * p_prev;
                }
            }
            current_knots.insert(k + 1, t);
            c = new_c;
        };

        // Insertion des extrémités jusqu'à multiplicité `p`
        for _ in mult_left..p {
            insert_knot(u_left);
        }
        for _ in mult_right..p {
            insert_knot(u_right);
        }

        // Le polygone de Bézier correspond aux p+1 points bornés par les nœuds u_right
        let mut first_right_idx = 0;
        for (idx, &k_val) in current_knots.iter().enumerate() {
            if (k_val - u_right).abs() < 1e-9 {
                first_right_idx = idx;
                break;
            }
        }

        let start_row = first_right_idx.saturating_sub(p + 1);
        let mut final_c = Array2::<f64>::zeros((p + 1, p + 1));

        for r in 0..=p {
            for col in 0..=p {
                final_c[[r, col]] = c[[start_row + r, col]];
            }
        }

        Ok(final_c)
    }

    fn new_controle_points(
        &self,
        idx: usize,
        start_idx: usize,
        end_idx: usize,
    ) -> Result<(Array2<f64>, Array1<f64>), String> {
        let knot_insertion_matrix = self.compute_knot_insertion_matrix(idx)?;

        let local_ctrl_pt = &self.control_points.slice(s![start_idx..end_idx, ..]);
        let local_weights = &self.weights.slice(s![start_idx..end_idx]);

        let weighted_points = &local_weights.insert_axis(Axis(1)) * local_ctrl_pt;

        let extracted_weighted_points = knot_insertion_matrix.dot(&weighted_points);
        let extracted_weights: Array1<f64> = knot_insertion_matrix.dot(local_weights);

        let final_points: Array2<f64> =
            &extracted_weighted_points / &extracted_weights.view().insert_axis(Axis(1));

        Ok((final_points, extracted_weights))
    }

    /// Converte a NURBS curve to Bezier curves
    pub fn to_bezier(&self) -> Result<Vec<BezierCurve>, String> {
        let mut bezier_curves: Vec<BezierCurve> = Vec::new();
        for i in self.degree..self.control_points.nrows() {
            if self.knots.as_slice()[i] == self.knots.as_slice()[i + 1] {
                continue;
            }
            let ctrl_pt_start_idx = i - self.degree;
            let ctrl_pt_end_idx = i + 1;

            if ctrl_pt_end_idx >= self.control_points.len() {
                continue;
            }

            let (control_points, weights) =
                self.new_controle_points(i, ctrl_pt_start_idx, ctrl_pt_end_idx)?;

            let bezier_curve = BezierCurve::new_with_weights(self.degree, control_points, weights)?;
            bezier_curves.push(bezier_curve);
        }
        Ok(bezier_curves)
    }
}

#[cfg(test)]
mod tests {
    use crate::{core::knot::KnotVector, geometry::spline::SplineCurve};
    use ndarray::Array2;

    // ==========================================
    // 4. Bezier Extraction & Knot Insertion Tests
    // ==========================================

    #[test]
    /// Tests that the extraction matrix has the correct (degree + 1) x (degree + 1) shape.
    fn test_compute_knot_insertion_matrix_dimension() {
        // Degree 2. 5 points (n=4). Knots = n+p+1 = 4+2+1=7 -> 8 knots
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0]).unwrap();
        let degree = 2;
        let segment_index = 3;
        let ctrl_pts = Array2::zeros((5, 3));
        let spline = SplineCurve::new(degree, ctrl_pts, knots).unwrap();

        let matrix = spline
            .compute_knot_insertion_matrix(segment_index)
            .expect("The function returned an error instead of the matrix");

        assert_eq!(
            matrix.nrows(),
            degree + 1,
            "The number of rows is incorrect"
        );
        assert_eq!(
            matrix.ncols(),
            degree + 1,
            "The number of columns is incorrect"
        );
    }

    #[test]
    /// Test that an error is returned if the segment_index is out of bounds
    fn test_compute_knot_insertion_matrix_invalid_segment_index() {
        // 2 points, degree 1 -> m = 1+1+1 = 3 (4 knots)
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let degree = 1;
        let invalid_segment_index = 3;
        let ctrl_pts = Array2::zeros((2, 3));
        let spline = SplineCurve::new(degree, ctrl_pts, knots).unwrap();

        let result = spline.compute_knot_insertion_matrix(invalid_segment_index);

        assert!(
            result.is_err(),
            "The function should have returned an error (Err)!"
        );
        assert!(result.unwrap_err().contains("segment_index"));
    }

    #[test]
    /// Tests that a spline is correctly broken down into the proper amount of Bezier segments.
    fn test_to_bezier_segment_count() {
        // Spline with 2 valid segments: [0.0, 0.5] and [0.5, 1.0]
        // 4 points, degree 2. 7 knots total
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));
        let spline = SplineCurve::new(2, ctrl_pts, knots).unwrap();

        let bezier_curves = spline.to_bezier().unwrap();
        assert_eq!(
            bezier_curves.len(),
            2,
            "There should be exactly 2 distinct Bezier segments generated"
        );
    }
}
