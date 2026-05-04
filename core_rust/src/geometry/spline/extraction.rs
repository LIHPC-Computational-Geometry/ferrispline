use crate::geometry::bezier::BezierCurve;
use crate::geometry::spline::SplineCurve;
use ndarray::{Array1, Array2, Axis, concatenate, s};

impl SplineCurve {
    fn compute_knot_insertion_matrix(&self, segment_index: usize) -> Result<Array2<f64>, String> {
        let num_knots = self.knots.as_slice().len();

        if segment_index >= num_knots.saturating_sub(1) {
            return Err(format!(
                "segment_index ({}) is out of bound for knots of length {}",
                segment_index, num_knots
            ));
        }

        let mut extraction_matrix = Array2::<f64>::eye(1);

        for degree_step in 1..=self.degree {
            let start_idx = segment_index.saturating_sub(degree_step);
            let end_idx = (segment_index + degree_step + 1).min(num_knots);

            let local_knots = &self.knots.as_slice()[start_idx..=end_idx];

            let mut tmp_matrix_a = Array2::<f64>::zeros((degree_step, degree_step + 1));
            let mut tmp_matrix_b = Array2::<f64>::zeros((degree_step, degree_step + 1));

            for row in 0..degree_step {
                let knot_start = row + 1;
                let knot_end = knot_start + degree_step;

                // NOTE: May be we want return an error and not 0
                let distance = if knot_end < local_knots.len() && knot_start < local_knots.len() {
                    local_knots[knot_end] - local_knots[knot_start]
                } else {
                    0.0
                };

                let (alpha, beta) = if distance != 0.0 {
                    let alpha_val = (local_knots[degree_step] - local_knots[knot_start]) / distance;
                    let beta_val =
                        (local_knots[degree_step + 1] - local_knots[knot_start]) / distance;
                    (alpha_val, beta_val)
                } else {
                    (0.0, 0.0)
                };

                tmp_matrix_a[[row, row]] = 1.0 - alpha;
                tmp_matrix_a[[row, row + 1]] = alpha;

                tmp_matrix_b[[row, row]] = 1.0 - beta;
                tmp_matrix_b[[row, row + 1]] = beta;
            }

            let upper_half = extraction_matrix.dot(&tmp_matrix_a);

            let last_row_idx = extraction_matrix.nrows() - 1;
            let last_row_matrix = extraction_matrix.slice(s![last_row_idx.., ..]);

            let lower_half = last_row_matrix.dot(&tmp_matrix_b);

            extraction_matrix = concatenate(Axis(0), &[upper_half.view(), lower_half.view()])
                .map_err(|e| format!("Failed to concatenate matrices: {}", e))?;
        }

        Ok(extraction_matrix)
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

    /// Converte a NURBS curve to a Bezier curve
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
