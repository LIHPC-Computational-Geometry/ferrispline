use crate::traits::ParametricCurve;
use crate::{core::knot::KnotVector, geometry::bezier::BezierCurve};
use nalgebra::Vector3;
use ndarray::{Array1, Array2, Axis, concatenate, s};

//
// Main SplineCurve Struct
//
#[derive(Debug)]
pub struct SplineCurve {
    pub degree: usize,
    pub controle_points: Array2<f64>,
    pub weights: Array1<f64>,
    pub knots: KnotVector,
}

//
// SplineCurveBuilder and its implementation
//
pub struct SplineCurveBuilder {
    degree: usize,
}

impl SplineCurveBuilder {
    pub fn new(degree: usize) -> Self {
        Self { degree }
    }

    fn default() -> Self {
        Self { degree: 0 }
    }

    pub fn degree(mut self, degree: usize) -> Self {
        self.degree = degree;
        self
    }

    pub fn build_bspline(
        self,
        controle_points: Array2<f64>,
        knots: KnotVector,
    ) -> Result<SplineCurve, String> {
        self.validate(controle_points.nrows(), &knots)?;
        let weights = Array1::ones(controle_points.nrows());
        Ok(SplineCurve {
            degree: self.degree,
            controle_points,
            weights,
            knots,
        })
    }

    pub fn build_nurbs(
        self,
        controle_points: Array2<f64>,
        weights: Array1<f64>,
        knots: KnotVector,
    ) -> Result<SplineCurve, String> {
        self.validate(controle_points.nrows(), &knots)?;
        if controle_points.nrows() != weights.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                controle_points.len(),
                weights.len()
            ));
        }
        Ok(SplineCurve {
            degree: self.degree,
            controle_points,
            weights,
            knots,
        })
    }

    fn validate(&self, num_ctrl_points: usize, knots: &KnotVector) -> Result<(), String> {
        if self.degree == 0 {
            return Err("Degree must be at least 1.".to_string());
        }
        knots.lenght_check(num_ctrl_points, self.degree)
    }
}

//
// Main impl block for SplineCurve
//
impl SplineCurve {
    pub fn builder() -> SplineCurveBuilder {
        SplineCurveBuilder::default()
    }

    pub fn eval_nurbs_curve(&self, sample: usize) -> Result<Array2<f64>, String> {
        let domain = self.domain();
        let u_vals = Array1::linspace(domain.0, domain.1, sample);

        let mut points: Array2<f64> = Array2::zeros((3, sample));
        for (idx, u) in u_vals.iter().enumerate() {
            let mut numerator = Vector3::zeros();
            let mut denominator = 0.0;

            for i in 0..self.controle_points.nrows() {
                let n = self.cox_de_boor(i, self.degree, *u)?;
                let weight_n = self.weights[i] * n;

                let cp_row = self.controle_points.row(i);
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
    fn cox_de_boor(&self, i: usize, degree: usize, u: f64) -> Result<f64, String> {
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

        let local_ctrl_pt = &self.controle_points.slice(s![start_idx..end_idx, ..]);
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
        for i in self.degree..self.controle_points.nrows() {
            if self.knots.as_slice()[i] == self.knots.as_slice()[i + 1] {
                continue;
            }
            let ctrl_pt_start_idx = i - self.degree;
            let ctrl_pt_end_idx = i + 1;

            if ctrl_pt_end_idx >= self.controle_points.len() {
                continue;
            }

            let (controle_points, weights) =
                self.new_controle_points(i, ctrl_pt_start_idx, ctrl_pt_end_idx)?;

            let bezier_curve =
                BezierCurve::new_with_weights(self.degree, controle_points, weights)?;
            bezier_curves.push(bezier_curve);
        }
        Ok(bezier_curves)
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
    use ndarray::{Array1, Array2, array};

    // ==========================================
    // 1. Construction & Validation Tests
    // ==========================================

    #[test]
    /// Tests that a B-Spline is successfully created and weights are defaulted to 1.0.
    fn test_builder_bspline_success() {
        // 4 control points, degree 3. Required knots (m = n + p + 1): 3 + 3 + 1 = 7 (so 8 elements)
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));

        let spline = SplineCurve::builder()
            .degree(3)
            .build_bspline(ctrl_pts, knots)
            .expect("Should build successfully");

        assert_eq!(
            spline.weights.len(),
            4,
            "Weights should be generated for each point"
        );
        assert_eq!(spline.weights[0], 1.0, "Default weights should be 1.0");
    }

    #[test]
    /// Tests that the builder fails if the number of weights doesn't match the control points.
    fn test_builder_nurbs_weight_mismatch() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));
        let invalid_weights = Array1::ones(3); // 3 weights for 4 points -> Error!

        let result = SplineCurve::builder()
            .degree(3)
            .build_nurbs(ctrl_pts, invalid_weights, knots);

        assert!(result.is_err(), "Should fail due to weight mismatch");
        assert!(result.unwrap_err().contains("Weight count mismatch"));
    }

    #[test]
    /// Tests that a degree of 0 is rejected by the builder.
    fn test_builder_degree_zero() {
        let knots = KnotVector::new(vec![0.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((1, 3));

        let result = SplineCurve::builder()
            .degree(0)
            .build_bspline(ctrl_pts, knots);

        assert!(result.is_err(), "Degree 0 should be rejected");
        assert!(result.unwrap_err().contains("Degree must be at least 1"));
    }

    // ==========================================
    // 2. Core Algorithm Tests (Cox-De Boor)
    // ==========================================

    #[test]
    /// Tests the basis step function (degree 0) of the Cox-De Boor algorithm.
    fn test_cox_de_boor_degree_0() {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((2, 3));
        let spline = SplineCurve::builder()
            .degree(1)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

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
        let spline = SplineCurve::builder()
            .degree(2)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

        let u_val = 0.5;
        let mut sum = 0.0;
        for i in 0..spline.controle_points.nrows() {
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
    fn test_eval_nurbs_curve_clamped_properties() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(); // Clamped
        let ctrl_pts = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let spline = SplineCurve::builder()
            .degree(2)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

        let samples = 10;
        let evaluated_points = spline.eval_nurbs_curve(samples).unwrap();

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
    fn test_eval_nurbs_curve_dimensions() {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((2, 3));
        let spline = SplineCurve::builder()
            .degree(1)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

        let samples = 100;
        let points = spline.eval_nurbs_curve(samples).unwrap();

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
        let spline = SplineCurve::builder()
            .degree(degree)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

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
        let spline = SplineCurve::builder()
            .degree(degree)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

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
        let spline = SplineCurve::builder()
            .degree(2)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

        let bezier_curves = spline.to_bezier().unwrap();
        assert_eq!(
            bezier_curves.len(),
            2,
            "There should be exactly 2 distinct Bezier segments generated"
        );
    }

    // ==========================================
    // 5. ParametricCurve Trait Tests
    // ==========================================

    #[test]
    /// Tests that the domain correctly ignores clamped bounds according to the degree.
    fn test_domain_extraction() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]).unwrap();
        let ctrl_pts = Array2::zeros((4, 3));
        let spline = SplineCurve::builder()
            .degree(2)
            .build_bspline(ctrl_pts, knots)
            .unwrap();

        let (u_min, u_max) = spline.domain();

        // p = 2, n = 3. knots[p] = knots[2] = 0.0. knots[n] = knots[3] = 0.5.
        // Oh wait, knots array is length 7. n = 4 points - 1 = 3.
        // formula: n = knots.len() - p - 1 = 7 - 2 - 1 = 4.
        // Let's rely on the evaluation!
        assert_eq!(u_min, 0.0, "Domain start should be 0.0");
        assert_eq!(u_max, 1.0, "Domain end should be 1.0");
    }
}
