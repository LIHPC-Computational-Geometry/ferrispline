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

        for degree_step in 1..(self.degree + 1) {
            let start_idx = segment_index.saturating_sub(degree_step);
            let end_idx = (segment_index + degree_step + 2).min(num_knots);

            let local_knots = &self.knots.as_slice()[start_idx..end_idx];

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

    // NOTE: Check que la fonction to_bezier fonction
    // NOTE: Faire des tests unitaire
    // NOTE: faire un cargo check et fmt
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
        for i in self.degree..(self.knots.as_slice().len() - self.degree - 1) {
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

            let bezier_curve = BezierCurve::new_with_weights(self.degree, controle_points, weights);
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
    use ndarray::Array2;

    use crate::{core::knot::KnotVector, geometry::spline::SplineCurve};

    #[test]
    /// Test that the extraction matrix has the correct (degree + 1) x (degree + 1) shape
    fn test_compute_knot_insertion_matrix_dimension() {
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

        let error_msg = result.unwrap_err();
        assert!(
            error_msg.contains("segment_index"),
            "The error message does not contain 'segment_index'"
        );
        assert!(
            error_msg.contains("is out of bound for knots of length"),
            "The error message does not contain the expected text"
        );
    }
}
