use nalgebra::{DMatrix, Point3};

// use crate::traits::PointTrait;

#[derive(Debug)]
pub struct BezierCurve {
    pub(crate) _degree: usize,
    pub(crate) _controle_points: Vec<Point3<f64>>,
}

impl BezierCurve {
    pub fn new(_degree: usize, _controle_points: Vec<Point3<f64>>) -> Self {
        Self {
            _degree,
            _controle_points,
        }
    }
}

// TODO: create an error type
pub fn compute_knot_insertion_matrix(
    knots: &[f64],
    degree: usize,
    segment_index: usize,
) -> Result<DMatrix<f64>, String> {
    let num_knots = knots.len();

    if segment_index >= num_knots.saturating_sub(1) {
        return Err(format!(
            "segment_index ({}) is out of bound for knots of length {}",
            segment_index, num_knots
        ));
    }

    let mut extraction_matrix = DMatrix::<f64>::identity(1, 1);

    for degree_step in 1..(degree + 1) {
        let start_idx = segment_index.saturating_sub(degree);
        let end_idx = (segment_index + degree_step + 2).min(num_knots);

        let local_knots = &knots[start_idx..end_idx];

        let mut tmp_matrix_a = DMatrix::<f64>::zeros(degree_step, degree_step + 1);
        let mut tmp_matrix_b = DMatrix::<f64>::zeros(degree_step, degree_step + 1);

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
                let beta_val = (local_knots[degree_step + 1] - local_knots[knot_start]) / distance;
                (alpha_val, beta_val)
            } else {
                (0.0, 0.0)
            };

            tmp_matrix_a[(row, row)] = 1.0 - alpha;
            tmp_matrix_a[(row, row + 1)] = alpha;

            tmp_matrix_b[(row, row)] = 1.0 - beta;
            tmp_matrix_b[(row, row + 1)] = beta;
        }

        let upper_half = &extraction_matrix * &tmp_matrix_a;

        let last_row = extraction_matrix.row(extraction_matrix.nrows() - 1);
        let last_row_matrix = DMatrix::from_row_slice(
            1,
            extraction_matrix.ncols(),
            last_row.into_owned().as_slice(),
        );

        let lower_half = &last_row_matrix * &tmp_matrix_b;

        let total_rows = upper_half.nrows() + lower_half.nrows();
        let mut next_extraction = DMatrix::zeros(total_rows, upper_half.ncols());

        next_extraction
            .view_mut((0, 0), upper_half.shape())
            .copy_from(&upper_half);

        next_extraction
            .view_mut((upper_half.nrows(), 0), lower_half.shape())
            .copy_from(&lower_half);

        extraction_matrix = next_extraction;
    }

    Ok(extraction_matrix)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Test that the extraction matrix has the correct (degree + 1) x (degree + 1) shape
    fn test_compute_knot_insertion_matrix_dimension() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        let degree = 2;
        let segment_index = 3;

        let matrix = compute_knot_insertion_matrix(&knots, degree, segment_index)
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
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let degree = 1;

        let invalid_segment_index = 3;

        let result = compute_knot_insertion_matrix(&knots, degree, invalid_segment_index);

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
