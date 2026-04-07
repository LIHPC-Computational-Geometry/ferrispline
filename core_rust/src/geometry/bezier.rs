use nalgebra::{DMatrix, Point3};
use ndarray::{Array1, linspace};
use num_integer::binomial;

use crate::core::controle_point::Point3Ext;

#[derive(Debug)]
pub struct BezierCurve {
    pub degree: usize,
    pub weights: Array1<f64>,
    pub controle_points: Array1<Point3<f64>>,
}

impl BezierCurve {
    pub fn new_with_weights(
        degree: usize,
        controle_points: Array1<Point3<f64>>,
        weights: Array1<f64>,
    ) -> Self {
        Self {
            degree,
            weights,
            controle_points,
        }
    }

    pub fn new(degree: usize, controle_points: Array1<Point3<f64>>) -> Self {
        Self {
            degree,
            weights: Array1::from(vec![1.0; controle_points.len()]),
            controle_points,
        }
    }

    pub fn bernstein(&self, v: usize, t: Vec<f64>) -> Vec<f64> {
        t.iter()
            .map(|t| {
                binomial(self.degree, v) as f64
                    * t.powf(v as f64)
                    * (1.0 - t).powf((self.degree - v) as f64)
            })
            .collect()
    }

    // NOTE: Rester en dimension 1 avec uniquement des curve pour le moment
    // NOTE: Faire des fonctions de modification des courbes
    // a implémenter dans le projet de Frank
    /// Evaluate Bezier curve for a number of points equal to `sample`
    pub fn evalutate(&self, sample: usize) -> Array1<Point3<f64>> {
        let t: Vec<f64> = linspace(0.0, 1.0, sample).collect();
        let mut points: Array1<Point3<f64>> =
            Array1::from(vec![Point3::new(0.0, 0.0, 0.0); sample]);
        for i in 0..(self.degree + 1) {
            let forces = self.bernstein(i, t.clone());
            points = points
                .into_iter()
                .zip(forces)
                .map(|(point, force)| point.add_points(&(self.controle_points[i] * force)))
                .collect();
        }
        points
    }

    /// Evaluate Rationnal Bezier curve for a number of points equal to `sample`
    pub fn evaluate_rational(&self, _sample: usize) -> Vec<Point3<f64>> {
        todo!()
    }

    #[allow(dead_code)]
    fn rational_basis(&self) -> Vec<f64> {
        todo!()
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
        let start_idx = segment_index.saturating_sub(degree_step);
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
    use ndarray::array;

    use super::*;

    #[test]
    fn test_bernstein_extremes() {
        let degree = 3;

        // At t=0, only the first polynomial (v=0) is 1.0
        let curve_v0: BezierCurve = BezierCurve::new(degree, array![]);
        assert_eq!(curve_v0.bernstein(0, vec![0.0]), vec![1.0]);
        let curve_v1: BezierCurve = BezierCurve::new(degree, array![]);
        assert_eq!(curve_v1.bernstein(1, vec![0.0]), vec![0.0]);

        // At t=1, only the last polynomial (v=degree) is 1.0
        let curve: BezierCurve = BezierCurve::new(degree, array![]);
        assert_eq!(curve.bernstein(degree, vec![1.0]), vec![1.0]);
        let curve: BezierCurve = BezierCurve::new(degree, array![]);
        assert_eq!(curve.bernstein(0, vec![1.0]), vec![0.0]);
    }

    #[test]
    fn test_bernstein_partition_of_unity() {
        let degree = 3;

        let mut t_vals = Vec::new();
        for i in 0..10 {
            t_vals.push(i as f64 / 9.0);
        }

        let mut totals = vec![0.0; 10];

        for v in 0..=degree {
            let curve: BezierCurve = BezierCurve::new(degree, array![]);
            let results = curve.bernstein(v, t_vals.clone());
            for (i, &res) in results.iter().enumerate() {
                totals[i] += res;
            }
        }

        for total in totals {
            assert!((total - 1.0).abs() < 1e-10);
        }
    }

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

    #[test]
    fn test_bezier_evaluate_simple() {
        let degree = 2;
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 2.0, 0.0);
        let p2 = Point3::new(2.0, 0.0, 0.0);
        let control_points = array![p0, p1, p2];
        let curve = BezierCurve::new(degree, control_points);

        // Evaluate with 3 samples: t=0.0, t=0.5, t=1.0
        let points = curve.evalutate(3);

        assert_eq!(points.len(), 3);

        // At t=0.0, point should be p0
        assert!((points[0].x - 0.0).abs() < 1e-6);
        assert!((points[0].y - 0.0).abs() < 1e-6);
        assert!((points[0].z - 0.0).abs() < 1e-6);

        // At t=0.5, point should be 0.25*p0 + 0.5*p1 + 0.25*p2 = (1.0, 1.0, 0.0)
        assert!((points[1].x - 1.0).abs() < 1e-6);
        assert!((points[1].y - 1.0).abs() < 1e-6);
        assert!((points[1].z - 0.0).abs() < 1e-6);

        // At t=1.0, point should be p2
        assert!((points[2].x - 2.0).abs() < 1e-6);
        assert!((points[2].y - 0.0).abs() < 1e-6);
        assert!((points[2].z - 0.0).abs() < 1e-6);
    }
}
