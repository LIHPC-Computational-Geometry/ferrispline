use crate::geometry::bezier::BezierCurve;
use ndarray::{Array1, Array2};

impl BezierCurve {
    /// Subdivide the Bezier curve at parameter `t` (0..=1).
    ///
    /// This uses de Casteljau in homogeneous coordinates, so it works for both
    /// polynomial and rational Bezier curves (weights).
    pub fn subdivide(&self, t: f64) -> Result<(BezierCurve, BezierCurve), String> {
        if !(0.0..=1.0).contains(&t) {
            return Err(format!("subdivide: t must be within [0,1], got {t}"));
        }

        let n = self.degree;
        if self.control_points.ncols() != 3 {
            return Err(format!(
                "subdivide: expected control_points shape (N,3), got (N,{})",
                self.control_points.ncols()
            ));
        }

        // Homogeneous control points: [x*w, y*w, z*w, w]
        let mut level: Array2<f64> =
            BezierCurve::homogeneous(&self.control_points, &self.weights, n);

        let mut left: Array2<f64> = Array2::zeros((n + 1, 4));
        let mut right: Array2<f64> = Array2::zeros((n + 1, 4));

        left.row_mut(0).assign(&level.row(0));
        right.row_mut(n).assign(&level.row(n));

        for r in 1..=n {
            for i in 0..=(n - r) {
                for k in 0..4 {
                    let a = level[[i, k]];
                    let b = level[[i + 1, k]];
                    level[[i, k]] = (1.0 - t) * a + t * b;
                }
            }
            left.row_mut(r).assign(&level.row(0));
            right.row_mut(n - r).assign(&level.row(n - r));
        }

        let (left_pts, left_w) = BezierCurve::dehomogenize(&left)?;
        let (right_pts, right_w) = BezierCurve::dehomogenize(&right)?;

        let left_curve = BezierCurve::new_with_weights(n, left_pts, left_w)?;
        let right_curve = BezierCurve::new_with_weights(n, right_pts, right_w)?;
        Ok((left_curve, right_curve))
    }

    /// Subdivide the curve into `segments` uniform parameter spans.
    ///
    /// Note: uniform in parameter space, not arc-length.
    pub fn subdivide_uniform(&self, segments: usize) -> Result<Vec<BezierCurve>, String> {
        if segments == 0 {
            return Err("subdivide_uniform: segments must be >= 1".to_string());
        }
        if segments == 1 {
            return Ok(vec![self.clone()]);
        }

        // Repeatedly subdivide at a normalized parameter for the remaining curve.
        let mut out: Vec<BezierCurve> = Vec::with_capacity(segments);
        let mut cur = self.clone();
        for i in 1..segments {
            // We want to cut the original [0,1] into segments parts.
            // Given that `cur` represents [i-1, 1] (in original space), the next cut
            // is at t = 1/(segments-(i-1)) in the current curve's parameter.
            let remaining = segments - (i - 1);
            let t = 1.0 / (remaining as f64);
            let (a, b) = cur.subdivide(t)?;
            out.push(a);
            cur = b;
        }
        out.push(cur);
        Ok(out)
    }

    fn homogeneous(pts: &Array2<f64>, w: &Array1<f64>, n: usize) -> Array2<f64> {
        let mut level: Array2<f64> = Array2::zeros((n + 1, 4));
        for i in 0..=n {
            let w = w[i];
            level[[i, 0]] = pts[[i, 0]] * w;
            level[[i, 1]] = pts[[i, 1]] * w;
            level[[i, 2]] = pts[[i, 2]] * w;
            level[[i, 3]] = w;
        }
        level
    }

    fn dehomogenize(h: &Array2<f64>) -> Result<(Array2<f64>, Array1<f64>), String> {
        let n = h.nrows();
        let mut pts: Array2<f64> = Array2::zeros((n, 3));
        let mut w: Array1<f64> = Array1::zeros(n);
        for i in 0..n {
            let wi = h[[i, 3]];
            if wi.abs() < 1e-12 {
                return Err(
                    "subdivide: zero weight encountered during dehomogenization".to_string()
                );
            }
            w[i] = wi;
            pts[[i, 0]] = h[[i, 0]] / wi;
            pts[[i, 1]] = h[[i, 1]] / wi;
            pts[[i, 2]] = h[[i, 2]] / wi;
        }
        Ok((pts, w))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn subdivide_preserves_endpoints_polynomial() {
        let degree = 3;
        let ctrl = array![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [3.0, 0.0, 0.0]
        ];
        let curve = BezierCurve::new(degree, ctrl.clone()).unwrap();

        let (l, r) = curve.subdivide(0.37).unwrap();

        // left start == original start
        for k in 0..3 {
            assert!((l.control_points[[0, k]] - ctrl[[0, k]]).abs() < 1e-10);
        }
        // right end == original end
        for k in 0..3 {
            assert!((r.control_points[[degree, k]] - ctrl[[degree, k]]).abs() < 1e-10);
        }
    }

    #[test]
    fn subdivide_uniform_counts_segments() {
        let degree = 2;
        let ctrl = array![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(degree, ctrl).unwrap();
        let segs = curve.subdivide_uniform(5).unwrap();
        assert_eq!(segs.len(), 5);
        for s in segs {
            assert_eq!(s.degree, degree);
        }
    }
}
