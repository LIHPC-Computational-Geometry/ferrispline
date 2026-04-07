use crate::core::knot::KnotVector;
use crate::traits::ParametricCurve;
use nalgebra::{Point3, Vector3};
use ndarray::linspace;

//
// Main SplineCurve Struct
//
#[derive(Debug)]
pub struct SplineCurve {
    pub degree: usize,
    pub controle_points: Vec<(Point3<f64>, f64)>,
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
        ctrl_pts: Vec<Point3<f64>>,
        knots: KnotVector,
    ) -> Result<SplineCurve, String> {
        self.validate(ctrl_pts.len(), &knots)?;
        let controle_points = ctrl_pts.into_iter().map(|pt| (pt, 1.0)).collect();
        Ok(SplineCurve {
            degree: self.degree,
            controle_points,
            knots,
        })
    }

    pub fn build_nurbs(
        self,
        ctrl_points: Vec<Point3<f64>>,
        ctrl_point_weight: Vec<f64>,
        knots: KnotVector,
    ) -> Result<SplineCurve, String> {
        self.validate(ctrl_points.len(), &knots)?;
        if ctrl_points.len() != ctrl_point_weight.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                ctrl_points.len(),
                ctrl_point_weight.len()
            ));
        }
        let controle_points = ctrl_points.into_iter().zip(ctrl_point_weight).collect();
        Ok(SplineCurve {
            degree: self.degree,
            controle_points,
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

    pub fn eval_nurbs_curve(&self, sample: usize) -> Result<Vec<Point3<f64>>, String> {
        let domain = self.domain();
        let u_vals = linspace(domain.0, domain.1, sample);

        let mut points: Vec<Point3<f64>> = Vec::new();
        for u in u_vals {
            let mut numerator = Vector3::zeros();
            let mut denominator = 0.0;

            for i in 0..self.controle_points.len() {
                let n = self.cox_de_boor(i, self.degree, u)?;
                let weight_n = self.controle_points[i].1 * n;
                numerator += self.controle_points[i].0.coords * weight_n;
                denominator += weight_n;
            }
            let point = if denominator.abs() < 1e-9 {
                Vector3::zeros()
            } else {
                numerator / denominator
            };

            points.push(Point3::from(point));
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
