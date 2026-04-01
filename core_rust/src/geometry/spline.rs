use crate::core::knot::KnotVector;
use crate::core::point::WeightedPoint;
use crate::geometry::bezier::BezierCurve;
use crate::traits::{ParametricCurve, PointTrait};
use nalgebra::Point3;

//
// Main SplineCurve Struct
//
#[derive(Debug)]
pub struct SplineCurve<P: PointTrait> {
    pub(crate) degree: usize,
    pub(crate) controle_points: Vec<P>,
    pub(crate) knots: KnotVector,
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
        controle_points: Vec<Point3<f64>>,
        knots: KnotVector,
    ) -> Result<SplineCurve<Point3<f64>>, String> {
        self.validate(controle_points.len(), &knots)?;
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
    ) -> Result<SplineCurve<WeightedPoint>, String> {
        self.validate(ctrl_points.len(), &knots)?;
        if ctrl_points.len() != ctrl_point_weight.len() {
            return Err(format!(
                "Weight count mismatch: {} points vs {} weights",
                ctrl_points.len(),
                ctrl_point_weight.len()
            ));
        }
        let controle_points = ctrl_points
            .into_iter()
            .zip(ctrl_point_weight)
            .map(|(p, w)| WeightedPoint(p, w))
            .collect();
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
impl<P: PointTrait> SplineCurve<P> {
    pub fn builder() -> SplineCurveBuilder {
        SplineCurveBuilder::default()
    }

    fn _cox_de_boor(&self, i: usize, degree: usize, u: f64) -> Result<f64, String> {
        let n = self.controle_points.len() - 1;

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
                first_part = (u - self.knots.as_slice()[i]) / denom1
                    * self._cox_de_boor(i, degree - 1, u)?;
            }
            if i + degree + 1 < n {
                let denom2 = self.knots.as_slice()[i + degree + 1] - self.knots.as_slice()[i + 1];
                if denom2 != 0.0 {
                    second_part = (self.knots.as_slice()[i + degree + 1] - u) / denom2
                        * self._cox_de_boor(i + 1, degree - 1, u)?;
                }
            }
        }
        Ok(first_part + second_part)
    }

    pub fn converte_to_bezier(&self, _sample: usize) -> Result<Vec<BezierCurve<P>>, String> {
        // This logic is still incomplete, returning an empty vec for now.
        todo!()
    }
}

//
// ParametricCurve Trait Implementation
//
impl<P: PointTrait> ParametricCurve for SplineCurve<P> {
    fn domain(&self) -> (f64, f64) {
        let n = self.controle_points.len() - 1;
        let p = self.degree;
        (self.knots.as_slice()[p], self.knots.as_slice()[n + 1])
    }
}
