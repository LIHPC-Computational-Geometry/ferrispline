use nalgebra::{Point3, Vector4};

use crate::traits::PointTrait;

impl PointTrait for Point3<f64> {
    fn to_homogeneous(self) -> Vector4<f64> {
        Vector4::new(self.x, self.y, self.z, 1.0)
    }

    fn weights(self) -> f64 {
        1.0
    }

    fn zeros() -> Self {
        Point3::new(0.0, 0.0, 0.0)
    }
}

#[derive(Debug, Copy, Clone)]
pub struct WeightedPoint(pub Point3<f64>, pub f64);

impl PointTrait for WeightedPoint {
    fn to_homogeneous(self) -> Vector4<f64> {
        Vector4::new(
            self.0.x * self.1,
            self.0.y * self.1,
            self.0.z * self.1,
            self.1,
        )
    }

    fn weights(self) -> f64 {
        self.1
    }

    fn zeros() -> Self {
        WeightedPoint(Point3::zeros(), 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_homogeneous() {
        let point = Point3::new(1.0, 2.0, 3.0);
        let weighted_point = WeightedPoint(point, 4.0);

        let homogeneous_point = point.to_homogeneous();
        let homogeneous_weighted_point = weighted_point.to_homogeneous();
        assert_eq!(homogeneous_point, Vector4::new(1.0, 2.0, 3.0, 1.0));
        assert_eq!(
            homogeneous_weighted_point,
            Vector4::new(4.0, 8.0, 12.0, 4.0)
        );
    }

    #[test]
    fn test_weights() {
        let point = Point3::new(1.0, 2.0, 3.0);
        let weighted_point = WeightedPoint(point, 4.0);

        assert_eq!(point.weights(), 1.0);
        assert_eq!(weighted_point.weights(), 4.0);
    }
}
