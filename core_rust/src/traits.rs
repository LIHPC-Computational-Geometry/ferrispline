use nalgebra::Vector4;

/// Defines the behavior of a parametric curve.
pub trait ParametricCurve {
    /// Returns the valid parameter range, or domain, of the curve.
    /// For a standard clamped NURBS curve, this is typically `(0.0, 1.0)`.
    fn domain(&self) -> (f64, f64);
}

pub trait PointTrait: Copy + Clone {
    /// Converts the point to homogeneous coordinates.
    fn to_homogeneous(self) -> Vector4<f64>;
    fn weights(self) -> f64;
    fn zeros() -> Self;
}
