/// Defines the behavior of a parametric curve.
pub trait ParametricCurve {
    /// Returns the valid parameter range, or domain, of the curve.
    /// For a standard clamped NURBS curve, this is typically `(0.0, 1.0)`.
    fn domain(&self) -> (f64, f64);
}
