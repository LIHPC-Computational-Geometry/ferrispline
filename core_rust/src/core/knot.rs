#[derive(Debug, PartialEq, Clone)]
pub struct KnotVector(pub Vec<f64>);

impl KnotVector {
    pub fn new(knots: Vec<f64>) -> Result<Self, String> {
        if knots.windows(2).any(|w| w[0] > w[1]) {
            return Err("Knots must be in non-decreasing order".to_string());
        }
        Ok(KnotVector(knots))
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.0
    }

    pub fn lenght_check(&self, num_control_points: usize, degree: usize) -> Result<&Self, String> {
        let n = num_control_points - 1;
        let m = self.as_slice().len() - 1;
        if m != n + degree + 1 {
            return Err(format!(
                "Knot vector length mismatch: m = n + p + 1 relation not hold. m={}, n={}, p={}",
                m, n, degree
            ));
        }
        Ok(self)
    }

    /// Calculates the multiplicity of a given knot value.
    pub fn multiplicity(&self, knot_value: f64) -> usize {
        self.as_slice()
            .iter()
            .filter(|&&k| (k - knot_value).abs() < 1e-9)
            .count()
    }
}
