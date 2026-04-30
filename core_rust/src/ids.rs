use std::fmt;

use uuid::Uuid;

/// Stable string ID for curves, suitable for cross-process addressing.
///
/// Format: `curve-<uuid>`
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct CurveId(String);

impl CurveId {
    pub fn new() -> Self {
        Self(format!("curve-{}", Uuid::new_v4()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Debug for CurveId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("CurveId").field(&self.0).finish()
    }
}

impl fmt::Display for CurveId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Stable string ID for curve control points.
///
/// Format: `curve-<uuid>.cp-<uuid>`
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ControlPointId(String);

impl ControlPointId {
    pub fn new(curve_id: &CurveId) -> Self {
        Self(format!("{}.cp-{}", curve_id.as_str(), Uuid::new_v4()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Debug for ControlPointId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("ControlPointId").field(&self.0).finish()
    }
}

impl fmt::Display for ControlPointId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_id_format_is_prefixed() {
        let id = CurveId::new();
        assert!(id.as_str().starts_with("curve-"));
    }

    #[test]
    fn control_point_id_is_namespaced_by_curve_id() {
        let curve = CurveId::new();
        let cp = ControlPointId::new(&curve);
        assert!(cp.as_str().starts_with(&format!("{}.cp-", curve.as_str())));
    }
}

