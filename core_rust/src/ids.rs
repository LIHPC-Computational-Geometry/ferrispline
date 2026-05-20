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

    pub fn try_from_str(s: &str) -> Result<Self, String> {
        let prefix = "curve-";
        if let Some(uuid_str) = s.strip_prefix(prefix) {
            if Uuid::parse_str(uuid_str).is_ok() {
                return Ok(Self(s.to_string()));
            }
        }
        Err(format!("Invalid CurveId: {}", s))
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

impl Default for CurveId {
    fn default() -> Self {
        Self::new()
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

    #[test]
    fn try_from_str_valid_curve_id() {
        let id = CurveId::new();
        let parsed = CurveId::try_from_str(id.as_str());
        assert!(parsed.is_ok());
        assert_eq!(id, parsed.unwrap());
    }

    #[test]
    fn try_from_str_invalid_curve_id() {
        assert!(CurveId::try_from_str("invalid-prefix").is_err());
        assert!(CurveId::try_from_str("curve-invalid-uuid").is_err());
    }
}
