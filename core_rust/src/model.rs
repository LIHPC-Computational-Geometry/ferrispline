use std::collections::HashMap;

use ndarray::{Array1, Array2};

use crate::geometry::{bezier::BezierCurve, spline::SplineCurve};
use crate::ids::CurveId;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    Bezier,
    Nurbs,
}

#[derive(Debug)]
pub enum ModelError {
    CurveNotFound { curve_id: CurveId },
    MutationFailed { curve_id: CurveId, message: String },
    MutationsFailed { curves_id: Array1<CurveId>, message: String },
}

/// A curve owned by the model store.
#[derive(Debug)]
pub enum Curve {
    Bezier(BezierCurve),
    Nurbs(SplineCurve),
}

impl Curve {
    /// Return the curve's kind.
    pub fn kind(&self) -> CurveKind {
        match self {
            Curve::Bezier(_) => CurveKind::Bezier,
            Curve::Nurbs(_) => CurveKind::Nurbs,
        }
    }

    /// Pure evaluation. Output shape is currently delegated to underlying implementations.
    pub fn evaluate(&self, sample: usize) -> Result<Array2<f64>, String> {
        match self {
            Curve::Bezier(c) => Ok(c.evaluate(sample)),
            Curve::Nurbs(c) => c.eval_nurbs_curve(sample),
        }
    }

    pub fn to_bezier(&self) -> Result<Vec<BezierCurve>, String> {
        match self {
            Curve::Bezier(_c) => todo!("fix: clone() with '&'"),
            Curve::Nurbs(c) => c.to_bezier(),
        }
    }

    pub fn to_nurbs(&self) -> Result<SplineCurve, String> {
        match self {
            Curve::Bezier(_c) => todo!(),
            Curve::Nurbs(_c) => todo!("fix: clone() with '&'"),
        }
    }
}

/// A wrapper around a curve that tracks its modification state.
/// The `dirty` flag is set to `true` whenever the curve is modified.
/// This simply tells the system: "I have changed! My geometry needs to be
/// recalculated before the next time you draw me."
#[derive(Debug)]
struct CurveEntry {
    curve: Curve,
    dirty: bool,
}

/// In-memory store of editable curves.
#[derive(Debug, Default)]
pub struct Model {
    curves: HashMap<CurveId, CurveEntry>,
}

impl Model {
    pub fn new() -> Self {
        Self {
            curves: HashMap::new(),
        }
    }

    // -----------------------------
    // Creation / deletion
    // -----------------------------

    pub fn create_bezier(&mut self, degree: usize, control_points: Array2<f64>) -> Result<CurveId, String> {
        let curve = BezierCurve::new(degree, control_points)?;
        let id = CurveId::new();
        self.curves.insert(
            id.clone(),
            CurveEntry {
                curve: Curve::Bezier(curve),
                dirty: true,
            },
        );
        Ok(id)
    }

    pub fn create_bezier_with_weights(
        &mut self,
        degree: usize,
        control_points: Array2<f64>,
        weights: Array1<f64>,
    ) -> Result<CurveId, String> {
        let curve = BezierCurve::new_with_weights(degree, control_points, weights)?;
        let id = CurveId::new();
        self.curves.insert(
            id.clone(),
            CurveEntry {
                curve: Curve::Bezier(curve),
                dirty: true,
            },
        );
        Ok(id)
    }

    pub fn insert_nurbs(&mut self, id: CurveId, curve: SplineCurve) {
        self.curves.insert(
            id,
            CurveEntry {
                curve: Curve::Nurbs(curve),
                dirty: true,
            },
        );
    }

    pub fn delete_curve(&mut self, curve_id: &CurveId) -> bool {
        self.curves.remove(curve_id).is_some()
    }

    // -----------------------------
    // Read-only access (pure)
    // -----------------------------

    pub fn curve_kind(&self, curve_id: &CurveId) -> Result<CurveKind, ModelError> {
        Ok(self.get_curve(curve_id)?.kind())
    }

    pub fn evaluate(&self, curve_id: &CurveId, sample: usize) -> Result<Array2<f64>, String> {
        self.get_curve(curve_id)
            .map_err(|e| format!("{e:?}"))?
            .evaluate(sample)
    }

    pub fn is_dirty(&self, curve_id: &CurveId) -> Result<bool, ModelError> {
        Ok(self
            .curves
            .get(curve_id)
            .ok_or_else(|| ModelError::CurveNotFound {
                curve_id: curve_id.clone(),
            })?
            .dirty)
    }

    // -----------------------------
    // Mutating access (marks dirty)
    // -----------------------------


    /// Entry point for modifying a curve.
    /// This function locates the curve to be modified, calls the 'f' closure, that performs the modification,
    /// and indicates that the curve has been modified using the 'dirty' variable.
    ///
    /// Exemple:
    /// ```
    ///  pub fn move_point_on_curve(&mut self, curve_id: &CurveId, index: usize, new_pos: Array1<f64>) -> Result<(), ModelError> {
    ///    self.with_curve_mut(curve_id, |curve| {
    ///        curve.move_control_point(index, new_pos)
    ///    })
    /// }
    /// ```
    // NOTE: `|curve|` correspond à la variable entry dans la fonction `with_curve_mut`
    pub fn with_curve_mut<R>(
        &mut self,
        curve_id: &CurveId,
        f: impl FnOnce(&mut Curve) -> Result<R, String>,
    ) -> Result<R, ModelError> {
        let entry = self
            .curves
            .get_mut(curve_id)
            .ok_or_else(|| ModelError::CurveNotFound {
                curve_id: curve_id.clone(),
            })?;

        let out = f(&mut entry.curve).map_err(|message| ModelError::MutationFailed {
            curve_id: curve_id.clone(),
            message,
        })?;
        entry.dirty = true;
        Ok(out)
    }

    pub fn with_curves_mut<R>(&mut self, curves_id: &Array1<CurveId>,f: impl FnOnce(Vec<(CurveId, &mut CurveEntry)>) -> Result<R, String>) -> Result<R, ModelError> {
        let ids_to_fetch: Vec<CurveId> = curves_id.to_vec();

        let entries: Vec<(CurveId, &mut CurveEntry)> = self.curves.iter_mut()
            .filter(|(id, _)| ids_to_fetch.contains(id))
            .map(|(id, entry)| {
                entry.dirty = true;
                (id.clone(), entry)
            })
            .collect();

        let out = f(entries).map_err(|message| {
            ModelError::MutationsFailed {
                curves_id: curves_id.clone().to_owned(), message,
             }
            })?;
        Ok(out)
    }

    pub fn clear_dirty(&mut self, curve_id: &CurveId) -> Result<(), ModelError> {
        let entry = self
            .curves
            .get_mut(curve_id)
            .ok_or_else(|| ModelError::CurveNotFound {
                curve_id: curve_id.clone(),
            })?;
        entry.dirty = false;
        Ok(())
    }

    fn get_curve(&self, curve_id: &CurveId) -> Result<&Curve, ModelError> {
        Ok(&self
            .curves
            .get(curve_id)
            .ok_or_else(|| ModelError::CurveNotFound {
                curve_id: curve_id.clone(),
            })?
            .curve)
    }

    pub fn set_degree(&mut self, curve_id: &CurveId, _degree: usize) -> Result<(), ModelError> {
        self.with_curve_mut(curve_id, |curve| {
            match curve {
                Curve::Bezier(_curve) => {
                    todo!("create a function set_degree with degree elevation and degree reduction")
                }
                Curve::Nurbs(_curve) => {
                    todo!("create a function set_degree with degree elevation and degree reduction")
                }
            }
        })
    }

    pub fn convert(&mut self, curves_id: &Array1<CurveId>, new_kind: CurveKind) -> Result<(), ModelError> {
        self.with_curves_mut(curves_id, |_curves| {
            match new_kind {
                CurveKind::Bezier => {
                    todo!("curve.to_bezier()")
                }
                CurveKind::Nurbs => {
                    todo!("curves[0].to_nurbs(cuvres) or curve.to_nurbs()")
                }
            }
        })
    }

    pub fn move_control_point(&mut self, curve_id: &CurveId, _index: usize, _new_pos: Array1<f64>) -> Result<(), ModelError> {
        self.with_curve_mut(curve_id, |curve| {
            match curve {
                Curve::Bezier(_curve) => {
                    todo!("curve.move_control_point(index, new_pos)")
                }
                Curve::Nurbs(_curve) => {
                    todo!("curve.move_control_point(index, new_pos)")
                }
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn can_create_and_evaluate_bezier_through_model() {
        let mut model = Model::new();
        let ctrl = array![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let id = model.create_bezier(2, ctrl).unwrap();

        let pts = model.evaluate(&id, 5).unwrap();
        assert_eq!(pts.ncols(), 5);
        assert_eq!(pts.nrows(), 3);
        assert!(model.is_dirty(&id).unwrap());
    }

    #[test]
    fn delete_curve_returns_true_when_present() {
        let mut model = Model::new();
        let ctrl = Array2::<f64>::zeros((2, 3));
        let id = model.create_bezier(1, ctrl).unwrap();
        assert!(model.delete_curve(&id));
        assert!(!model.delete_curve(&id));
    }
}
