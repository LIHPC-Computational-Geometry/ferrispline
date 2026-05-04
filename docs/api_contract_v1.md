# Ferrispline API Contract v1 (Core Rust + Python `nurbslib`)


## 0. Scope and layers

- **Core Rust (workspace member `core_rust/`)**
  - Owns all persistent geometry data (single source of truth).
  - Implements algorithms, validation, mutation, caching, and undo/redo primitives.
- **Python bindings (`nurbslib/`, PyO3)**
  - Exposes a Python-first API to the core.
  - Uses NumPy arrays for inputs/outputs.
  - Prefer **zero-copy output** when possible (`IntoPyArray`) and avoid extra conversions.
- **Bot**
  - Must not store authoritative curve data.
  - Only holds ephemeral render meshes and UI state.

---

## 1. IDs (string) — stable cross-process addressing

### 1.1 ID types

All IDs are **strings**.

- `CurveId`: `curve-<uuid>`
- `ControlPointId`: `curve-<uuid>.cp-<cp_uuid_or_index>`
- `KnotId` (optional v1): `curve-<uuid>.knot-<knot_uuid_or_index>`

### 1.2 ID stability policy (decision)

We must choose **one** of these strategies for CP/knots:

- **Strategy S1 (index-based IDs)**: `cp-0`, `cp-1`, ...
  - Pros: simple.
  - Cons: IDs change when inserting/deleting control points or knots → UI selections can become stale.
- **Strategy S2 (UUID-based CP IDs)**: CPs/knots carry stable UUIDs; ordering is separate.
  - Pros: stable selection, robust undo/redo and insertion/deletion.
  - Cons: slightly more bookkeeping.

**Decision for v1**: adopt **S2 (UUID-based child IDs)** for future-proofing.
This aligns with cross-process selection stability and undo/redo requirements.

---

## 2. Data model (minimal v1)

### 2.1 Curve kinds

Ferrispline v1 supports:

- **BezierCurve** (rational optional via weights)
- **SplineCurve** (B-spline / NURBS via knots + weights)

### 2.2 Storage model

Core provides a store (name not mandated yet) that:

- creates/removes curves
- resolves curves and components by **ID**
- invalidates and serves baked geometry
- applies undo/redo commands (command pattern)

---

## 3. Mutability contract (Rust `&self` vs `&mut self`)

### 3.1 Pure methods (`&self`)

Must be side-effect free:

- evaluate / sample
- read control points / weights / knots / degree
- domain
- conversion helpers (e.g. `to_bezier`) as **pure** transforms
- baked getters **if cache hit** (cache read is allowed but must not change geometry)

### 3.2 Mutating methods (`&mut self`)

These must:

- validate invariants
- mutate geometry
- mark baked geometry as **dirty**
- be compatible with undo/redo (either via commands or returning “before” data)

Examples:

- move control point
- set weight(s)
- knot insertion / knot move
- degree elevation
- delete operations (optional v1)

---

## 4. NumPy formats (canonical shapes and dtypes)

### 4.1 Coordinate conventions

- Coordinates are 3D cartesian: `[x, y, z]`.
- Use `float64` internally in Rust for math; expose `float64` in Python for API correctness.
- For **render/baked buffers**, expose `float32` by default (GPU friendly), unless precision matters.

### 4.2 Canonical shapes

#### Control points

- **Positions**: `np.ndarray[np.float64] shape (N, 3)`
- **Weights**: `np.ndarray[np.float64] shape (N,)`

#### Knot vector

- `np.ndarray[np.float64] shape (M,)` or `list[float]` (v1 may accept Python list for convenience)

#### Evaluated curve samples (math output)

Pick one canonical shape (must be consistent across all curves):

- **Option A (row-major points)**: `(sample, 3)`
- **Option B (axis-major)**: `(3, sample)`

**Contract v1 recommendation**: use **Option A `(sample, 3)`** in the public API because it matches
typical Python expectations and direct consumption in render pipelines.

> Note: existing implementation currently returns shapes that may differ (e.g. `(3, sample)` in some paths).
> A1/A2 must normalize this.

We use the option A.

#### Baked render buffers

Minimum baked payload for each curve:

- `vertices`: `np.ndarray[np.float32] shape (V, 3)`  (polyline sample points)
- `indices` or `edges`:
  - `edges`: `np.ndarray[np.int32] shape (E, 2)` for line segments **or**
  - `indices`: `np.ndarray[np.int32] shape (I,)` for indexed line list
- `bounds`: `min/max/center/size` as `np.ndarray[np.float32] shape (3,)` or Python lists

Optional baked payload (for interactive editing):

- `control_points`: `np.ndarray[np.float32] shape (N, 3)` (handle positions)
- `control_point_ids`: `list[str]` length N (stable IDs)

---

## 5. Python API surface (v1) — signatures

> Names are indicative. The contract is the *behavior and types*.

### 5.1 Model / store

- `class nurbslib.Model:`
  - `create_bezier(degree: int, control_points: np.ndarray[(N,3), float64], weights: Optional[np.ndarray[(N,), float64]] = None) -> str`
  - `create_nurbs(degree: int, control_points: np.ndarray[(N,3), float64], weights: np.ndarray[(N,), float64], knots: np.ndarray[(M,), float64] | list[float]) -> str`
  - `delete_curve(curve_id: str) -> None`

  - `get_curve_kind(curve_id: str) -> str  # "bezier" | "nurbs"`
  - `get_degree(curve_id: str) -> int`
  - `get_control_points(curve_id: str) -> np.ndarray[(N,3), float64]`
  - `get_weights(curve_id: str) -> np.ndarray[(N,), float64]`
  - `get_knots(curve_id: str) -> np.ndarray[(M,), float64]`  (nurbs only)

  - `evaluate(curve_id: str, sample: int, rational: Optional[bool] = None) -> np.ndarray[(sample,3), float64]`

  - `move_control_point(curve_id: str, control_point_id: str | int, new_pos: np.ndarray[(3,), float64]) -> None`
  - `set_weight(curve_id: str, control_point_id: str | int, weight: float) -> None`
  - `degree_elevation(curve_id: str, new_degree: int) -> None`  (bezier v1, nurbs optional)

  - `get_baked(curve_id: str, sample: int) -> dict`
    - returns the baked payload described in section 4.2

### 5.2 Undo / redo

- `model.undo() -> bool` (returns False if nothing to undo)
- `model.redo() -> bool`

### 5.3 Optional “object wrappers” (Python ergonomics)

Optionally expose lightweight wrappers that store only the `curve_id` and forward calls to a shared `Model`.
This is compatible with Bot’s “wrapper without data” direction.

---

## 6. Error contract (Python)

Bindings must map core errors into consistent Python exceptions.

- **Unknown IDs**: `KeyError`
  - unknown `curve_id`
  - unknown `control_point_id`
- **Invalid values**: `ValueError`
  - degree mismatch
  - knots not non-decreasing
  - wrong array shapes
  - invalid sample count (`sample <= 1`)
- **Unsupported operation** (v1 constraints): `NotImplementedError`
  - deletions if not supported yet
  - knot moving if not supported yet

Error messages should include:

- the failing ID
- the expected constraints (shape, range, etc.)

---

## 7. Performance contract (v1)

- `evaluate(...)` is allowed to allocate its output.
- `get_baked(...)` should reuse cached baked buffers when geometry is unchanged and `sample` matches.
- Python outputs should be:
  - **zero-copy** when consuming an owned Rust `ndarray::Array` via `IntoPyArray`
  - avoid copies on the input path when reading `PyReadonlyArray` views, except where unavoidable.

---

## 8. Gap analysis vs current repo (informational)

Already present today:

- `PyBezierCurve` / `PySplineCurve` exist in `nurbslib`.
- Bezier already distinguishes pure vs mutating methods (`evaluate(&self)` vs `degree_elevation(&mut self)`).

Missing for v1 contract:

- Global **Model / store** that owns multiple curves and addresses them by ID.
- End-to-end **string IDs** (curve + components).
- Standard **baked geometry** payload (vertices/edges/bounds + CP handles).
- Command-pattern undo/redo at model/store level.
