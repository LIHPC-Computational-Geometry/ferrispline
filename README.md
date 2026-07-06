# FerriSpline

High-performance library for creating, evaluating, and manipulating NURBS and Bézier curves. A pure **Rust** computational core (`core_rust`) provides memory-safe geometry algorithms; a **Python** extension (`ferrispline`) exposes them through PyO3 and Maturin.

The sole public API is **`ferrispline.PyModel`**: a multi-curve store with stable string IDs and dirty-flag invalidation for downstream consumers.

## Features

**Implemented**

- Multi-curve store with UUID-based IDs (`curve-<uuid>`)
- Bézier curve evaluation (Bernstein basis), subdivision (De Casteljau), degree elevation
- NURBS/B-spline evaluation (Cox–De Boor, rational weighting)
- NURBS → Bézier segment extraction (Boehm knot-insertion matrix)
- Per-curve dirty tracking (`is_dirty`) for cache invalidation
- Zero-copy NumPy output on evaluation (`IntoPyArray`)

**In progress**

- Knot insertion and removal on NURBS curves
- Spline control-point mutation and degree change
- Bézier → NURBS conversion
- `clear_dirty` Python binding (available in Rust, not yet exposed)

## Quick start

```bash
cd python
maturin develop
```

```python
import numpy as np
import ferrispline

model = ferrispline.PyModel()
cp = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
curve_id = model.create_bezier(2, cp)
assert model.is_dirty(curve_id)
points = model.evaluate(curve_id, 50)  # shape (50, 3)
```

## Installation

### Prerequisites

- [Rust toolchain](https://rustup.rs/) (2024 edition)
- Python 3.8 or higher (the Makefile prefers `python3.14` when available)
- [Maturin](https://www.maturin.rs/)

For the sandbox visualisation environment, OpenGL libraries are required (`libgl1` on Debian/Ubuntu).

### Development build (fast iteration)

```bash
cd python
maturin develop
```

### Full build (release wheel + sandbox venv)

```bash
make build
```

This creates `sandbox_python/.venv`, compiles a release wheel with Maturin, and installs it into the venv.

## Running tests

```bash
make test
```

Runs `cargo test` (Rust core) and `pytest sandbox_python/tests` (Python sandbox).

Other Makefile targets: `make venv`, `make run <file.vtk> [samples]`, `make clean`, `make rebuild`.

## Repository structure

```text
ferrispline/
├── Cargo.toml              # Workspace: core_rust, python
├── core_rust/              # Pure Rust geometry engine
│   └── src/
│       ├── model.rs        # Multi-curve store, dirty tracking
│       ├── ids.rs          # CurveId, ControlPointId
│       └── geometry/       # bezier/, spline/
├── python/                 # PyO3 bindings → import ferrispline
│   ├── pyproject.toml      # Maturin configuration
│   └── src/
│       ├── lib.rs          # Module entry point
│       └── model.rs        # PyModel
├── sandbox_python/         # Reference Python impl + VTK demos
│   ├── src/nurbs_math/
│   └── tests/
├── docs/
│   └── technical_reference_v1.md
├── Makefile
├── Dockerfile
└── CONTRIBUTING.md
```

## Documentation

- [Technical Reference V1](docs/technical_reference_v1.md) — architecture, API, dirty invalidation, CI/CD
- [Contributing](CONTRIBUTING.md) — development workflow and quality gates

## Docker

Build a development image:

```bash
BRANCH_SLUG="$(git rev-parse --abbrev-ref HEAD | sed 's/[^a-zA-Z0-9._-]/-/g')"
docker build --build-arg USER="${USER}" \
             --network=host \
             --tag "ghcr.io/LIHPC-Computational-Geometry/ferrispline:${BRANCH_SLUG}" \
             --target dev .
```

The `latest` Docker target is used by CI; the `dev` target adds an interactive shell with Rust and Python tooling.

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
