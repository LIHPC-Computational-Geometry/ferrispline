# Contributing to FerriSpline

Thank you for your interest in contributing. This guide describes the development workflow, quality gates, and API policies for the project.

For architecture and API details, see [docs/technical_reference_v1.md](docs/technical_reference_v1.md).

## Getting set up

```bash
git clone <repository-url>
cd ferrispline
make build
make test
```

For fast iteration on the Python bindings:

```bash
cd python
maturin develop
```

Install pre-commit hooks (recommended):

```bash
pip install pre-commit
pre-commit install
```

## Branch naming

Use descriptive branch prefixes:

| Prefix | Use case |
|---|---|
| `feature/` | New functionality |
| `fix/` | Bug fixes |
| `docs/` | Documentation only |
| `test/` | Test additions or fixes |
| `refactor/` | Code restructuring without behaviour change |

Example: `feature/clear-dirty-binding`

## Commit style

Follow [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` — new feature
- `fix:` — bug fix
- `docs:` — documentation
- `test:` — tests
- `refactor:` — refactoring
- `chore:` — tooling, dependencies
- `perf:` — performance improvement

Example: `feat: expose clear_dirty on PyModel`

## Pre-submit checklist

Before opening a pull request:

- [ ] `pre-commit run --all-files` passes
- [ ] `make test` passes (`cargo test` + `pytest sandbox_python/tests`)
- [ ] Documentation updated if the public API changed
- [ ] New features exposed through `PyModel` only (see API policy below)
- [ ] Stubbed operations are not documented as available

## CI expectations

Every push and pull request to `main` triggers `.github/workflows/ci-pipeline.yml`:

1. **Build** — Docker image built and pushed to `ghcr.io/<repository>:<branch-slug>` (target `latest`)
2. **Quality and tests** — runs inside the pushed image:
   - `pre-commit run --all-files`
   - `make test`

Your pull request must pass both jobs before merge.

## Public API policy

`ferrispline.PyModel` is the sole public API. All new functionality must be added as methods on `PyModel` in `python/src/model.rs`, backed by `core_rust::model::Model`.

Do **not**:

- Add new public methods to `PyBezierCurve` or `PySplineCurve`
- Document standalone curve classes as supported integration points
- Describe unimplemented (`todo!()`) operations as available

Refer to the [implementation status table](docs/technical_reference_v1.md#15-implementation-status) in the technical reference when documenting or implementing features.

## Code quality

### Rust (`core_rust/`, `python/`)

- Format with `cargo fmt`
- Lint with `clippy`
- Add `#[test]` modules for new algorithms and store logic
- Prefer pure reads (`&self`) for evaluation and queries; mutations (`&mut self`) must set the dirty flag via `with_curve_mut`

### Python (`sandbox_python/`)

- Lint and format with `ruff` (enforced by pre-commit)
- Add tests under `sandbox_python/tests/` for reference implementations and integration scenarios

## Pull request process

1. Fork the repository and create a branch from `main`.
2. Make your changes following the checklist above.
3. Open a pull request with a clear description of the change and its motivation.
4. Ensure CI passes.
5. A maintainer will review and provide feedback.

Do not push directly to `main`. The `no-commit-to-branch` pre-commit hook enforces this locally.

## Questions

For architectural decisions and API design, consult [docs/technical_reference_v1.md](docs/technical_reference_v1.md) first. Open an issue for discussion if the reference does not cover your case.
