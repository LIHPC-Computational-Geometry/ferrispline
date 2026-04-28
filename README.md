# FerriSpline

## Overview
This repository hosts a high-performance library for generating, manipulating, and computing hexahedral meshes. It leverages a pure **Rust** core for intensive mathematical operations to ensure memory safety and optimal execution speed. It provides a native **Python API** via PyO3, designed for seamless integration into meshing workflows, blocking algorithms, and rendering pipelines (e.g., Panda3D, gmsh).

## Repository Structure
This project is organized as a Cargo Workspace (monorepo) to strictly separate the computational logic from the Python bindings.

```text
mon_projet_maillage/
├── Cargo.toml                 # Global workspace configuration
├── core_rust/                 # Pure Rust mathematical core
│   ├── Cargo.toml
│   └── src/lib.rs             # NURBS algorithms, mesh generation, geometry
├── nurbslib/                  # Rust-to-Python bindings (PyO3)
│   ├── Cargo.toml
│   ├── pyproject.toml         # Maturin build configuration
│   └── src/lib.rs             # Exposes `core_rust` functions to Python
└── sandbox_python/            # Sandbox python with the same content as `cargo_rust`
    ├── pyproject.toml
    ├── src/                   # Python logic (VTK conversion, nurbs conversion to bezier curve)
    └── tests/                 # Unit tests
```

## Prerequisites
- [Rust toolchain](https://rustup.rs/) (2024 edition recommended)
- Python 3.8 or higher
- [Maturin](https://www.maturin.rs/) (Build system for Rust-based Python extensions)

## Building for Development
To compile the Rust core and install the Python extension directly into your active virtual environment, run:

```bash
cd nurbslib
maturin develop
```

## Image for dev

One can run:

```bash
cd /path/to/ferrispline/
BRANCH_SLUG="$(git rev-parse --abbrev-ref HEAD | sed 's/[^a-zA-Z0-9._-]/-/g')"
IMAGE_NAME='ghcr.io/LIHPC-Computational-Geometry/ferrispline'
docker build --build-arg USER=${USER} \
             --network=host \
             --tag $IMAGE_NAME:$BRANCH_SLUG \
             --target dev ./
```

## Run tests

One can run:

```bash
xhost +local:docker
docker run --interactive \
           --network=host \
           --rm \
           --tty \
           --volume ./:/home/${USER}/ferrispline/ \
           $CI_REGISTRY_IMAGE:$CI_COMMIT_REF_SLUG
cd ferrispline/
pre-commit run
cargo test
./run.sh test
```
