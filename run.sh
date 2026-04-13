#!/bin/bash

set -e

echo "🟢 Step 1 : Check the virtual environnement..."

if [ -n "$VIRTUAL_ENV" ]; then
    echo "   -> The virtual environnement is already active ($VIRTUAL_ENV)."
elif [ -d ".venv" ]; then
    echo "   -> Any environnement detected. .venv activation..."
    source .venv/bin/activate
else
    echo "  -> No .venv detected. Using global system environement (CI mode)."
fi

echo "⚙️ Step 2 : Compilation of the Rust's librairy (nurbslib)..."
cd nurbslib

if [ -n "$VIRTUAL_ENV" ]; then
    maturin develop
else
    # NOTE: build the rust package (.so)
    maturin build --release --out dist
    # NOTE: Python can't ready .so file so the package is packaging into an archive standard file wheel (.whl)
    pip install --break-system-packages dist/*.whl --force-reinstall
fi
cd ..

echo "🚀 Step 3 : Launch of nurbs-convert..."
cd sandbox_python
if [ "$1" == "test" ]; then
    echo "Lancement des tests avec Pytest..."
    pytest
else
    nurbs-convert "$@"
fi

echo "✅ Done !"
