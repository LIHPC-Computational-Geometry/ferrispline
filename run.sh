#!/bin/bash

set -e

echo "🟢 Step 1 : Check the virtual environnement..."

if [ -z "$VIRTUAL_ENV" ]; then
    echo "   -> Any environnement detected. .venv activation..."
    source .venv/bin/activate
else
    echo "   -> The virtual environnement is already active ($VIRTUAL_ENV)."
fi

echo "⚙️ Step 2 : Compilation of the Rust's librairy (nurbslib)..."
cd nurbslib
maturin develop
cd ..

echo "🚀 Step 3 : Launch of nurbs-convert..."
cd sandbox_python
nurbs-convert "$@"

echo "✅ Done !"
