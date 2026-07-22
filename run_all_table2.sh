#!/usr/bin/env bash

# Stop immediately if a command fails.
set -euo pipefail

# Locate the repository root, regardless of the current directory.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=================================================="
echo "Running all Table 2 experiments"
echo "Repository: ${ROOT_DIR}"
echo "=================================================="

echo
echo "[1/5] Running CIHNP-CSURF..."
sage "${ROOT_DIR}/CIHNP_CSURF.sage"

echo
echo "[2/5] Running MIHNP..."
sage "${ROOT_DIR}/MIHNP.sage"

echo
echo "[3/5] Running ECHNP..."
sage "${ROOT_DIR}/ECHNP.sage"

echo
echo "[4/5] Running LCG..."
sage "${ROOT_DIR}/LCG.sage"

echo
echo "[5/5] Running LIPH-POKE..."
sage "${ROOT_DIR}/LIPH_POKE.sage"

echo
echo "=================================================="
echo "All Table 2 experiments completed."
echo "=================================================="
