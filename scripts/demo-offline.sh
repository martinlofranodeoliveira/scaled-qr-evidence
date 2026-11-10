#!/bin/sh
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_ROOT}"
mkdir -p dist

if [ ! -x "dist/scaled-qr-evidence" ]; then
    BUILD_TMP="$(mktemp -d /tmp/sqe-demo.XXXXXX)"
    trap 'rm -rf "${BUILD_TMP}"' EXIT
    gfortran -O2 -ffloat-store -J"${BUILD_TMP}" \
        src/sqe_types.f90 \
        src/sqe_validation.f90 \
        src/sqe_solver.f90 \
        src/sqe_verifier.f90 \
        src/scaled_qr_evidence.f90 \
        app/main.f90 \
        -o dist/scaled-qr-evidence
fi

if [ "$#" -eq 0 ]; then
    ./dist/scaled-qr-evidence demo --assert
elif [ "$1" = "--assert" ] || [ "$1" = "--json" ]; then
    ./dist/scaled-qr-evidence demo "$@"
else
    ./dist/scaled-qr-evidence "$@"
fi
