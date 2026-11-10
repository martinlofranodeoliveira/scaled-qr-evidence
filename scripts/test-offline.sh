#!/bin/sh
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

WORK_TMP="$(mktemp -d /tmp/sqe-test.XXXXXX)"
trap 'rm -rf "${WORK_TMP}"' EXIT

mkdir -p "${WORK_TMP}/pkg" "${WORK_TMP}/installed/lib" "${WORK_TMP}/installed/include" "${WORK_TMP}/consumer"
cp -R "${PROJECT_ROOT}/fpm.toml" "${PROJECT_ROOT}/src" "${PROJECT_ROOT}/app" "${PROJECT_ROOT}/test" "${WORK_TMP}/pkg/"

cd "${WORK_TMP}/pkg"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
fpm test --profile debug

# Compile static library and module files into an isolated install prefix (no LAPACK!)
mkdir -p "${WORK_TMP}/obj"
gfortran -O2 -ffloat-store -J"${WORK_TMP}/installed/include" -c \
    src/sqe_types.f90 \
    src/sqe_validation.f90 \
    src/sqe_solver.f90 \
    src/sqe_verifier.f90 \
    src/scaled_qr_evidence.f90
ar rcs "${WORK_TMP}/installed/lib/libscaled_qr_evidence.a" \
    sqe_types.o sqe_validation.o sqe_solver.o sqe_verifier.o scaled_qr_evidence.o

# Build distributed production CLI without LAPACK
mkdir -p "${PROJECT_ROOT}/dist"
gfortran -O2 -ffloat-store -I"${WORK_TMP}/installed/include" \
    app/main.f90 "${WORK_TMP}/installed/lib/libscaled_qr_evidence.a" \
    -o "${PROJECT_ROOT}/dist/scaled-qr-evidence"

# Verify production binary does NOT link liblapack or libblas
if ldd "${PROJECT_ROOT}/dist/scaled-qr-evidence" | grep -E "lapack|blas"; then
    echo "ERROR: production binary accidentally linked LAPACK/BLAS!" >&2
    exit 1
fi

# Compile and run an outside-checkout Fortran consumer against the installed library
cat > "${WORK_TMP}/consumer/consumer.f90" <<'EOF'
program consumer
  use, intrinsic :: iso_fortran_env, only: real64
  use scaled_qr_evidence
  implicit none
  real(real64) :: a(4, 2), b(4, 1)
  type(sqe_budget_t) :: budget
  type(sqe_result_t) :: res

  a(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
  a(:, 2) = [1.0_real64, 1.0_real64, 2.0_real64, 1.0_real64]
  b(:, 1) = 2.0_real64 * a(:, 1) - 3.0_real64 * a(:, 2)

  call solve_scaled_qr(a, b, 4, 2, 1, 25, budget, res)
  if (res%status /= STATUS_COMPUTED) stop 1
  if (abs(res%x(1, 1) - 2.0_real64) > 1.0e-12_real64) stop 1
  if (abs(res%x(2, 1) + 3.0_real64) > 1.0e-12_real64) stop 1
end program consumer
EOF

gfortran -O2 -ffloat-store -I"${WORK_TMP}/installed/include" \
    "${WORK_TMP}/consumer/consumer.f90" "${WORK_TMP}/installed/lib/libscaled_qr_evidence.a" \
    -o "${WORK_TMP}/consumer/consumer_bin"
"${WORK_TMP}/consumer/consumer_bin"
