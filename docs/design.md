# Design — Scaled QR Evidence

## Architecture

1. **`sqe_types.f90`**: Fixed-bound result and budget types (`sqe_result_t`, `sqe_budget_t`).
2. **`sqe_validation.f90`**: Pre-allocation IEEE finiteness, magnitude range $[2^{-100}, 2^{100}]$, and dimension validation.
3. **`sqe_solver.f90`**: Scaled Euclidean norms (`scaled_vec_norm2`), column equilibration, column-pivoted Householder QR, explicit orthogonal $Q$ accumulation, up to 3 iterative refinement steps, and first-order optimality residual calculation.
4. **`sqe_verifier.f90`**: Independent Frobenius reconstruction ($\|A_s P - Q R\|_F / \|A_s\|_F$) and orthogonality ($\|Q^T Q - I\|_F$) checker.
5. **`test/test_suite.f90`**: Native Fortran test suite linking LAPACK `DGELSD` as an external SVD oracle while keeping the production library and CLI free of LAPACK/BLAS dependencies.
