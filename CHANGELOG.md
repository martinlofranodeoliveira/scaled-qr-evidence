# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - Unreleased

### Added
- Original Fortran 2008 scaled column-pivoted Householder QR solver (`sqe_solver.f90`).
- Iterative refinement and first-order residual diagnostics.
- Independent QR reconstruction and orthogonality verifier (`sqe_verifier.f90`).
- Native LAPACK `DGELSD` oracle test suite (`test/test_suite.f90`).
