# Scaled QR Evidence

Development version `0.1.0` — a bounded dense least-squares library and CLI in Fortran 2008 (`gfortran` / `fpm`) with original scaled column-pivoted Householder QR, up to 3 iterative refinement corrections, independent reconstruction diagnostics, and test-time verification against native LAPACK `DGELSD`.

Built by [Lofrano Tecnologia](https://lofranotecnologia.com.br). [Versão em português](README.pt-BR.md).

## Overview

`scaled-qr-evidence` solves full-column-rank overdetermined linear least-squares problems $\min_X \|AX - B\|_F$ for finite `real64` arrays ($1 \le n \le 32$, $n \le m \le 128$, $1 \le k_{\text{rhs}} \le 4$):
- **Column equilibration & column-pivoted Householder QR**: scales columns by their overflow-safe Euclidean norms, selects pivots by maximum residual norm, and checks $|R_{jj}| > 2^{-k} |R_{11}|$ ($k \in [10, 40]$).
- **Iterative refinement**: performs up to 3 normal-equation refinement steps on the equilibrated triangular factors.
- **Independent diagnostics**: checks $\|A_s P - Q(:,1:n)R\|_F / \|A_s\|_F \le 5 \times 10^{-12}$, $\|Q^T Q - I\|_F \le 5 \times 10^{-12}$, and normalized first-order optimality $\le 5 \times 10^{-11}$.
- **Explicit status taxonomy**: `COMPUTED`, `INVALID_INPUT`, `NUMERIC_UNRESOLVED`, `UNKNOWN_BUDGET`.
- **Zero runtime LAPACK dependency**: LAPACK `DGELSD` is linked only in the test oracle suite.

## Quick Start

```sh
/bin/sh scripts/test-offline.sh
/bin/sh scripts/demo-offline.sh --assert
/bin/sh scripts/demo-offline.sh --assert --json
```

Example output:

```text
PASS status=COMPUTED rel_recon_err= 1.50083E-16 orth_err= 4.67892E-16 normal_eq_rel_residual= 1.24019E-16
```

## Limitations and Scope Exclusions

- Rank screening via $2^{-k}$ is a numerical policy, not an exact algebraic rank certificate.
- Does not solve underdetermined ($m < n$) or minimum-norm rank-deficient systems.

## License

MIT — Copyright (c) 2026 Lofrano Tecnologia. See [LICENSE](LICENSE).
