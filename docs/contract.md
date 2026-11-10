# Contract — Scaled QR Evidence (v0.1.0)

## Input Domain

- Dimensions: $1 \le n \le 32$, $n \le m \le 128$, $1 \le k_{\text{rhs}} \le 4$.
- Entries: finite IEEE `real64` values; every nonzero magnitude $|x|$ must satisfy $2^{-100} \le |x| \le 2^{100}$.
- Rank-screening parameter: integer $k \in [10, 40]$, defining pivot threshold $\tau = 2^{-k}$.
- Resource budgets: `max_updates = 2,000,000`, `max_scratch_bytes = 1,048,576`, `max_output_bytes = 65,536`.

## Outcomes

- `COMPUTED` (`0`): All pivots satisfy $|R_{jj}| > 2^{-k} |R_{11}|$, solution components are finite, and normalized first-order residual satisfies $\eta \le 5 \times 10^{-11}$.
- `INVALID_INPUT` (`1`): Dimension out of range, `NaN`/`Inf` entry, or nonzero magnitude outside $[2^{-100}, 2^{100}]$.
- `NUMERIC_UNRESOLVED` (`2`): Zero column norm, screened small pivot, or residual diagnostic above threshold.
- `UNKNOWN_BUDGET` (`3`): Arithmetic update or scratch byte budget exceeded.
