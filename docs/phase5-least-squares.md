# Phase 5 scope — checked QR least squares

Native-only extension of the existing Householder QR. Preserve LU, solve, QR,
browser adapters, license and Phase4 artifacts. Original authored mathematics;
no new REA findings or original-binary execution.

Contract before implementation: finite real binary64, column-major A/B with
validated leading dimensions; 1 <= n <= m <=256 and 1..16 RHS. Solve each RHS's
Euclidean least-squares problem by existing thin QR, Q^T B, checked back solve.
Inputs are immutable. Underdetermined/minimum-norm problems are unsupported.

Default rejection tolerance is epsilon*max(m,n); reject if any |Rii|/max|R|
is <= tolerance. A caller may supply a finite tolerance in[0,1). This unpivoted
diagonal test is a conservative acceptance policy, not an SVD rank estimate;
it is order/scale sensitive and may reject full-rank ill-conditioned systems.
No solution is returned on rejection; pivoted QR/SVD are future work.

Phase4 review: the browser solver's infinity-norm scaled residual and recorded
backward-error figures remain unchanged. They describe a residual relative to
input/solution scales, not a forward-error guarantee, a condition estimate or
a rank test. Batched RHS reporting does not compute a minimal common structured
perturbation of A. The new least-squares diagnostics use a different, explicitly
documented norm and are not substituted into those historical measurements.

Diagnostics: per-RHS residual2 norm, scaled residual
||Ax-b||2/(||A||F||x||2+||b||2), and normalized stationarity
||A^T(Ax-b)||2/(||A||F||Ax-b||2). Nonzero residual is normal for inconsistent
least squares, so this scaled residual is not a least-squares backward-error
certificate. Stationarity is an optimality diagnostic, not forward accuracy;
it can be dominated by roundoff when the residual is near machine precision.
Optional cond_inf(R) estimate uses checked triangular inverse solves; it is
not cond2(A), a rank proof or a forward-error bound. Nonfinite intermediates
and diagnostic overflow throw overflow_error, including cases a differently
scaled formulation could handle. No normal-equation solve is introduced.

Independent validation: NumPy lstsq (SVD reference), residual/orthogonality
checks, padded/multiple RHS, ill-conditioned/singular/overflow cases. Benchmark
end-to-end API latency against NumPy, with deterministic inputs, warmup, raw
samples, hardware/compiler/flags/versions and explicit timing boundaries.
References: https://netlib.org/lapack/lug/node27.html and
https://numpy.org/doc/2.3/reference/generated/numpy.linalg.lstsq.html . These are
mathematical/reference context, distinct from recorded REA binary evidence.
