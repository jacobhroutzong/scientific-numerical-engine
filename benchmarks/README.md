# Phase 5 reproducible measurements

`phase5-windows-msvc.json` records seed20261104, all50 input-case hashes,
dimensions/RHS counts, policies, residual/reference comparisons, CPU/compiler/
flags, NumPy/BLAS versions, executable hash and seven raw timing samples.
Generate inputs with tests/least_squares_reference.py, not a private dataset.

```sh
python tests/least_squares_reference.py --probe build/Release/least_squares_probe.exe --benchmark --compiler "MSVC 19.44 x64 Release" --flags "/O2 /Ob2 /DNDEBUG /MD /W4 /fp:strict; C++20" --cpu "your actual CPU" --output benchmarks/reproduction.json
```

Recorded AMD Ryzen9 8945HS/MSVC19.44/NumPy2.3.5/OpenBLAS0.3.30 example:

| Case | Native median us | NumPy median us |
| --- | ---: | ---: |
| tall11x4,1 RHS | 3.781 | 21.903 |
| tall35x16,3 RHS | 107.856 | 85.684 |
| square32,3 RHS | 239.350 | 256.350 |
| square64,3 RHS | 1293.550 | 831.850 |

Nine accepted cases are timed; singular and near-singular rejection cases are
recorded but not assigned a misleading solve-speed comparison. Native timing
excludes subprocess/serialization and includes validation, QR, explicit Q,
projection, solve and diagnostics. Optional R condition computation is off.
NumPy times its lstsq API (SVD, different outputs/contract). Three warmup calls,
seven samples,32 calls/sample at n<=16 and4 calls/sample otherwise. One BLAS
thread was requested through environment variables before importing NumPy;
the actual thread count was not independently measured. No performance
superiority, portable speed guarantee or isolated-kernel equivalence is claimed.

Validation:48 accepted systems,2 policy rejections. Maximum sampled solution
difference relative to the SVD reference was8.984806e-8 on Hilbert8 with
cond2(A)≈1.525758e10; Hilbert4 difference1.953924e-13. Tolerances are recorded
per case and conditioning-dependent, not a guaranteed forward-error bound.
Existing LU timings/benchmark artifacts are preserved separately.
