# Numerical Computing Engine

## Phase 5 checked least squares

Native-only `sce::least_squares` extends the existing thin Householder QR to
full-column-rank square/overdetermined systems with multiple RHS. It returns
explicit fit diagnostics and an optional estimate of cond_inf(R), with a
documented unpivoted rank-rejection policy and checked arithmetic. See
[contract and limitations](docs/phase5-least-squares.md) and
[reproducible benchmarks](benchmarks/README.md). The solver rejects nonzero
subnormal active inputs because the QR calculation can lose significant
precision at that scale. The native and reference suite has 29 entries
(26 GoogleTest and 3 NumPy reference tests).

Configure optional references with `-DSCE_REFERENCE_PYTHON=/path/to/python`
using an environment containing NumPy2.3.5. New probe/reference commands are
documented above; no configure-time downloads are added. CI runs both native
platforms with the independent references and preserves the existing WASM smoke.

The following sections retain the Phase4 baseline and its dated test counts;
the Phase 5 counts above describe the current native API.

C++20 recursive LU with partial pivoting, an independent checked multiple-RHS triangular solver, and a separate independent Householder QR library. Existing LU/solve source is preserved from the validated Scientific Computing workspace. No BLAS/LAPACK runtime dependency; scalar kernels favor clarity over throughput.

## Build and test

CMake 3.24+, C++20 compiler and external GoogleTest 1.17.0 (commit 52eb8108c5bdec04579160ae17225d66034bd723):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSCE_GTEST_SOURCE=/path/to/googletest
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Windows: Visual Studio 2022 x64, use `-A x64` on configure. Library-only builds may set `BUILD_TESTING=OFF`. Configure does not download dependencies. Optional NumPy 2.3.5 reference tests require an external Python and `SCE_REFERENCE_PYTHON`. Native CI runs 21 GoogleTest entries (19 LU/solve + 2 QR); enabling both Python references yields 23 CTest entries. The preserved original workspace separately passes its 21 entries (19 native + 2 Python).

## API and evidence

`include/sce/lu.hpp`: column-major rectangular factorization, packed unit-lower/upper factors, zero-based sequential pivot swaps and one-based first zero-diagonal INFO. `include/sce/solve.hpp`: checked square solves, immutable factors/RHS, dense column-major result. Invalid shape/storage/pivots, singular factors and nonfinite solve arithmetic have explicit exceptions. Finite LU input does not guarantee finite intermediates. Small residual is not forward accuracy for an ill-conditioned problem.

`include/sce/qr.hpp`, target `sce_qr`: independent Householder QR, thin Q/R, 1 <= columns <= rows <= 256. No rank estimate, least-squares solver or original-binary equivalence claim. QR tests check reconstruction and orthogonality, including dependent columns, zero matrices, scaling and padded input. Arithmetic overflow is rejected; extreme cases are not guaranteed to complete even if an equivalent scaled formulation could do so.

REA/Ghidra observed recursive DGETRF2 in the LAPACK DLL distributed with GNU Octave 11.3.0. Its control flow informed LU; scalar kernels, checked boundaries and triangular solving are new C++. This is not a whole-Octave reconstruction. Original-DLL differential validation has not run. No original DLL or decompiled excerpts are included. See [provenance](docs/provenance.md), [browser build](docs/browser.md) and [release review](PUBLICATION.md).

## License and repository

[scientific-numerical-engine](https://github.com/jacobhroutzong/scientific-numerical-engine) — BSD-3-Clause with the retained full LAPACK upstream notice. See [third-party notices](THIRD_PARTY_NOTICES.md). Original research binaries and raw decompilation are excluded.
