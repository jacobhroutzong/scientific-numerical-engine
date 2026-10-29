# Third-party notices and source provenance

## Project source
Original project additions use the root BSD-3-Clause license. The recursive LU control flow in src/lu.cpp is conservatively treated as a LAPACK adaptation and retains the complete upstream LAPACK3.12.0 notice in LICENSES/LAPACK-3.12.0.txt, including its additional disclaimer. DGETRF2 authorship identifies the Universities of Tennessee, California Berkeley and Colorado Denver, and NAG Ltd. No endorsement is implied.

Review compared src/lu.cpp to official Reference-LAPACK v3.12.0 SRC/dgetrf2.f: recursive split, panel permutations, triangular/update calls, pivot offsets, safe-minimum base case and INFO continuation. The recovered flow is covered by that permissively licensed upstream implementation; checked spans, zero-based pivots, scalar kernels, solve and QR are authored adaptations/additions. The exact compiled DLL is not asserted byte-identical to upstream source and is not distributed. This review does not depend on redistributing that DLL or Octave GPL wrapper code. No wrapper or other GPL numerical implementation is present in this source release.

Source-license references, separate from REA evidence:
- https://github.com/Reference-LAPACK/lapack/blob/v3.12.0/LICENSE
- https://github.com/Reference-LAPACK/lapack/blob/v3.12.0/SRC/dgetrf2.f

## External testing and toolchain dependencies

GoogleTest1.17.0 is fetched separately by CI at commit52eb8108c5bdec04579160ae17225d66034bd723. Its BSD-3-Clause license remains with that external checkout; it is not vendored in this source release. NumPy2.3.5 (numerical references) and SymPy1.14.0/mpmath1.3.0 (symbolic references) are separately acquired test tools, not engine runtime dependencies or bundled source. Preserve their notices if you distribute those dependencies yourself.

WebAssembly build uses Emscripten6.0.13. Generated binaries include permissively licensed runtime portions. The SDK's Emscripten LICENSE (MIT/NCSA and Node notices), musl COPYRIGHT, LLVM libc++/libc++abi/compiler-rt/libunwind LICENSE.TXT, and allocator notices are copied without truncation into the portfolio's runtime notice directory. No optional ports or original research DLLs are linked. Before distributing a different build, recheck its dependency/notice inventory.
