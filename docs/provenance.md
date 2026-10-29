# Provenance categories

Binary observation: PE x86-64 LAPACK DLL shipped with Octave 11.3.0; SHA-256 5132bafddc40afdedaa6e4ad7e72cce27595f4b395cff9cbc3332f4b19cc0c9d. Static DGETRF2 VA 0x2e066ace0, DLAMCH 0x2e07e4770, DLASWP 0x2e06c2a80. DGETRF2 evidence ID ev_d5506c16009862d2a35427eecb73ebf07fecb1eaef88690f16ac9501a6316a14. These are recorded prior REA findings, not a fresh binary import or source-code evidence. Private raw evidence is retained outside this proposed repository.

Implementation: recursive LU layout/control flow is informed by the observed binary and may be derivative of LAPACK. Scalar BLAS replacements, span checking, zero-based API and solve logic are project implementations. Treat LU conservatively for licensing; assembly translation does not erase upstream rights. LAPACK's published BSD-style terms are a candidate basis, not proof of the selected DLL's complete provenance. Package-to-upstream source mapping and target distribution notices must be retained/verified before release.

Source reference: https://github.com/Reference-LAPACK/lapack/blob/v3.12.0/LICENSE describes that upstream version's BSD-style conditions. Octave itself is GPL software; do not infer that all bundled libraries share its license or that Octave wrapper code can be relicensed. No Octave wrapper source is included here.

Theory: PA=LU and triangular substitution describe the mathematics independently of implementation evidence. Validation against NumPy checks numerical samples; it does not prove exact agreement with the selected DLL. Original dependency closure remains unresolved.

## Release source-license review

The previous conservative gate is resolved for this source-only release by reviewing the entire LU implementation against official permissively licensed DGETRF2 v3.12.0 and retaining its full upstream notice. See THIRD_PARTY_NOTICES.md for comparison scope and license decision. Exact DLL-to-source equivalence and original-DLL differential execution remain unknown; neither is claimed or required to redistribute an independently compiled adaptation under the retained upstream permissions. This does not authorize the original DLL redistribution.
