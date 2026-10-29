# C++ WebAssembly adapter

Development preview, 28 October 2026. No JavaScript replacement algorithms.

Install official emsdk manager commit `d7b4130688bfa348392b921ab04015aa15c4f154`, install and activate SDK **6.0.13**. This is an explicit toolchain installation, not a configure-time download. Run from this repository:

```sh
python tools/build_wasm.py --emsdk /path/to/emsdk --engine numerical --source . --out build-wasm
node tests/browser-smoke.mjs build-wasm
```

The script compiles the same native core sources with C++20, exceptions and contraction disabled. Embind exports one JSON-returning function. Manifest records core/adapter source hashes, output hashes and flags. Memory starts at 16 MiB and is capped at 128 MiB; stack 1 MiB; virtual filesystem disabled. The portfolio hosts modules in a dedicated worker, with a 20-second load deadline and 5-second computation deadline; timeout/cancel terminates the worker. This bounds work, not a guarantee against every browser resource failure.

Numerical adapter accepts row-major textual input and converts to the engine's column-major storage; order 1–16, RHS 1–4, finite magnitude <=1e100, 32 KiB wire limit. Exact singularity is rejected; backward error can overflow and causes an explicit diagnostic. Small residual is not a forward-error or condition guarantee. Factors and zero-based sequential pivots retain native semantics.

Symbolic adapter limits ASCII expressions to 256 bytes, variables to 32 characters, literals to 64 digits, trees to 64 nodes/depth 32, direct integer powers -16..16, conservative arithmetic-growth cost 512 and output 64 KiB. Core domain restrictions remain intact, including x^0 requiring x!=0. Browser limits are intentionally narrower than the native API; rejection is not a claim that the native mathematics is invalid.

The smoke suite instantiates the actual generated module in Node. Cross-engine native/WASM comparison is additionally run in the private lab integration suite. No original DLL behavior is exercised. Generated runtimes require their own notices before publication; no WASM files or emsdk distribution are committed.
