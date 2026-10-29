"""Independent numerical reference, NOT differential validation of Octave's DLL."""
import argparse
import contextlib
import io
import json
import os
import pathlib
import platform
import hashlib
import sys
import datetime
import statistics
import subprocess
import time

# Set before importing NumPy. No package is downloaded or installed.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np

# Assertions are the validation oracle: fail closed if Python disables them.
if not __debug__:
    raise RuntimeError("Numerical validation requires Python without -O or PYTHONOPTIMIZE")

parser = argparse.ArgumentParser()
parser.add_argument("--probe", required=True)
parser.add_argument("--output")
parser.add_argument("--benchmark", action="store_true")
args = parser.parse_args()
rng = np.random.default_rng(20261025)
cases = []
for n in range(21):
    for _ in range(5):
        cases.append((rng.uniform(-1, 1, size=(n, n)), 1, True))
for a in (np.zeros((4, 4)), np.array([[1., 2.], [2., 4.]]),
          np.eye(4)*1e-200, np.eye(4)*1e200):
    cases.append((a, 1, True))
if args.benchmark:
    for n in (32, 64, 128, 256):
        cases.append((rng.uniform(-1, 1, size=(n, n)), 7, False))
protocol = "".join(
    f"{len(a)} {reps} {int(packed)}\n" +
    " ".join(format(float(x), ".17g") for x in a.ravel(order="F")) + "\n"
    for a, reps, packed in cases)
run = subprocess.run([args.probe], input=protocol, text=True, capture_output=True,
                     timeout=120, check=True)
results = [json.loads(line) for line in run.stdout.splitlines()]
assert len(results) == len(cases), "Missing probe results"
max_log_error = 0.
max_scaled_residual = 0.
benchmarks = []
for (a, reps, packed_output), result in zip(cases, results):
    sign, logarithm = np.linalg.slogdet(a)
    assert int(sign) == result["sign"], (a.shape, sign, result)
    if sign:
        difference = abs(float(logarithm)-result["logabsdet"])
        max_log_error = max(max_log_error, difference)
        assert difference <= 1e-9*max(1., abs(float(logarithm))), (a.shape, difference)
    else:
        assert result["info"] > 0
    n = len(a)
    if packed_output and n:
        lu = np.array(result["packed"]).reshape((n,n), order="F")
        permuted = a.copy()
        for row, pivot in enumerate(result["pivots"]):
            assert row <= pivot < n
            permuted[[row,pivot],:] = permuted[[pivot,row],:]
        product = (np.tril(lu,-1)+np.eye(n)) @ np.triu(lu)
        norm = max(float(np.linalg.norm(a, ord=np.inf)), np.finfo(float).tiny)
        residual = float(np.linalg.norm(permuted-product, ord=np.inf))/norm
        max_scaled_residual = max(max_scaled_residual,residual)
        assert residual <= 128*np.finfo(float).eps*max(1,n), (a.shape,residual)
    if not packed_output:
        # API-level timing: native probe includes copying, checks, LU and logdet;
        # NumPy includes Python dispatch, its internal copy, LU and logdet.
        np.linalg.slogdet(a)  # Warm up; no startup/import time in measurements.
        samples = []
        for _ in range(reps):
            start = time.perf_counter_ns()
            np.linalg.slogdet(a)
            samples.append(time.perf_counter_ns()-start)
        benchmarks.append({"n":n, "repetitions":reps,
                           "cpp_median_ns":result["median_ns"],
                           "numpy_median_ns":statistics.median(samples),
                           "cpp_samples_ns":result["samples_ns"],
                           "numpy_samples_ns":samples})
config = io.StringIO()
with contextlib.redirect_stdout(config):
    np.show_config()
report = {"recorded_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "seed":20261025, "python_version":sys.version,
          "platform":platform.platform(), "machine":platform.machine(),
          "processor":platform.processor(), "logical_cpu_count":os.cpu_count(),
          "probe_path":str(pathlib.Path(args.probe).resolve()),
          "probe_sha256":hashlib.sha256(pathlib.Path(args.probe).read_bytes()).hexdigest(),
          "script_sha256":hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
          "requested_threads":{"OPENBLAS_NUM_THREADS":os.environ["OPENBLAS_NUM_THREADS"],
                               "OMP_NUM_THREADS":os.environ["OMP_NUM_THREADS"]},
          "reference":"NumPy.linalg.slogdet, independent of the Octave target",
          "numpy_version":np.__version__, "numpy_path":np.__file__,
          "cases_checked":len(cases), "max_absolute_logdet_error":max_log_error,
          "max_relative_infinity_residual":max_scaled_residual,
          "benchmarks":benchmarks, "numpy_build_config":config.getvalue(),
          "limitations":["Not original-Octave differential tests or packed-factor equality",
                         "API-level timings have different language/validation overheads",
                         "No claims of bitwise equivalence, speed parity or universal stability"]}
if args.output:
    pathlib.Path(args.output).write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps({k:v for k,v in report.items() if k not in ("numpy_build_config","numpy_path")},indent=2))
