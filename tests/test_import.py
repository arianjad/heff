"""Gate A8 -- import cost and dependency hygiene.

Uniquely catches: a heavy or plotting dependency creeping into package import
time. This is the exact regression C2V-Molecules has today (6.5 s import from
dead pandas/sklearn/matplotlib/sympy imports at atm_core/config.py:8-20, spec
S3.7 A8). FAIL is reachable in one line: add `import sympy` to the top of
heff/__init__.py and this test fails on the `forbidden` assertion.
"""
import subprocess
import sys

FORBIDDEN = ("matplotlib", "pandas", "sklearn", "sympy", "torch")


def _probe():
    code = (
        "import time, sys; t = time.perf_counter(); import heff; "
        "dt = time.perf_counter() - t; "
        "bad = [m for m in ('matplotlib','pandas','sklearn','sympy','torch') if m in sys.modules]; "
        "print(dt); print(','.join(bad))"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    lines = out.stdout.splitlines()
    bad_line = lines[1] if len(lines) > 1 else ""
    return float(lines[0]), [m for m in bad_line.split(",") if m]


def test_import_is_under_one_second():
    dt, _ = _probe()
    assert dt < 1.0, f"import heff took {dt:.2f} s"


def test_import_pulls_no_heavy_dependency():
    _, bad = _probe()
    assert bad == [], f"import heff pulled {bad}; every one of {FORBIDDEN} must be lazy"
