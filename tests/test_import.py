"""Package imports must not initialize optional or heavy dependencies."""
import subprocess
import sys

FORBIDDEN = ("matplotlib", "pandas", "sklearn", "sympy", "torch")


def test_import_pulls_no_heavy_dependency():
    code = (
        "import sys; import heff; "
        f"bad = [m for m in {FORBIDDEN!r} if m in sys.modules]; "
        "assert not bad, f'Package import initialized {bad}'"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
