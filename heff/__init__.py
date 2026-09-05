"""heff -- effective Hamiltonians for molecules, term-matrix first.

Import stays light on purpose (gate A8): numpy only. sympy is imported inside
heff.wigner's kernels, scipy inside heff.track's assignment call, matplotlib
never (plotting lives in notebooks, spec S3.8).
"""

__version__ = "0.1.0"

from . import wigner  # noqa: F401

__all__ = ["wigner"]
