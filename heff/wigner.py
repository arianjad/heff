"""3j/6j/9j behind one interface.

Backend: sympy.physics.wigner under functools.lru_cache -- exact, and already
what both source repos use (Molecule-Structure Source Code/matrix_elements.py:18-23;
C2V-Molecules atm_core/physics.py:12-27, which records a 99.9 % hit rate). The
matrix-first architecture pays the Wigner cost once per basis rather than once
per parameter set, which is why a fast-but-fragile dependency is not needed:
`wigners` (Rust) ships wheels but has no 6j/9j, and `py3nj` has 6j/9j but ships
no wheels at all (spec S3.8). This module is an interface so a fast 3j path can
be swapped in behind gate A6 later.

sympy is imported INSIDE the kernels: `import heff` must stay under 1 s and must
not pull sympy (gate A8).

Arguments are half-integers as Python floats. They are converted to exact
sympy Rationals via round(2*x)/2 -- passing a float straight to sympy risks a
Float argument where the routine wants a Rational. Gate A7 guarantees every
value coming out of the enumerator is an exact multiple of 0.5, so the rounding
never loses information.
"""
from functools import lru_cache

__all__ = ["w3j", "w6j", "w9j", "cache_info", "clear_cache"]


def _r(x):
    from sympy import Rational

    return Rational(int(round(2.0 * float(x))), 2)


@lru_cache(maxsize=None)
def _w3j(j1, j2, j3, m1, m2, m3):
    from sympy.physics.wigner import wigner_3j

    return float(wigner_3j(_r(j1), _r(j2), _r(j3), _r(m1), _r(m2), _r(m3)))


@lru_cache(maxsize=None)
def _w6j(j1, j2, j3, j4, j5, j6):
    from sympy.physics.wigner import wigner_6j

    return float(wigner_6j(_r(j1), _r(j2), _r(j3), _r(j4), _r(j5), _r(j6)))


@lru_cache(maxsize=None)
def _w9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    from sympy.physics.wigner import wigner_9j

    args = [_r(x) for x in (j1, j2, j3, j4, j5, j6, j7, j8, j9)]
    return float(wigner_9j(*args, prec=64))


def _key(x):
    return round(2.0 * float(x)) / 2.0


def w3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (j1 j2 j3 / m1 m2 m3) as a float."""
    return _w3j(*(_key(x) for x in (j1, j2, j3, m1, m2, m3)))


def w6j(j1, j2, j3, j4, j5, j6):
    """Wigner 6j symbol {j1 j2 j3 / j4 j5 j6} as a float."""
    return _w6j(*(_key(x) for x in (j1, j2, j3, j4, j5, j6)))


def w9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    """Wigner 9j symbol as a float."""
    return _w9j(*(_key(x) for x in (j1, j2, j3, j4, j5, j6, j7, j8, j9)))


def cache_info():
    """Cache statistics for all three kernels."""
    return {"w3j": _w3j.cache_info(), "w6j": _w6j.cache_info(), "w9j": _w9j.cache_info()}


def clear_cache():
    """Drop every cached symbol (used by gate A6 to measure hit rate)."""
    for f in (_w3j, _w6j, _w9j):
        f.cache_clear()
