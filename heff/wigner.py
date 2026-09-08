"""3j/6j/9j behind one interface.

The backend is ``sympy.physics.wigner`` behind ``functools.lru_cache``. Matrix
construction reuses each symbol within a basis. The small public interface
keeps the backend replaceable.

SymPy is imported inside the kernels so importing ``heff`` does not require it.

Arguments are half-integers as Python floats. They are converted to exact
sympy Rationals via round(2*x)/2 -- passing a float straight to sympy risks a
Float argument where the routine wants a Rational. Basis enumeration produces
values that are exact multiples of 0.5, so this conversion preserves them.
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
    """Drop every cached symbol."""
    for f in (_w3j, _w6j, _w9j):
        f.cache_clear()
