"""Gate A6 -- Wigner backend vs uncached sympy on a random half-integer sample.

Uniquely catches: backend precision loss and half-integer argument-convention
mismatch (a fast 3j library that takes 2j instead of j, or that silently
truncates a half-integer to an int). Sampling INCLUDES half-integers on
purpose -- an int-only sample passes with a backend that cannot do half-integer
arguments at all. FAIL is reachable: change `_r` in heff/wigner.py to
`Rational(int(round(float(x))), 1)` and every half-integer case disagrees.
"""
import random

import pytest

from heff import wigner


def _half(rng, lo, hi):
    return rng.randrange(2 * lo, 2 * hi + 1) / 2.0


def test_w3j_matches_uncached_sympy():
    from sympy import Rational
    from sympy.physics.wigner import wigner_3j

    rng = random.Random(20260905)
    checked = 0
    for _ in range(200):
        j1, j2, j3 = (_half(rng, 0, 5) for _ in range(3))
        m1, m2 = _half(rng, -5, 5), _half(rng, -5, 5)
        m3 = -(m1 + m2)
        R = lambda x: Rational(int(round(2 * x)), 2)
        try:
            ref = float(wigner_3j(R(j1), R(j2), R(j3), R(m1), R(m2), R(m3)))
        except ValueError:
            continue
        assert wigner.w3j(j1, j2, j3, m1, m2, m3) == pytest.approx(ref, abs=1e-13)
        checked += 1
    assert checked > 50, f"sample was degenerate: only {checked} legal 3j drawn"


def test_w6j_matches_uncached_sympy():
    from sympy import Rational
    from sympy.physics.wigner import wigner_6j

    rng = random.Random(1729)
    checked = 0
    for _ in range(120):
        js = [_half(rng, 0, 4) for _ in range(6)]
        R = lambda x: Rational(int(round(2 * x)), 2)
        try:
            ref = float(wigner_6j(*[R(x) for x in js]))
        except ValueError:
            continue
        assert wigner.w6j(*js) == pytest.approx(ref, abs=1e-13)
        checked += 1
    assert checked > 10, f"sample was degenerate: only {checked} legal 6j drawn"


def test_w9j_matches_uncached_sympy():
    from sympy import Rational
    from sympy.physics.wigner import wigner_9j

    R = lambda x: Rational(int(round(2 * x)), 2)
    cases = [
        (1, 1, 2, 1, 1, 2, 2, 2, 2),
        (0.5, 0.5, 1, 0.5, 0.5, 1, 1, 1, 2),
        (1.5, 0.5, 1, 1, 1, 2, 0.5, 0.5, 1),
    ]
    for js in cases:
        ref = float(wigner_9j(*[R(x) for x in js], prec=64))
        assert wigner.w9j(*js) == pytest.approx(ref, abs=1e-12)


def test_cache_is_actually_used():
    wigner.clear_cache()
    wigner.w3j(1, 1, 2, 0, 0, 0)
    wigner.w3j(1, 1, 2, 0, 0, 0)
    info = wigner.cache_info()["w3j"]
    assert info.hits >= 1 and info.misses == 1


def test_half_integer_arguments_are_exact_not_rounded():
    """A 3j with half-integer j must not equal the same call with j rounded."""
    a = wigner.w3j(0.5, 0.5, 1, 0.5, -0.5, 0)
    b = wigner.w3j(1.0, 1.0, 1, 0.5, -0.5, 0)
    assert a != b
