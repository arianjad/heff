"""Shared test helpers. One home, so a fix lands everywhere at once.

`_elem` and `_find` are used by every element test; `_gamma` and `_kappa` are the
closed-form projection factors that [HAM] V2 (Stark) and V5 (g-factors) compare
against, written out once here instead of once per test module.
"""
import numpy as np

from heff.terms import REGISTRY


def _elem(name, basis, ctx, i, j):
    """<basis[i]| <the registered term `name`> |basis[j]>, parameter-free."""
    return REGISTRY[name].fn(basis[i], basis[j], ctx)


def _find(basis, J, Om, F, mF):
    """Index of the one ket with these four quantum numbers."""
    hit = np.flatnonzero((basis["J"] == J) & (basis["Om"] == Om)
                         & (basis["F"] == F) & (basis["mF"] == mF))
    assert len(hit) == 1, f"no unique ket for J={J} Om={Om} F={F} mF={mF}"
    return int(hit[0])


def _gamma(J, F, *, I):
    """[J(J+1)+F(F+1)-I(I+1)] / [2F(F+1)J(J+1)] -- Leanhardt Eqs. 20-21.

    I is keyword-only with no default: a nuclear spin is a physical quantity,
    and the whole point of these closed forms is their I dependence.
    """
    return (J * (J + 1) + F * (F + 1) - I * (I + 1)) / (2 * F * (F + 1) * J * (J + 1))


def _kappa(J, F, *, I):
    """[F(F+1)-J(J+1)+I(I+1)] / [2F(F+1)] -- the I-on-F projection factor.

    I is keyword-only with no default, as for _gamma.
    """
    return (F * (F + 1) - J * (J + 1) + I * (I + 1)) / (2 * F * (F + 1))
