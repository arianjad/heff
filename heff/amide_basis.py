"""Basis enumeration for C2v M-NH2 effective-Hamiltonian models.

This module defines only the coupled-state bookkeeping.  It deliberately
contains no isotope constants, molecular parameters, or operators.  The
equivalent-proton exchange filter is ported from
``C2V-Molecules/atm_core/physics.py::generate_states_and_blocks``; this module
records that source convention and does not make an independent molecular
symmetry claim.
"""

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np


AMIDE_KET = np.dtype([
    ("N", "f8"),
    ("K", "f8"),
    ("J", "f8"),
    ("F_N", "f8"),
    ("I_T", "f8"),
    ("F_core", "f8"),
    ("F", "f8"),
    ("mF", "f8"),
])


def _half_integer(name, value):
    if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite nonnegative integer or half-integer")
    doubled = 2.0 * float(value)
    if value < 0 or doubled != round(doubled):
        raise ValueError(f"{name} must be a finite nonnegative integer or half-integer")
    return float(value)


def _integer(name, value, *, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    value = int(value)
    if nonnegative and value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _coupled_values(first, second):
    """Allowed totals ``|first-second|, ..., first+second`` in unit steps."""
    two_first, two_second = int(round(2.0 * first)), int(round(2.0 * second))
    return tuple(
        two_total / 2.0
        for two_total in range(abs(two_first - two_second), two_first + two_second + 1, 2)
    )


@dataclass(frozen=True)
class AmideBasisSpec:
    """Inputs for ``|((((N S) J I_N) F_N I_T) F_core I_M) F mF>``.

    ``I_T`` is the total spin of the equivalent proton pair, while ``I_M`` is
    one metal nuclear spin.  ``N_range`` is an inclusive pair of nonnegative
    integer endpoints.  ``K_values`` are absolute requested magnitudes; input
    signs are normalized and each nonzero value produces both ``+K`` and ``-K``.
    """

    S: float
    I_N: float
    i_H: float
    I_M: float
    N_range: tuple
    K_values: tuple
    vibronic_sign: int
    M: str = "all"
    frame: str = "lab"

    def __post_init__(self):
        for name in ("S", "I_N", "i_H", "I_M"):
            object.__setattr__(self, name, _half_integer(name, getattr(self, name)))

        try:
            n_range = tuple(self.N_range)
        except TypeError as exc:
            raise ValueError("N_range must be a pair of nonnegative integer endpoints") from exc
        if len(n_range) != 2:
            raise ValueError("N_range must be a pair of nonnegative integer endpoints")
        n_range = tuple(_integer("N_range", value, nonnegative=True) for value in n_range)
        if n_range[0] > n_range[1]:
            raise ValueError("N_range must satisfy N_range[0] <= N_range[1]")
        object.__setattr__(self, "N_range", n_range)

        try:
            requested_K = tuple(self.K_values)
        except TypeError as exc:
            raise ValueError("K_values must contain integer magnitudes") from exc
        if not requested_K:
            raise ValueError("K_values must not be empty")
        K_values = tuple(sorted({abs(_integer("K_values", value)) for value in requested_K}))
        object.__setattr__(self, "K_values", K_values)

        if isinstance(self.vibronic_sign, bool) or self.vibronic_sign not in (-1, 1):
            raise ValueError("vibronic_sign must be +1 or -1")
        object.__setattr__(self, "vibronic_sign", int(self.vibronic_sign))
        if self.M not in {"all", "blocks"}:
            raise ValueError("M must be 'all' or 'blocks'")
        if self.frame != "lab":
            raise ValueError("frame must be 'lab'")


def enumerate_amide_kets(spec: AmideBasisSpec):
    """Return the full named ``float64`` basis for an amide specification.

    ``M='all'`` and ``M='blocks'`` intentionally return the same full basis;
    callers can apply the existing named-``mF`` blocking layer after
    enumeration.  The metal spin is appended only after the legacy core chain:
    ``N + S = J``, ``J + I_N = F_N``, ``F_N + I_T = F_core``, and
    ``F_core + I_M = F``.
    """
    if not isinstance(spec, AmideBasisSpec):
        raise TypeError("spec must be an AmideBasisSpec")

    K_values = tuple(
        K
        for magnitude in spec.K_values
        for K in ((0,) if magnitude == 0 else (-magnitude, magnitude))
    )
    rows = []
    required_exchange_sign = (-1) ** int(round(2.0 * spec.i_H))
    for K in K_values:
        for N in range(spec.N_range[0], spec.N_range[1] + 1):
            if N < abs(K):
                continue
            for I_T in _coupled_values(spec.i_H, spec.i_H):
                nuclear_exchange_sign = (-1) ** int(round(2.0 * spec.i_H - I_T))
                if (spec.vibronic_sign * (-1) ** K * nuclear_exchange_sign
                        != required_exchange_sign):
                    continue
                for J in _coupled_values(N, spec.S):
                    for F_N in _coupled_values(J, spec.I_N):
                        for F_core in _coupled_values(F_N, I_T):
                            for F in _coupled_values(F_core, spec.I_M):
                                for two_mF in range(-int(round(2.0 * F)),
                                                    int(round(2.0 * F)) + 1, 2):
                                    rows.append((float(N), float(K), J, F_N, I_T,
                                                 F_core, F, two_mF / 2.0))
    return np.array(rows, dtype=AMIDE_KET)


__all__ = ["AMIDE_KET", "AmideBasisSpec", "enumerate_amide_kets"]
