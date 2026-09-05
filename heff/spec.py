"""State spec, ket dtype, case (c) enumerator, and the blocking layer.

The spec is data; the enumerator is the only code that knows the coupling case
(spec S3.1). Kets are a structured numpy array with named float64 fields -- it
keeps C2V-Molecules' float64 discipline (their (n,7) basis MUST be float64 or
numpy >= 2 NEP-50 leaks ~1e-7 into every matrix element, atm_core/physics.py:161-169)
without C2V's positional-7-tuple opacity.

v1 basis (docs/thf-plus-x3delta1-effective-hamiltonian.md S1.1): Hund's case (c),
|J, Omega = +-1, F, m_F>, F = J + I with I(19F) = 1/2, J = 1..4. Lambda and Sigma
never appear; the electronic structure enters only through the effective
constants A_par, omega_ef, d_mf, G_par, E_eff.
"""
from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

KET_C = np.dtype([("J", "f8"), ("Om", "f8"), ("F", "f8"), ("mF", "f8")])

_M_MODES = ("none", "blocks", "all")


@dataclass(frozen=True)
class ElecState:
    """One electronic state. No defaults for Omega, S, Lambda (spec S3.2)."""
    label: str
    Omega: float
    S: float
    Lam: float
    T0: float = 0.0


@dataclass(frozen=True)
class StateSpec:
    case: str
    electronic: tuple
    I: float
    J_range: tuple
    v: int = 0
    M: str = "blocks"
    frame: str = "rotating"

    def __post_init__(self):
        if self.M not in _M_MODES:
            raise ValueError(f"M must be one of {_M_MODES}, got {self.M!r}")
        if self.J_range[0] > self.J_range[1]:
            raise ValueError(f"J_range must be (lo, hi) with lo <= hi, got {self.J_range}")


def thf_spec():
    """The v1 ThF+ X 3Delta1 spec: Omega = +-1, J = 1-4, I(19F) = 1/2.

    docs/thf-plus-x3delta1-effective-hamiltonian.md S1.1 (Ng thesis App. C.3.1,
    PDF p. 323). Lambda = +2, S = 1, Sigma = -1 are carried as case-(a)
    bookkeeping labels only; no case (c) element uses them except the parity phase.
    """
    return StateSpec(
        case="c",
        electronic=(ElecState(label="X3Delta1", Omega=1.0, S=1.0, Lam=2.0, T0=0.0),),
        I=0.5,
        J_range=(1, 4),
        v=0,
        M="blocks",
        frame="rotating",
    )


def enumerate_kets(spec):
    """Enumerate |J, Omega, F, m_F> for a case (c) spec.

    Deterministic order: sorted by (J, F, m_F, Omega), so the two Omega partners
    of one (J, F, m_F) level are adjacent -- which is what makes the Omega-flip
    permutation used by the parity operator a simple pairwise swap.
    """
    if spec.case != "c":
        raise NotImplementedError(f"no enumerator registered for case {spec.case!r}")
    if len(spec.electronic) != 1:
        raise NotImplementedError("v1 enumerates exactly one electronic state")
    Om0 = abs(spec.electronic[0].Omega)
    I = spec.I
    rows = []
    for J in range(spec.J_range[0], spec.J_range[1] + 1):
        for Om in (-Om0, Om0):
            if abs(Om) > J:
                continue
            for twoF in range(int(round(2 * abs(J - I))), int(round(2 * (J + I))) + 1, 2):
                for twomF in range(-twoF, twoF + 1, 2):
                    rows.append((float(J), float(Om), twoF / 2.0, twomF / 2.0))
    kets = np.array(rows, dtype=KET_C)
    kets.sort(order=["J", "F", "mF", "Om"])
    return kets


def check_basis_invariants(kets, spec):
    """Raise ValueError on any basis invariant violation (gate A7).

    Structural checks raise; nothing here is a hash or a snapshot.
    """
    for name in KET_C.names:
        v = np.asarray(kets[name], dtype=float)
        if not np.array_equal(2.0 * v, np.round(2.0 * v)):
            raise ValueError(f"field {name!r} holds a value that is not a multiple of 0.5")
    J, Om, F, mF = (np.asarray(kets[n], dtype=float) for n in ("J", "Om", "F", "mF"))
    if np.any(np.abs(Om) > J):
        raise ValueError("|Omega| > J for at least one ket")
    if np.any(F < np.abs(J - spec.I) - 1e-12) or np.any(F > J + spec.I + 1e-12):
        raise ValueError("F violates the triangle condition |J - I| <= F <= J + I")
    if np.any(np.abs(mF) > F + 1e-12):
        raise ValueError("|m_F| > F for at least one ket")
    rows = {tuple(k) for k in kets}
    if len(rows) != len(kets):
        raise ValueError("duplicate kets in the basis")


@dataclass(frozen=True)
class Blocking:
    """A named, invertible partition of the basis, with index maps back into it.

    Blocking is a first-class object, not a hidden loop, because two v1
    requirements cross blocks: an E1 transition dipole connects different m_F,
    and a transverse or rotating field couples Delta m_F = +-1 (spec S3.1).
    """
    kind: str
    labels: tuple
    index: Mapping
    n_total: int

    def merge(self, labels):
        """Concatenated index array for several blocks, in the order given.

        Reproduces C2V-Molecules' concatenated-with-offsets full basis
        (atm_core/transverse.py:83 build_full_basis) without stitching cached
        blocks -- the caller re-assembles, which avoids a class of index bug.
        """
        return np.concatenate([np.asarray(self.index[label], dtype=int) for label in labels])


def block_by_mF(kets):
    """One block per SIGNED m_F. Exact for collinear static E || B || z."""
    values = np.asarray(kets["mF"], dtype=float)
    labels = tuple(sorted({float(v) for v in values}))
    index = {label: np.flatnonzero(values == label) for label in labels}
    return Blocking(kind="signed_mF", labels=labels, index=index, n_total=len(kets))
