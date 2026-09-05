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
from dataclasses import dataclass
from typing import Mapping

import numpy as np

KET_C = np.dtype([("J", "f8"), ("Om", "f8"), ("F", "f8"), ("mF", "f8")])

# v2 (spec-v2 S2.2): one added field, not a variable-width dtype. The coupling
# order is inner-spin-first, |((J I_1) F1, I_2) F, m_F>; a different order is a
# different basis (9j-related, different spectator phases), so it is fixed here
# and nowhere else.
KET_C2 = np.dtype([("J", "f8"), ("Om", "f8"), ("F1", "f8"), ("F", "f8"), ("mF", "f8")])

_M_MODES = ("none", "blocks", "all")

_ISOTOPOLOGUES = ("232", "229", "227")


@dataclass(frozen=True)
class ElecState:
    """One electronic state. No defaults for Omega, S, Lambda (spec S3.2)."""
    label: str
    Omega: float
    S: float
    Lam: float
    T0: float = 0.0


@dataclass(frozen=True)
class Spin:
    """One coupled nuclear spin. No defaults (spec S3.2).

    `couple_to` is the label of the angular momentum this spin adds to: 'J' for
    the innermost spin, then 'F1', 'F2', ... for each subsequent one. It is
    data; the enumerator is the only code that reads it.
    """
    label: str
    I: float
    couple_to: str


@dataclass(frozen=True)
class StateSpec:
    case: str
    electronic: tuple
    I: float
    J_range: tuple
    v: int = 0
    M: str = "blocks"
    # A DOCUMENTED LABEL ONLY: no code path reads `frame` yet -- it is stamped
    # into the term-matrix manifest and nothing else. The trigger for wiring it
    # is the rotating-frame term (omega_rot F_x), which is the first operator
    # whose matrix elements differ between 'rotating' and 'lab' (spec S3.1).
    frame: str = "rotating"
    # The coupled nuclear spins, inner-first. Empty is the v1 basis, where `I`
    # alone carries the single (19F) spin; when `spins` is non-empty it is the
    # authority and `I` is bookkeeping.
    spins: tuple = ()

    def __post_init__(self):
        if self.M not in _M_MODES:
            raise ValueError(f"M must be one of {_M_MODES}, got {self.M!r}")
        if self.J_range[0] > self.J_range[1]:
            raise ValueError(f"J_range must be (lo, hi) with lo <= hi, got {self.J_range}")
        for i, s in enumerate(self.spins):
            want = "J" if i == 0 else f"F{i}"
            if s.couple_to != want:
                raise ValueError(
                    f"spins must form a chain coupled inner-first: spins[{i}] "
                    f"({s.label!r}) must have couple_to == {want!r}, got {s.couple_to!r}")


def thf_spec(isotopologue=None, *, J_max=4):
    """The ThF+ X 3Delta1 spec: Omega = +-1, J = 1..J_max, I(19F) = 1/2.

    docs/thf-plus-x3delta1-effective-hamiltonian.md S1.1 (Ng thesis App. C.3.1,
    PDF p. 323). Lambda = +2, S = 1, Sigma = -1 are carried as case-(a)
    bookkeeping labels only; no case (c) element uses them except the parity phase.

    `isotopologue` is None (the v1 default) or one of '232' | '229' | '227'.
    None and '232' return the same object: 232Th is spin-0, so the basis is the
    v1 |J, Omega, F, m_F>. '229' adds I(229Th) = 5/2 and '227' adds I(227Th) =
    1/2, both as an inner spin coupled to J, giving |((J I_Th) F1, I_F) F, m_F>.

    The 227Th value is the ENSDF-tentative (1/2+) ground-state assignment. It
    rests on decay-scheme systematics, not on any moment or hyperfine
    measurement -- no measured or estimated mu(227Th) exists in Stone's
    compilations or the IAEA NDS moments database
    (docs/lit/lookup-227th-nuclear-moment.md). The 9.3 keV level, 5/2+, is the
    live alternative: if the assignment moves there, I_Th = 2.5 here.

    `J_max` is a keyword because the truncation is a physics choice -- the Th
    DeltaJ = +-1 hyperfine is a ~2 GHz off-diagonal element -- not a constant.
    """
    if isotopologue is not None and isotopologue not in _ISOTOPOLOGUES:
        raise ValueError(
            f"isotopologue must be None or one of {_ISOTOPOLOGUES}, got {isotopologue!r}")
    I_Th = {"229": 2.5, "227": 0.5}.get(isotopologue)
    spins = () if I_Th is None else (
        Spin(label=f"{isotopologue}Th", I=I_Th, couple_to="J"),
        Spin(label="19F", I=0.5, couple_to="F1"),
    )
    return StateSpec(
        case="c",
        electronic=(ElecState(label="X3Delta1", Omega=1.0, S=1.0, Lam=2.0, T0=0.0),),
        I=0.5,
        J_range=(1, J_max),
        v=0,
        M="blocks",
        frame="rotating",
        spins=spins,
    )


def enumerate_kets(spec):
    """Enumerate |J, Omega, F, m_F> for a case (c) spec.

    Deterministic order: sorted by (J, F, m_F, Omega), so the two Omega partners
    of one (J, F, m_F) level are adjacent -- which is what makes the Omega-flip
    permutation used by the parity operator a simple pairwise swap.

    With two coupled spins on the spec the basis is |((J I_1) F1, I_2) F, m_F>,
    dtype KET_C2, sorted by (J, F1, F, m_F, Omega) -- the same guarantee, one
    key deeper. At I_1 = 0 the F1 column is identically J and the order
    collapses to the v1 order, which is what makes the v1 basis a special case
    of this one rather than a different basis (spec-v2 S2.1).
    """
    if spec.case != "c":
        raise NotImplementedError(f"no enumerator registered for case {spec.case!r}")
    if len(spec.electronic) != 1:
        raise NotImplementedError("v1 enumerates exactly one electronic state")
    if len(spec.spins) == 2:
        return _enumerate_kets_two_spin(spec)
    if len(spec.spins) > 2:
        raise NotImplementedError(
            "no enumerator for more than two coupled spins; add a third dtype when "
            "a species with three coupled nuclear spins arrives")
    Om0 = abs(spec.electronic[0].Omega)
    # `spins` is the authority when it is there; `I` is the v1 spelling of the
    # same single spin. Identical for spins=() -- the v1 path is untouched.
    I = spec.spins[0].I if spec.spins else spec.I
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


def _enumerate_kets_two_spin(spec):
    """|((J I_1) F1, I_2) F, m_F>. Half-integers as integer 2x counters, as v1 does."""
    Om0 = abs(spec.electronic[0].Omega)
    twoI1 = int(round(2 * spec.spins[0].I))
    twoI2 = int(round(2 * spec.spins[1].I))
    rows = []
    for J in range(spec.J_range[0], spec.J_range[1] + 1):
        twoJ = 2 * J
        for Om in (-Om0, Om0):
            if abs(Om) > J:
                continue
            for twoF1 in range(abs(twoJ - twoI1), twoJ + twoI1 + 1, 2):
                for twoF in range(abs(twoF1 - twoI2), twoF1 + twoI2 + 1, 2):
                    for twomF in range(-twoF, twoF + 1, 2):
                        rows.append((float(J), float(Om), twoF1 / 2.0,
                                     twoF / 2.0, twomF / 2.0))
    kets = np.array(rows, dtype=KET_C2)
    kets.sort(order=["J", "F1", "F", "mF", "Om"])
    return kets


def check_basis_invariants(kets, spec):
    """Raise ValueError on any basis invariant violation (gate A7).

    Structural checks raise; nothing here is a hash or a snapshot.
    """
    for name in kets.dtype.names:
        v = np.asarray(kets[name], dtype=float)
        if not np.array_equal(2.0 * v, np.round(2.0 * v)):
            raise ValueError(f"field {name!r} holds a value that is not a multiple of 0.5")
    J, Om, F, mF = (np.asarray(kets[n], dtype=float) for n in ("J", "Om", "F", "mF"))
    if np.any(np.abs(Om) > J):
        raise ValueError("|Omega| > J for at least one ket")
    if "F1" in kets.dtype.names:
        # v2: F couples to F1, not to J, and F1 carries its own triangle.
        I1, I2 = spec.spins[0].I, spec.spins[1].I
        F1 = np.asarray(kets["F1"], dtype=float)
        if np.any(F1 < np.abs(J - I1) - 1e-12) or np.any(F1 > J + I1 + 1e-12):
            raise ValueError(
                "F1 violates the triangle condition |J - I_1| <= F1 <= J + I_1")
        lo, hi, parent = np.abs(F1 - I2), F1 + I2, "F1"
    else:
        lo, hi, parent = np.abs(J - spec.I), J + spec.I, "J"
    if np.any(F < lo - 1e-12) or np.any(F > hi + 1e-12):
        raise ValueError(
            f"F violates the triangle condition |{parent} - I| <= F <= {parent} + I")
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


def blocks_for(spec, kets):
    """The Blocking a spec's `M` mode asks for (spec S3.1).

    'blocks' -> one block per signed m_F, the collinear fast path and the
    default; 'all' -> a single 'all' block holding the whole basis, which is
    what a Delta-m_F != 0 term (transverse or rotating field) must be built on;
    'none' -> the m-free field-free basis, not implemented.

    `block_by_mF` keeps its own signature for the callers that want the
    per-m_F partition regardless of the spec.
    """
    if spec.M == "blocks":
        return block_by_mF(kets)
    if spec.M == "all":
        return Blocking(kind="all", labels=("all",),
                        index={"all": np.arange(len(kets))}, n_total=len(kets))
    raise NotImplementedError(
        "M='none' is the field-free spectra path; add when a consumer needs an "
        "m-free basis")
