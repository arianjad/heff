"""The convention contract, in exactly one module.

Every phase in the package comes from here (spec S3.2). Defaults are the v1
block agreed with Arian on 2026-09-05 and recorded in the plan's Global
Constraints; each has a live fork in the ThF+ literature, listed with its
open-item number in docs/thf-plus-x3delta1-effective-hamiltonian.md S7.
"""
from dataclasses import dataclass

import numpy as np

_ALLOWED = {
    "n_hat": ("F_to_Th", "Th_to_F"),
    "ef_rule": ("brown1975", "thesis_S_half"),
    "zeeman_sign": ("plus_Gpar", "minus_Gpar"),
    "zeeman_energy": ("E_minus_g_muB_B_mF",),
    "dg_def": ("Delta", "delta"),
    "edm_factor": ("ng", "leanhardt_half"),
    "dipole_origin": ("center_of_mass", "heavy_nucleus"),
    "formalism": ("R2", "N2"),
}


@dataclass(frozen=True)
class Conventions:
    """Which fork of each documented ThF+ convention this calculation uses.

    n_hat        'F_to_Th' = JILA (Leanhardt p.15, Ng thesis p.318): Omega = +1
                 for Lambda = +2, d_mf = +3.37 D. 'Th_to_F' = Petrov/Skripnikov
                 (arXiv:2503.02840 p.3), which flips Omega, signed d and E_eff.
                 [HAM] OPEN-11.
    ef_rule      'brown1975' = B&C PDF p.283: e <=> parity +(-1)^J for integer J,
                 +(-1)^(J-1/2) for half-integer J. 'thesis_S_half' = the thesis
                 rule P = (-1)^(J-S-l), which is the S = 1/2 specialisation and
                 INVERTS every label at S = 1. [HAM] OPEN-2.
    zeeman_sign  'plus_Gpar' = +G_par mu_B (J.n)(n.B) - g_N mu_N I.B. Ng thesis
                 Eq. C.6 prints a minus; the minus contradicts the g_F formula
                 three lines below it, his own G_par = 0.048(2), and the measured
                 |g| = 0.0149. [HAM] S2.8, OPEN-3.
    dg_def       'Delta' = g^u - g^l (Petrov, Ng 2022 paper); 'delta' = half that
                 (Ng thesis App. C.4). [HAM] OPEN-12.
    edm_factor   'ng' = -(d_e E_eff + W_TP k_TP) Omega/|Omega| (Ng Eq. C.8, no
                 Leanhardt 1/2). [HAM] S2.12.
    dipole_origin 'center_of_mass' -- d_mf is origin-dependent for an ION; the
                 Th-nucleus origin differs by 0.72 D, a 20 % Stark error.
                 [HAM] S2.7.
    """
    n_hat: str = "F_to_Th"
    ef_rule: str = "brown1975"
    zeeman_sign: str = "plus_Gpar"
    zeeman_energy: str = "E_minus_g_muB_B_mF"
    dg_def: str = "Delta"
    edm_factor: str = "ng"
    dipole_origin: str = "center_of_mass"
    formalism: str = "R2"
    version: str = "thf-v1"

    def __post_init__(self):
        for name, allowed in _ALLOWED.items():
            got = getattr(self, name)
            if got not in allowed:
                raise ValueError(f"{name}={got!r} is not one of {allowed}")

    def stamp(self):
        """The dict that every result object carries (spec S3.5, S3.6)."""
        return {f: getattr(self, f) for f in _ALLOWED} | {"version": self.version}


def n_hat_sign(conv):
    """+1 for the JILA n_hat (F -> Th), -1 for Petrov's (Th -> F).

    Multiplies the signed molecule-frame dipole and flips the sign of Omega's
    physical meaning. [HAM] OPEN-11.
    """
    return 1.0 if conv.n_hat == "F_to_Th" else -1.0


def parity_phase(J, S, *, ell, s):
    """The composite case-(a) parity phase (-1)^(J - S - l + s).

    E* = sigma_xz R_y(pi); thesis Eq. A.15 = B&C Eq. (6.234) (PDF p.283) in the
    form (-1)^(J-S+s). For 3Delta (S = 1, l = 0, s = 0) it is (-1)^(J-1).
    [HAM] S2 convention contract.

    ell and s are required keyword-only arguments -- no default, so a Sigma-
    or a bending state cannot silently inherit a linear-molecule parity.
    """
    return float((-1.0) ** (J - S - ell + s))


def parity_operator(kets, S, *, ell, s):
    """The parity operator in the signed-Omega primitive basis.

    P |J, Omega, F, m_F> = (-1)^(J-S-l+s) |J, -Omega, F, m_F>. Real symmetric,
    P^2 = 1. Used by gate B2 / [HAM] V8 to check that a field-free term set is
    parity-conserving and that the Stark and eEDM operators are parity-odd.

    ell and s are required keyword-only arguments -- no default, so a Sigma-
    or a bending state cannot silently inherit a linear-molecule parity.
    """
    n = len(kets)
    P = np.zeros((n, n))
    key = {(kets["J"][i], kets["Om"][i], kets["F"][i], kets["mF"][i]): i for i in range(n)}
    for i in range(n):
        j = key[(kets["J"][i], -kets["Om"][i], kets["F"][i], kets["mF"][i])]
        P[j, i] = parity_phase(kets["J"][i], S, ell=ell, s=s)
    return P


def superposition_parity(J, sym, S, *, ell, s):
    """Parity of (|+Omega> + sym |-Omega>)/sqrt(2), sym = +1 or -1.

    ell and s are required keyword-only arguments -- no default, so a Sigma-
    or a bending state cannot silently inherit a linear-molecule parity.
    """
    return int(round(sym * parity_phase(J, S, ell=ell, s=s)))


def ef_label(J, parity, *, rule, S=None, ell):
    """'e' or 'f' for a level of total angular momentum J and parity +-1.

    ell is a required keyword-only argument -- no default, so a Sigma- or a
    bending state cannot silently inherit a linear-molecule parity.
    """
    if rule == "brown1975":
        integral = abs(J - round(J)) < 1e-9
        ref = (-1.0) ** J if integral else (-1.0) ** (J - 0.5)
    elif rule == "thesis_S_half":
        if S is None:
            raise ValueError("rule 'thesis_S_half' needs S; it is P = (-1)^(J-S-l)")
        ref = (-1.0) ** (J - S - ell)
    else:
        raise ValueError(f"unknown ef_rule {rule!r}")
    return "e" if round(parity) == round(ref) else "f"
