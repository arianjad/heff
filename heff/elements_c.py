"""Case (c) matrix elements for |J, Omega, F, m_F>.

EVERY element in this file is copied from
docs/thf-plus-x3delta1-effective-hamiltonian.md with its primary-source citation
carried into the decorator. Nothing here was re-derived. If a check disagrees
with one of these forms, report both with citations -- do not adjust a sign
(Arian's standing rule; [HAM] S7 lists the fifteen open items).

Sign conventions in force (heff.conventions.Conventions defaults):
  n_hat = 'F_to_Th'  (JILA): Omega = +1 for Lambda = +2, d_mf = +3.37 D
  zeeman_sign = 'plus_Gpar': +G_par mu_B (J.n)(n.B), NOT as Ng Eq. C.6 prints it
  edm_factor = 'ng': no Leanhardt 1/2
  parity: E* |J, Omega> = (-1)^(J-S+s) |J, -Omega>  (heff.conventions.parity_phase)
"""
import numpy as np

from .terms import Rules, term


def _ph(x):
    """(-1)**x for an x that must be an integer.

    Guarded on purpose: every phase exponent in this file is an integer (a sum
    like F + J' + 1 + I with two half-integers), and `(-1.0) ** 2.5` in Python
    is a complex number, so a half-integer slipping in would silently poison the
    matrix instead of raising.
    """
    n = round(float(x))
    if abs(float(x) - n) > 1e-9:
        raise ValueError(f"phase exponent {x!r} is not an integer")
    return 1.0 if n % 2 == 0 else -1.0


def _same(bra, ket, *fields):
    return all(bra[f] == ket[f] for f in fields)


# ----------------------------------------------------------------- rotation

@term(name="rotation", param=("B0",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.7c p.322: H_rot = B_v J(J+1). Case (c), NO -Omega^2 "
           "term (Ng Ch. 2.4 Eq. 2.2a p.37) -- the Omega-dependent constant is "
           "absorbed into the band origin, which is why 4B = 29.09733(4) GHz is a "
           "literal statement about the J = 1 -> 2 interval. [HAM] S2.1")
def rotation(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "F", "mF"):
        return 0.0
    J = float(ket["J"])
    return J * (J + 1.0)


@term(name="centrifugal", param=("D0",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.7d p.322: -D [J(J+1)]^2. The minus lives in the "
           "matrix, so D0 stays positive (3.897 kHz). [HAM] S2.1")
def centrifugal(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "F", "mF"):
        return 0.0
    J = float(ket["J"])
    return -((J * (J + 1.0)) ** 2)


# ------------------------------------------------------------ Omega-doubling

@term(name="omega_doubling", param=("omega_ef",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(-2.0, 2.0), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.3 p.319, sign transposed to the heff ket phase "
           "(-1)^(J-S+s); see docs/open-questions.md OQ-A. Ng writes "
           "H = ((-1)^J / 2) hbar omega_ef Omega_x^(J) with "
           "Omega_x^(J) = [J(J+1)/2](|+1><-1| + |-1><+1|). In heff's convention "
           "the parity operator is E*|J,Om> = (-1)^(J-S+s)|J,-Om> = (-1)^(J-1)|J,-Om> "
           "for a 3Delta, and the physical invariant -- upper doublet component "
           "of parity (-1)^J at every J, i.e. e above f uniformly under Brown's "
           "1975 rule (Ng 2022 Fig. 2 read as an image, [HAM] S2.3; Gresh k'' < 0) "
           "-- requires the J-INDEPENDENT off-diagonal -omega_ef J(J+1)/4. Ng's "
           "(-1)^J is the same physics in a ket phase where E*|J,Om> = |J,-Om>. "
           "The splitting law omega_ef J(J+1)/2 is identical either way and is "
           "hard-gated. [HAM] S2.3")
def omega_doubling(bra, ket, ctx):
    if not _same(bra, ket, "J", "F", "mF"):
        return 0.0
    if bra["Om"] != -ket["Om"]:
        return 0.0
    J = float(ket["J"])
    return -J * (J + 1.0) / 4.0


# ----------------------------------------------------------------- hyperfine

@term(name="hyperfine_A_par", param=("A_par",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.2 p.319 = B&C Eq. (9.50) PDF p.636 / book p.604: "
           "A_par [F(F+1) - I(I+1) - J(J+1)] / [2 J(J+1)], with "
           "A_par = {a Lambda + (b_F + (2/3) c) Sigma} Omega = 2a - b_F - (2/3)c "
           "for 3Delta1 (Lambda=+2, Sigma=-1, Omega=+1). [HAM] S2.4")
def hyperfine_A_par(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "F", "mF"):
        return 0.0
    J, F, I = float(ket["J"]), float(ket["F"]), float(ctx.I)
    return (F * (F + 1.0) - I * (I + 1.0) - J * (J + 1.0)) / (2.0 * J * (J + 1.0))


@term(name="hyperfine_A_par_dJ1", param=("A_par",), cases=("c",),
      rules=Rules(dJ=(-1, 1), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="B&C Eq. (9.51) PDF p.636 / book p.604, <eta,J|H_hf|eta,J-1> = "
           "-{a Lambda + (b_F + (2/3)c) Sigma} (J^2-Omega^2)^(1/2) "
           "{(F-I+J)(F+I+J+1)(J+I-F)(F-J+I+1)}^(1/2) / [2J(4J^2-1)^(1/2)], "
           "with J the LARGER of the two and the brace = A_par / Omega. The "
           "term Ng App. C drops; 1.6-2.6 kHz, above this model's kHz floor. "
           "[HAM] S2.5")
def hyperfine_A_par_dJ1(bra, ket, ctx):
    if not _same(bra, ket, "Om", "F", "mF"):
        return 0.0
    if abs(float(bra["J"]) - float(ket["J"])) != 1.0:
        return 0.0
    J = max(float(bra["J"]), float(ket["J"]))
    Om, F, I = float(ket["Om"]), float(ket["F"]), float(ctx.I)
    rad = (F - I + J) * (F + I + J + 1.0) * (J + I - F) * (F - J + I + 1.0)
    if rad <= 0.0:
        return 0.0
    return -(1.0 / Om) * np.sqrt(J * J - Om * Om) * np.sqrt(rad) / (
        2.0 * J * np.sqrt(4.0 * J * J - 1.0))


@term(name="spin_rotation_cI", param=("c_I",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="B&C Eq. (8.7) PDF p.410 H_nsr = c_I T1(J).T1(I) with its coupled-basis "
           "element Eq. (8.20) PDF p.414, which reduces to "
           "c_I [F(F+1) - I(I+1) - J(J+1)] / 2. c_I is an ESTIMATE (~20 kHz, "
           "factor 3) from the CsF 19F anchor, and it may bias the fitted A_par "
           "by 40-120 kHz. [HAM] S2.6, OPEN-6")
def spin_rotation_cI(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "F", "mF"):
        return 0.0
    J, F, I = float(ket["J"]), float(ket["F"]), float(ctx.I)
    return (F * (F + 1.0) - I * (I + 1.0) - J * (J + 1.0)) / 2.0
