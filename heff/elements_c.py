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

from .conventions import n_hat_sign
from .terms import Rules, term
from .wigner import w3j, w6j


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
           "absorbed into the band origin. 4B is the rotational-only J=1->2 interval; "
           "with the centrifugal term it is 4B-32D. The stored B is interval/4 "
           "from Ng, not a joint refit; see the 2026-09-08 estimate audit. [HAM] S2.1")
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


# --------------------------------------------------- the shared E1 geometry

def dipole_geometry(bra, ket, I, p):
    """The three-factor geometry of a rank-1 molecule-frame operator along n_hat.

    Ng thesis Eq. C.5 p.320 (H_Stark = -d_mf n_hat . E, q = 0 because the dipole
    lies along the axis), independently re-derived in [HAM] S2 from B&C Eqs.
    (5.174) + (5.172) + (5.186) and found identical term for term, including
    every phase, with the two 6j symbols related by column exchange:

      <J',Om',F',m'| . |J,Om,F,m>
        = (-1)^(F + J' + 1 + I)  { J  F  I ;  F' J' 1 }
        x (-1)^(F' - m') sqrt((2F+1)(2F'+1)) ( F'  1  F ; -m'  p  m )
        x (-1)^(J' - Om') sqrt((2J+1)(2J'+1)) ( J'  1  J ; -Om' 0  Om )

    Primed = bra. Diagonal in Omega. Returns a dimensionless float; the caller
    supplies -d_mf E_p (Stark), +G_par mu_B B_p Omega (Zeeman) or the E1
    transition dipole (heff.spectra). On the diagonal it evaluates to
    <n_hat_z> = Omega m_F gamma_F ([HAM] S2.7 closed form, S2.8).

    Public because heff.spectra calls it with p = +-1 while the Stark term uses
    p = 0 -- one formula, one place, per spec S2.1(3).
    """
    if bra["Om"] != ket["Om"]:
        return 0.0
    J1, F1, m1 = float(ket["J"]), float(ket["F"]), float(ket["mF"])
    J2, F2, m2 = float(bra["J"]), float(bra["F"]), float(bra["mF"])
    Om = float(ket["Om"])
    if abs(m2 - m1 - p) > 1e-9:
        return 0.0
    six = w6j(J1, F1, I, F2, J2, 1)
    if six == 0.0:
        return 0.0
    a = _ph(F1 + J2 + 1.0 + I) * six
    b = (_ph(F2 - m2) * np.sqrt((2 * F1 + 1.0) * (2 * F2 + 1.0))
         * w3j(F2, 1, F1, -m2, p, m1))
    c = (_ph(J2 - Om) * np.sqrt((2 * J1 + 1.0) * (2 * J2 + 1.0))
         * w3j(J2, 1, J1, -Om, 0, Om))
    return a * b * c


# --------------------------------------------------------------------- Stark

@term(name="stark_z", param=("d_mf", "E_z"), cases=("c",),
      rules=Rules(dJ=(-1, 0, 1), dOm=(0.0,), dF=(-1, 0, 1), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.5 p.320, H_Stark = -d_mf n_hat . E; re-derived in "
           "[HAM] S2 from B&C (5.174)+(5.172)+(5.186) and found identical. "
           "PARITY-ODD: within a J it connects the two parity components of the "
           "Omega-doublet, which is the whole polarisation mechanism. The "
           "Delta-J = +-1 part is not a nuisance -- it is what generates delta-g "
           "([HAM] S2.7, S2.9). Odd in n_hat, so it carries n_hat_sign: with "
           "'Th_to_F' (Petrov/Skripnikov) the signed dipole flips, [HAM] OPEN-11.")
def stark_z(bra, ket, ctx):
    return -n_hat_sign(ctx.conventions) * dipole_geometry(bra, ket, ctx.I, 0)


# -------------------------------------------------------------------- Zeeman

@term(name="zeeman_Gpar", param=("G_par", "B_z"), cases=("c",),
      rules=Rules(dJ=(-1, 0, 1), dOm=(0.0,), dF=(-1, 0, 1), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.6 p.321 with the sign CORRECTED to "
           "+G_par mu_B (J.n_hat)(n_hat.B): the printed minus contradicts the "
           "g_F formula three lines below it, Ng's own G_par = 0.048(2) / "
           "-0.042(2) pair, and the measured |g_F=3/2| = 0.0149(3). "
           "[HAM] S2.8, OPEN-3. Same geometry as the Stark term with "
           "-d_mf E_p -> +G_par mu_B B_p Omega, as Ng notes himself. "
           "PARITY-EVEN and EVEN in n_hat (quadratic in n_hat), so unlike the "
           "Stark and PT-odd terms it does NOT carry n_hat_sign.")
def zeeman_Gpar(bra, ket, ctx):
    sign = 1.0 if ctx.conventions.zeeman_sign == "plus_Gpar" else -1.0
    return sign * ctx.mu_B * float(ket["Om"]) * dipole_geometry(bra, ket, ctx.I, 0)


@term(name="zeeman_nuclear", param=("g_N", "B_z"), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(-1, 0, 1), dmF=(0,)),
      hermitian=True, real=True,
      cite="H = -g_N mu_N I.B (Ng thesis Eq. C.6 p.321), with T1_p(I) from B&C's "
           "spectator theorem Eq. (5.175) on the I part: "
           "(-1)^(F'+J+I+1) sqrt(I(I+1)(2I+1)) {I F' J; F I 1} x "
           "(-1)^(F'-m') sqrt((2F+1)(2F'+1)) (F' 1 F; -m' p m). "
           "g_N mu_N = 4.008 kHz/G -- a kHz-scale effect at 1 G that must be "
           "kept. [HAM] S2.8")
def zeeman_nuclear(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "mF"):
        return 0.0
    I = float(ctx.I)
    J, F1, m1 = float(ket["J"]), float(ket["F"]), float(ket["mF"])
    F2, m2 = float(bra["F"]), float(bra["mF"])
    six = w6j(I, F2, J, F1, I, 1)
    if six == 0.0:
        return 0.0
    a = _ph(F2 + J + I + 1.0) * np.sqrt(I * (I + 1.0) * (2 * I + 1.0)) * six
    b = (_ph(F2 - m2) * np.sqrt((2 * F1 + 1.0) * (2 * F2 + 1.0))
         * w3j(F2, 1, F1, -m2, 0, m1))
    return -ctx.mu_N * a * b


# -------------------------------------------------------------------- PT-odd

def _pt_odd_matrix(bra, ket, ctx):
    """The shared -Omega/|Omega| matrix of [HAM] S2.12, opt-in via its knob.

    Odd in n_hat: flipping to 'Th_to_F' flips the sign of Omega's physical
    meaning and with it E_eff and W_TP ([HAM] OPEN-11, and Skripnikov & Titov
    2015 p.2, whose Omega = <Psi|J.n|Psi> is defined with n from Th to F).
    """
    if not _same(bra, ket, "J", "Om", "F", "mF"):
        return 0.0
    half = 0.5 if ctx.conventions.edm_factor == "leanhardt_half" else 1.0
    return -half * n_hat_sign(ctx.conventions) * float(np.sign(ket["Om"]))


@term(name="pt_odd_edm", param=("d_e", "E_eff"), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.8 p.322: H_eEDM = -d_e E_eff Omega/|Omega|, with NO "
           "Leanhardt 1/2 (Ng footnote 3 p.322: the 1/2 in Leanhardt Eq. 26 "
           "reflects a different E_eff convention and is removed for consistency "
           "with the cited value). The observable is f^BD = 2 d_e E_eff. "
           "Opt-in: d_e defaults to 0. [HAM] S2.12")
def pt_odd_edm(bra, ket, ctx):
    return _pt_odd_matrix(bra, ket, ctx)


@term(name="pt_odd_scalar_pseudoscalar", param=("W_TP", "k_TP"), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng Eq. C.8 form, same added matrix as the eEDM term with "
           "d_e E_eff -> W_TP k_TP (Skripnikov & Titov 2015 Eqs. 1-4; W_TP is "
           "tabulated in kHz so the shift is W_TP k_TP). Experimentally "
           "degenerate with the eEDM in a single species. Opt-in: k_TP defaults "
           "to 0. [HAM] S2.12")
def pt_odd_scalar_pseudoscalar(bra, ket, ctx):
    return _pt_odd_matrix(bra, ket, ctx)
