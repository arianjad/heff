"""Case-(c) elements for |J, Omega, F, m_F>; see the cited decorators and [HAM].

Defaults use F→Th n_hat, +G_zz mu_B (J.n)(n.B), Ng's eEDM normalization, and
E*|J,Omega> = (-1)^(J-S+s)|J,-Omega>.
"""
import numpy as np

from .conventions import n_hat_sign
from .terms import Rules, term
from .wigner import w3j, w6j


def _ph(x):
    """Return (-1)**x, rejecting noninteger phase exponents."""
    n = round(float(x))
    if abs(float(x) - n) > 1e-9:
        raise ValueError(f"phase exponent {x!r} is not an integer")
    return 1.0 if n % 2 == 0 else -1.0


def _same(bra, ket, *fields):
    return all(bra[f] == ket[f] for f in fields)



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



@term(name="omega_doubling", param=("omega_ef",), cases=("c",),
      rules=Rules(dJ=(0,), dOm=(-2.0, 2.0), dF=(0,), dmF=(0,)),
      hermitian=True, real=True,
      cite="Ng thesis Eq. C.3 p.319, sign transposed to the heff ket phase "
           "(-1)^(J-S+s); see docs/open-questions.md OQ-A. Ng writes "
           "H = ((-1)^J / 2) hbar omega_ef Omega_x^(J) with "
           "Omega_x^(J) = [J(J+1)/2](|+1><-1| + |-1><+1|). In heff's convention "
           "the parity operator is E*|J,Om> = (-1)^(J-S+s)|J,-Om> = (-1)^(J-1)|J,-Om> "
           "for a 3Delta, so the ORDERING -- e (parity +(-1)^J) BELOW f at every "
           "J -- requires the J-INDEPENDENT off-diagonal +omega_ef J(J+1)/4. That "
           "ordering is Petrov & Skripnikov arXiv:2503.02840's model: the 1Sigma+ "
           "state 314 cm-1 above mixes with the e component alone (0+ states carry "
           "only e levels) and pushes it down, the same mixing Petrov fits the "
           "Omega-doubling to. It is consistent with Gresh 2016 Table 1's k'' > 0 "
           "on the Omega = 0+ <- 3Delta1 bands. Ng's schematics assume an ordering "
           "his experiment did not determine (Ng thesis Fig. 1.4 caption p.30: 'We "
           "did not determine the energy ordering of the parity states'). The "
           "splitting law omega_ef J(J+1)/2 is identical either way and is "
           "hard-gated. [HAM] S2.3")
def omega_doubling(bra, ket, ctx):
    if not _same(bra, ket, "J", "F", "mF"):
        return 0.0
    if bra["Om"] != -ket["Om"]:
        return 0.0
    J = float(ket["J"])
    return J * (J + 1.0) / 4.0



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



def inner_axial(J_bra, Om_bra, J_ket, Om_ket):
    """<J',Om'||T1(n^)||J,Om> = delta_{Om Om'} (-1)^(J'-Om) [(2J+1)(2J'+1)]^(1/2)
    (J' 1 J; -Om 0 Om) -- B&C Eq. (5.186) PDF p.207 / book p.175.

    The Omega-diagonal early return lives HERE and not in dipole_geometry, so
    that a non-axial `inner` (the Delta-Omega = +-2 Zeeman channel) reaches the
    same spectator recoupling. [HAM] S2.7.
    """
    if Om_bra != Om_ket:
        return 0.0
    return (_ph(J_bra - Om_ket) * np.sqrt((2 * J_ket + 1.0) * (2 * J_bra + 1.0))
            * w3j(J_bra, 1, J_ket, -Om_ket, 0, Om_ket))


def inner_J(J_bra, Om_bra, J_ket, Om_ket):
    """<J',Om'||T1(J)||J,Om> = delta_JJ' delta_{Om Om'} [J(J+1)(2J+1)]^(1/2).

    B&C Eq. (5.174) PDF p.205 / book p.173 closed on the standard
    angular-momentum reduced element Eq. (5.179). This is the J-space factor of
    the lab operator T1_0(J), i.e. of B.J with B along z. [HAM] S2.8.
    """
    if J_bra != J_ket or Om_bra != Om_ket:
        return 0.0
    return np.sqrt(J_ket * (J_ket + 1.0) * (2 * J_ket + 1.0))


def inner_flip(J_bra, Om_bra, J_ket, Om_ket):
    """J-space factor of C = B_x J_x - B_y J_y, the Delta-Omega = +-2 channel.

    In body spherical components C = (1/sqrt(2)) [D1*_{0,-1} J^b_+ -
    D1*_{0,+1} J^b_-] (B&C Eq. (9.70) term (vii) and Eq. (9.71) PDF pp.652-653;
    body components of J have the ANOMALOUS commutation of B&C Eqs.
    (5.152)-(5.153) PDF p.200, so J^b_+ LOWERS Omega, and the operator is
    evaluated Brown-Howard style, Eqs. (5.155)-(5.162) pp.201-202). The relative
    minus is forced by Hermiticity, (D1*_{0,-1})^dagger = -D1*_{0,+1}; without
    it the operator is anti-Hermitian. The two factors in each product commute,
    so ladder-on-bra and ladder-on-ket orderings agree.

    Writing q = Om' - Om = +-2 and s = q/2, closure over the intermediate
    Omega + s leaves one D1 element and one ladder factor:

      -sign(q) sqrt(2) [J(J+1) - Om(Om+s)]^(1/2) (-1)^(J'-Om')
      [(2J+1)(2J'+1)]^(1/2) (J' 1 J; -Om' s Om+s)

    C is a lab rank-1 tensor times a lab scalar, so the k = 1 spectator chain in
    front of it is unchanged -- a lab rank-2 geometry would be quadrupolar in
    m_F and is the wrong operator, see the 2026-09-14 handoff audit S2.1. The
    overall sqrt(2) is the normalisation that makes <J,-/+1,M|C|J,+/-1,M> = +M
    exactly and J-independently, which is what fixes the meaning of G_Delta.
    """
    q = Om_bra - Om_ket
    if abs(q) != 2.0:
        return 0.0
    s = q / 2.0
    rad = J_ket * (J_ket + 1.0) - Om_ket * (Om_ket + s)
    if rad <= 0.0:
        return 0.0
    return (-np.sign(q) * np.sqrt(2.0) * np.sqrt(rad)
            * _ph(J_bra - Om_bra)
            * np.sqrt((2 * J_ket + 1.0) * (2 * J_bra + 1.0))
            * w3j(J_bra, 1, J_ket, -Om_bra, s, Om_ket + s))


def dipole_geometry(bra, ket, I, p, *, inner=inner_axial):
    """Dimensionless rank-1 geometry (Ng C.5; B&C 5.172, 5.174, 5.186).

    Primed labels are bra. `inner(J_bra, Om_bra, J_ket, Om_ket)` is the J-space
    reduced element; the spectator recoupling in front of it is the SAME for
    every lab rank-1 operator, which is why the Zeeman tensor's three terms
    differ only in this argument. See [HAM] S2.7.
    """
    J1, F1, m1 = float(ket["J"]), float(ket["F"]), float(ket["mF"])
    J2, F2, m2 = float(bra["J"]), float(bra["F"]), float(bra["mF"])
    if abs(m2 - m1 - p) > 1e-9:
        return 0.0
    six = w6j(J1, F1, I, F2, J2, 1)
    if six == 0.0:
        return 0.0
    a = _ph(F1 + J2 + 1.0 + I) * six
    b = (_ph(F2 - m2) * np.sqrt((2 * F1 + 1.0) * (2 * F2 + 1.0))
         * w3j(F2, 1, F1, -m2, p, m1))
    c = inner(J2, float(bra["Om"]), J1, float(ket["Om"]))
    return a * b * c



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



@term(name="zeeman_Gzz", param=("G_zz", "B_z"), cases=("c",),
      rules=Rules(dJ=(-1, 0, 1), dOm=(0.0,), dF=(-1, 0, 1), dmF=(0,)),
      hermitian=True, real=True,
      cite="The zz component of H_Z = mu_B B.G.J with G = diag(G_xx, G_yy, G_zz) "
           "in the molecule frame (z along n_hat): G_zz (B.n_hat)(J.n_hat). "
           "G_zz IS Ng's G_par -- Ng thesis Eq. C.6 p.321 with the sign "
           "CORRECTED to +G_zz mu_B (J.n_hat)(n_hat.B): the printed minus "
           "contradicts the g_F formula three lines below it, Ng's own "
           "G_par = 0.048(2) / -0.042(2) pair, and the measured "
           "|g_F=3/2| = 0.0149(3). [HAM] S2.8, OPEN-3. Same geometry as the "
           "Stark term with -d_mf E_p -> +G_zz mu_B B_p Omega, as Ng notes "
           "himself. PARITY-EVEN and EVEN in n_hat (quadratic in n_hat), so "
           "unlike the Stark and PT-odd terms it does NOT carry n_hat_sign.")
def zeeman_Gzz(bra, ket, ctx):
    sign = 1.0 if ctx.conventions.zeeman_sign == "plus_Gpar" else -1.0
    return sign * ctx.mu_B * float(ket["Om"]) * dipole_geometry(bra, ket, ctx.I, 0)


def _zeeman_perp(bra, ket, ctx, flip):
    """mu_B/2 [P +- C] -- the x or y component of mu_B B.G.J ([HAM] S2.8).

    P = B.J - (B.n_hat)(J.n_hat) is the isotropic perpendicular piece that
    B_x J_x and B_y J_y share; C = B_x J_x - B_y J_y is the anisotropic,
    Delta-Omega = +-2 piece they differ in. `flip` is +1 for xx and -1 for yy,
    so G_xx and G_yy enter as (G_xx + G_yy)/2 on P and (G_xx - G_yy)/2 on C.
    """
    Om = float(ket["Om"])
    P = (dipole_geometry(bra, ket, ctx.I, 0, inner=inner_J)
         - Om * dipole_geometry(bra, ket, ctx.I, 0))
    C = dipole_geometry(bra, ket, ctx.I, 0, inner=inner_flip)
    return 0.5 * ctx.mu_B * (P + flip * C)


_PERP_CITE = (
    "The {c} component of H_Z = mu_B B.G.J with G = diag(G_xx, G_yy, G_zz) in "
    "the molecule frame: B_{a} J_{a} = (1/2)[P {s} C] with "
    "P = B.J - (B.n_hat)(J.n_hat) and C = B_x J_x - B_y J_y. P is B&C Eq. "
    "(9.60) PDF p.638 (lab T1_0(J), inner_J) minus the axial piece the G_zz "
    "term already carries (inner_axial); C is B&C Eq. (9.70) term (vii) and "
    "Eq. (9.71) PDF pp.652-653, evaluated as inner_flip. Both sit under the "
    "SAME lab rank-1 spectator recoupling as stark_z and zeeman_Gzz, B&C "
    "(5.172) + (5.174) + (5.186), because every Zeeman term is linear in B: "
    "Delta J = 0, +-1 and Delta F = 0, +-1. Delta Omega = 0 from P and +-2 "
    "from C. PARITY-EVEN and even in n_hat (quadratic), so no n_hat_sign. "
    "The measured datum constrains only G_zz + G_perp with "
    "G_perp = (G_xx + G_yy)/2, so G_xx and G_yy separately are Petrov & "
    "Skripnikov arXiv:2503.02840 Eqs. (3)-(15) second-order estimates; "
    "G_Delta = (G_xx - G_yy)/2 < 0 is what splits g_e from g_f and it ADDS to "
    "the e level (the component that mixes with the Omega = 0+ states). "
    "docs/superpowers/reports/2026-09-14-zeeman-tensor-handoff-audit.md S2.2, "
    "S3; [HAM] S2.8, OPEN-10.")

_PERP_RULES = Rules(dJ=(-1, 0, 1), dOm=(-2.0, 0.0, 2.0), dF=(-1, 0, 1), dmF=(0,))


@term(name="zeeman_Gxx", param=("G_xx", "B_z"), cases=("c",),
      rules=_PERP_RULES, hermitian=True, real=True,
      cite=_PERP_CITE.format(c="xx", a="x", s="+"))
def zeeman_Gxx(bra, ket, ctx):
    return _zeeman_perp(bra, ket, ctx, +1.0)


@term(name="zeeman_Gyy", param=("G_yy", "B_z"), cases=("c",),
      rules=_PERP_RULES, hermitian=True, real=True,
      cite=_PERP_CITE.format(c="yy", a="y", s="-"))
def zeeman_Gyy(bra, ket, ctx):
    return _zeeman_perp(bra, ket, ctx, -1.0)


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



def _pt_odd_matrix(bra, ket, ctx):
    """Return the n_hat-odd -Omega/|Omega| PT-odd matrix ([HAM] S2.12)."""
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
