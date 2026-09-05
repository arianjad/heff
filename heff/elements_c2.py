"""Case (c) matrix elements for the two-nuclear-spin basis |((J I_Th) F1, I_F) F, m_F>.

EVERY element in this file is copied from docs/thf-plus-x3delta1-effective-
hamiltonian.md section 9 -- the single source for the v2 recoupling -- with its
Brown & Carrington equation number and PDF/book page carried into the decorator,
in exactly the style heff/elements_c.py uses for v1. Nothing here was
re-derived. If a check disagrees with one of these forms, report both with
citations -- do not adjust a sign (Arian's standing rule).

Notation contract, [HAM] S9 preamble, stated once because it is the commonest
way to get a phase backwards:
  * A PRIME MEANS THE BRA. B&C's own convention is the opposite (their primed
    labels are the ket in (5.172)-(5.176), (5.186), (9.50)-(9.53)); every form
    below is the primes-moved-onto-the-bra rewrite that [HAM] S9 prints.
  * Coupling scheme F1 = J + I_Th, F = F1 + I_F, inner-spin-first. I_F is
    ctx.I; I_Th is ctx.spins[0].I (Task 3's Ctx contract).
  * Molecule-frame component q = Om_bra - Om_ket ([HAM] S9 preamble).

Registry: this module registers into REGISTRY_C2, NOT into heff.terms.REGISTRY.
The v1 registry and heff/elements_c.py are untouched, so a v1 session behaves
identically whether or not this module was imported.
"""
import numpy as np

from .conventions import n_hat_sign
from .elements_c import _ph, _same
from .elements_c import (centrifugal as _v1_centrifugal,
                         omega_doubling as _v1_omega_doubling,
                         pt_odd_edm as _v1_pt_odd_edm,
                         pt_odd_scalar_pseudoscalar as _v1_pt_odd_sps,
                         rotation as _v1_rotation)
from .terms import Rules, term
from .wigner import w3j, w6j

REGISTRY_C2 = {}


def _term_c2(**kw):
    """heff.terms.term bound to REGISTRY_C2 and to cases=('c2',)."""
    return term(cases=("c2",), registry=REGISTRY_C2, **kw)


def _spins(ctx):
    """(I_Th, I_F) for the two-spin basis.

    ctx.I is the 19F spin and ctx.spins[0] is the inner (Th) spin -- Task 3's
    contract. A one-spin Ctx has no inner spin, which is a caller error here,
    not a zero: raising beats silently evaluating the whole v2 registry at
    I_Th = 0.
    """
    if not ctx.spins:
        raise ValueError(
            "case 'c2' needs a two-spin Ctx (ctx.spins is empty); build it from a "
            "StateSpec whose `spins` chain has the Th spin coupled to J")
    return float(ctx.spins[0].I), float(ctx.I)


# ------------------------------------------------- the two generic kernels

def axial_geometry(bra, ket, ctx, *, k, q=None, p):
    """Rank-k molecule-frame tensor with lab component p, in the two-spin basis.

    [HAM] S9.1 "The two-spectator axial geometry", the master element, copied
    verbatim (primes = bra):

      <J',Om',F1',F',m'| T^k_p(axial, molecule-frame component q) |J,Om,F1,F,m>
        = (-1)^(F'-m') ( F'  k  F ; -m'  p  m )                                 <- B&C (5.172), PDF p.205 / book p.173
        x (-1)^(F + F1' + k + I_F) sqrt((2F'+1)(2F+1)) { F1  F  I_F ; F' F1' k } <- B&C (5.174), I_F spectator, same pages
        x (-1)^(F1 + J' + k + I_Th) sqrt((2F1'+1)(2F1+1)) { J F1 I_Th ; F1' J' k } <- B&C (5.174) again, I_Th spectator
        x (-1)^(J'-Om') sqrt((2J'+1)(2J+1)) ( J'  k  J ; -Om'  q  Om )          <- B&C (5.186), PDF p.207 / book p.175

    with q = Om' - Om and Delta m_F = p forced by the first 3j. The 6j column
    orders are LITERAL: upper row (ket's inner, ket's total, spectator), lower
    row (bra's total, bra's inner, operator rank). [HAM] S9.1 measures the
    I_F-6j transposition at 0.65 in a quantity of order 1.

    B&C (5.175) is deliberately not used: both reductions have the operator on
    the FIRST constituent of their pair, which is the whole reason the coupling
    scheme is inner-first ([HAM] S9.1, "Why (5.174) twice and (5.175) never").

    At I_Th = 0 line 3 collapses to exactly 1 (a 6j with a zero in the upper
    right, {a a 0; c c f} = (-1)^(a+c+f)/sqrt((2a+1)(2c+1))) and the remaining
    three lines at k=1, q=0 are elements_c.dipole_geometry term for term -- an
    analytic identity, not an approximation ([HAM] S9.1 "Analytic collapse").

    `q` defaults to Om_bra - Om_ket. |q| > k raises: a rank-k operator has no
    such component, so a caller asking for one has a bug, not a zero.
    """
    Om1, Om2 = float(ket["Om"]), float(bra["Om"])
    if q is None:
        q = Om2 - Om1
    if abs(float(q)) > float(k) + 1e-9:
        raise ValueError(
            f"|q| = {abs(float(q))} exceeds the operator rank k = {k}; a rank-k "
            "molecule-frame tensor has no such component")
    I_Th, I_F = _spins(ctx)
    J1, G1, F1, m1 = (float(ket["J"]), float(ket["F1"]), float(ket["F"]),
                      float(ket["mF"]))
    J2, G2, F2, m2 = (float(bra["J"]), float(bra["F1"]), float(bra["F"]),
                      float(bra["mF"]))
    if abs(m2 - m1 - p) > 1e-9:
        return 0.0
    six_F = w6j(G1, F1, I_F, F2, G2, k)
    if six_F == 0.0:
        return 0.0
    six_Th = w6j(J1, G1, I_Th, G2, J2, k)
    if six_Th == 0.0:
        return 0.0
    line1 = _ph(F2 - m2) * w3j(F2, k, F1, -m2, p, m1)
    line2 = (_ph(F1 + G2 + k + I_F)
             * np.sqrt((2 * F2 + 1.0) * (2 * F1 + 1.0)) * six_F)
    line3 = (_ph(G1 + J2 + k + I_Th)
             * np.sqrt((2 * G2 + 1.0) * (2 * G1 + 1.0)) * six_Th)
    line4 = (_ph(J2 - Om2) * np.sqrt((2 * J2 + 1.0) * (2 * J1 + 1.0))
             * w3j(J2, k, J1, -Om2, q, Om1))
    return line1 * line2 * line3 * line4


def outer_spin_scalar(bra, ket, ctx, *, dJ):
    """The recoupled T1(rotational axial) . T1(I_F) kernel. [HAM] S9.3.

    T1(I_F) acts on the OUTER spin while the rotational factor lives inside F1,
    so B&C (5.176) does not apply and F1 goes off-diagonal: Delta F1 = 0, +-1.

    (a) Scalar product across the coupled pair -- B&C Eq. (5.173), PDF p.205 /
    book p.173, with j1 = F1, j2 = I_F, j12 = F, k = 1 and primes on the bra:

      s1 = (-1)^(F1 + F + I_F) { I_F  F1  F ;  F1'  I_F  1 }
           x <J',F1'||T1(A1)||J,F1> x sqrt(I_F(I_F+1)(2I_F+1))

    delta_{F F'} delta_{m m'}. THE PHASE CARRIES THE KET'S F1, not the bra's.
    That distinction is invisible to Hermiticity (both readings give a symmetric
    matrix, max|M - M^T| = 0) and was settled in [HAM] S9.3 by an independent
    rebuild in the fully decoupled |J,m_J>|I_Th,m1>|I_F,m2> basis with explicit
    Clebsch-Gordan coefficients -- a route using neither (5.173) nor (5.174).
    The bra-F1 variant is wrong on every Delta F1 = +-1 element.
    <I_F||T1(I_F)||I_F> = [I_F(I_F+1)(2I_F+1)]^(1/2) is B&C Eq. (5.179), PDF
    p.206 / book p.174 (note: [HAM] S2 miscites this as PDF p.205).

    (b) I_Th spectator -- B&C Eq. (5.174), PDF p.205 / book p.173, with
    j1 = J, j2 = I_Th, j12 = F1, k1 = 1:

      <J',F1'||T1(A1)||J,F1> = (-1)^(F1 + J' + 1 + I_Th) sqrt((2F1'+1)(2F1+1))
                               { J  F1  I_Th ;  F1'  J'  1 } x <J',Om'||T1(A1)||J,Om>

    with the innermost reduced element chosen by `dJ`:
      dJ=False -> A1 = J:   <J'||T1(J)||J> = delta_{JJ'} [J(J+1)(2J+1)]^(1/2)
                            -- B&C (5.179), PDF p.206 / book p.174.
      dJ=True  -> A1 = n^:  <J',Om'||T1(n^)||J,Om> = (-1)^(J'-Om') sqrt((2J'+1)(2J+1))
                            (J' 1 J; -Om' 0 Om) -- B&C (5.186), PDF p.207 /
                            book p.175, i.e. axial_geometry's line 4 at k=1,
                            q=0, with J' != J allowed.

    The two choices are IDENTICAL on the Delta J = 0 block once the A_par/Om
    normalisation is applied ([HAM] S9.3 term 3: the projection theorem gives
    <J,Om||T1(n^)||J,Om> = Om <J||T1(J)||J>/[J(J+1)] exactly), which is why the
    Delta J = +-1 term must exclude Delta J = 0 rather than add to it.
    """
    if not _same(bra, ket, "F", "mF"):
        return 0.0
    I_Th, I_F = _spins(ctx)
    J1, G1, Om1 = float(ket["J"]), float(ket["F1"]), float(ket["Om"])
    J2, G2, Om2 = float(bra["J"]), float(bra["F1"]), float(bra["Om"])
    F = float(ket["F"])
    if dJ:
        inner = (_ph(J2 - Om2) * np.sqrt((2 * J2 + 1.0) * (2 * J1 + 1.0))
                 * w3j(J2, 1, J1, -Om2, 0, Om1))
    else:
        if J1 != J2 or Om1 != Om2:
            return 0.0
        inner = np.sqrt(J1 * (J1 + 1.0) * (2 * J1 + 1.0))
    if inner == 0.0:
        return 0.0
    six_Th = w6j(J1, G1, I_Th, G2, J2, 1)
    if six_Th == 0.0:
        return 0.0
    six_F = w6j(I_F, G1, F, G2, I_F, 1)
    if six_F == 0.0:
        return 0.0
    reduced = (_ph(G1 + J2 + 1.0 + I_Th)
               * np.sqrt((2 * G2 + 1.0) * (2 * G1 + 1.0)) * six_Th * inner)
    return (_ph(G1 + F + I_F) * six_F * reduced
            * np.sqrt(I_F * (I_F + 1.0) * (2 * I_F + 1.0)))


# ------------------------------- the terms that need no recoupling at all

def _delegate(fn):
    """Field adapter: run a v1 element on a KET_C2 row, adding delta_{F1 F1'}.

    [HAM] S9.2: these five operators touch only J and Omega and are diagonal in
    F1, F and m_F, so their matrix elements are the v1 formulas unchanged. The
    v1 functions already read their arguments by field name, and KET_C2 is
    KET_C plus one field, so the only thing the adapter adds is the F1
    diagonality the v1 function cannot know about. One indirection, no new
    algebra (controller ruling: delegate, do not re-implement).
    """
    def wrapped(bra, ket, ctx):
        if bra["F1"] != ket["F1"]:
            return 0.0
        return fn(bra, ket, ctx)
    wrapped.__name__ = fn.__name__
    wrapped.__doc__ = f"Two-spin adapter over elements_c.{fn.__name__}. [HAM] S9.2."
    return wrapped


_DIAG = Rules(dJ=(0,), dOm=(0.0,), dF1=(0,), dF=(0,), dmF=(0,))

rotation = _term_c2(
    name="rotation", param=("B0",), rules=_DIAG, hermitian=True, real=True,
    cite="[HAM] S9.2: H_rot touches only J, so B&C Eq. (5.176) PDF p.205 / book "
         "p.173 (operators acting on the same inner part of a coupled scheme) "
         "makes it diagonal in F1, F and m_F and equal to the v1 element. The "
         "body IS elements_c.rotation through the field adapter -- Ng thesis Eq. "
         "C.7c p.322, H_rot = B_v J(J+1), case (c) with NO -Omega^2 term. "
         "[HAM] S2.1"
)(_delegate(_v1_rotation))

centrifugal = _term_c2(
    name="centrifugal", param=("D0",), rules=_DIAG, hermitian=True, real=True,
    cite="[HAM] S9.2 via B&C Eq. (5.176) PDF p.205 / book p.173: diagonal in F1, "
         "F, m_F and equal to the v1 element, elements_c.centrifugal -- Ng thesis "
         "Eq. C.7d p.322, -D [J(J+1)]^2, the minus in the matrix so D0 stays "
         "positive. [HAM] S2.1"
)(_delegate(_v1_centrifugal))

omega_doubling = _term_c2(
    name="omega_doubling", param=("omega_ef",),
    rules=Rules(dJ=(0,), dOm=(-2.0, 2.0), dF1=(0,), dF=(0,), dmF=(0,)),
    hermitian=True, real=True,
    cite="[HAM] S9.2 via B&C Eq. (5.176) PDF p.205 / book p.173: the Omega-flip "
         "operator acts inside F1, so it is diagonal in F1, F, m_F and equal to "
         "the v1 element, elements_c.omega_doubling -- Ng thesis Eq. C.3 p.319 "
         "transposed to heff's ket phase (-1)^(J-S+s), giving the J-INDEPENDENT "
         "off-diagonal -omega_ef J(J+1)/4. [HAM] S2.3, docs/open-questions.md OQ-A"
)(_delegate(_v1_omega_doubling))

pt_odd_edm = _term_c2(
    name="pt_odd_edm", param=("d_e", "E_eff"), rules=_DIAG,
    hermitian=True, real=True,
    cite="[HAM] S9.2 via B&C Eq. (5.176) PDF p.205 / book p.173: -Omega/|Omega| "
         "is built from the inner part alone, so it is diagonal in F1, F, m_F "
         "and equal to the v1 element, elements_c.pt_odd_edm -- Ng thesis Eq. "
         "C.8 p.322, H_eEDM = -d_e E_eff Omega/|Omega|, no Leanhardt 1/2. "
         "Opt-in: d_e defaults to 0. [HAM] S2.12"
)(_delegate(_v1_pt_odd_edm))

pt_odd_scalar_pseudoscalar = _term_c2(
    name="pt_odd_scalar_pseudoscalar", param=("W_TP", "k_TP"), rules=_DIAG,
    hermitian=True, real=True,
    cite="[HAM] S9.2 via B&C Eq. (5.176) PDF p.205 / book p.173: diagonal in F1, "
         "F, m_F and equal to the v1 element, elements_c.pt_odd_scalar_"
         "pseudoscalar -- Ng Eq. C.8 form with d_e E_eff -> W_TP k_TP "
         "(Skripnikov & Titov 2015 Eqs. 1-4). Opt-in: k_TP defaults to 0. "
         "[HAM] S2.12"
)(_delegate(_v1_pt_odd_sps))


# --------------------------------------------- the recoupled 19F operators

_F_SCALAR = Rules(dJ=(0,), dOm=(0.0,), dF1=(-1, 0, 1), dF=(0,), dmF=(0,))


@_term_c2(name="hyperfine_A_par_F", param=("A_par",), rules=_F_SCALAR,
          hermitian=True, real=True,
          cite="[HAM] S9.3 term 1: B&C Eq. (9.50) PDF p.636 / book p.604 in its "
               "Delta-J = 0 projection form, A_par^F T1(J).T1(I_F) / [J(J+1)] "
               "(the form Ng thesis Eq. C.2 p.319 uses), recoupled through B&C "
               "(5.173) + (5.174) + (5.179). 19F is the OUTER spin, so B&C "
               "(5.176) does NOT apply and F1 is off-diagonal: Delta F1 = 0, "
               "+-1, which is the declared rule and is exactly what [HAM] S9.3 "
               "gives. Cross-checked in [HAM] S9.3 Tables 1 and 2 against a "
               "6j-free double-projection closed form and against [TH] S4.5.")
def hyperfine_A_par_F(bra, ket, ctx):
    J = float(ket["J"])
    if J != float(bra["J"]):
        return 0.0
    return outer_spin_scalar(bra, ket, ctx, dJ=False) / (J * (J + 1.0))


@_term_c2(name="hyperfine_A_par_F_dJ1", param=("A_par",),
          rules=Rules(dJ=(-1, 1), dOm=(0.0,), dF1=(-1, 0, 1), dF=(0,), dmF=(0,)),
          hermitian=True, real=True,
          cite="[HAM] S9.3 term 3: the SAME scalar product with the axial "
               "reduced element, H_hf^F(axial) = (A_par^F / Omega) T1(n^).T1(I_F), "
               "B&C (5.173) + (5.174) + (5.186) PDF pp.205-207 / book pp.173-175. "
               "Its Delta-J = 0 block is identical to hyperfine_A_par_F (the "
               "projection theorem makes the two reduced elements equal there), "
               "so Delta J = 0 is EXCLUDED here to avoid double counting -- "
               "exactly as elements_c.hyperfine_A_par_dJ1 does. At I_Th = 0 it "
               "reproduces B&C Eq. (9.51) PDF p.636 / book p.604 to 5.6e-17 "
               "([HAM] S9.3 self-check). The 1/Omega is signed because B&C's "
               "brace is A_par/Omega. Delta F1 = 0, +-1 as for the Delta-J = 0 "
               "partner.")
def hyperfine_A_par_F_dJ1(bra, ket, ctx):
    if abs(float(bra["J"]) - float(ket["J"])) != 1.0:
        return 0.0
    Om = float(ket["Om"])
    return outer_spin_scalar(bra, ket, ctx, dJ=True) / Om


@_term_c2(name="spin_rotation_cI_F", param=("c_I",), rules=_F_SCALAR,
          hermitian=True, real=True,
          cite="[HAM] S9.3 term 2: H_nsr^F = c_I^F T1(J).T1(I_F), B&C Eq. (8.7) "
               "PDF p.410 / book p.378 -- the same recoupled skeleton as "
               "hyperfine_A_par_F with no 1/[J(J+1)]. Delta F1 = 0, +-1. At "
               "I_Th = 0 it reproduces the (8.20) PDF p.414 / book p.382 closed "
               "form [F(F+1)-I(I+1)-J(J+1)]/2 to 2.2e-16 ([HAM] S9.3 self-check). "
               "c_I is an ESTIMATE (~20 kHz, factor 3). [HAM] S2.6, OPEN-6")
def spin_rotation_cI_F(bra, ket, ctx):
    return outer_spin_scalar(bra, ket, ctx, dJ=False)


@_term_c2(name="zeeman_nuclear_F", param=("g_N", "B_z"),
          rules=Rules(dJ=(0,), dOm=(0.0,), dF1=(-1, 0, 1), dF=(-1, 0, 1),
                      dmF=(0,)),
          hermitian=True, real=True,
          cite="H = -g_N mu_N T1(I_F).T1(B) (Ng thesis Eq. C.6 p.321), the v1 "
               "elements_c.zeeman_nuclear element with the spectator J -> F1. "
               "[HAM] S9 does NOT write this term out (S9.1 notes B&C (5.175) is "
               "unused THERE because both hyperfine reductions put the operator "
               "on the first constituent); the lab-frame nuclear Zeeman is the "
               "one v1 operator that acts on the SECOND constituent of F = F1 + "
               "I_F, so it keeps v1's B&C Eq. (5.175) PDF p.205 / book p.173 "
               "spectator form with the spectator relabelled J -> F1, which is "
               "forced by the S9 coupling contract and by nothing else: "
               "(-1)^(F'+F1+I_F+1) sqrt(I_F(I_F+1)(2I_F+1)) {I_F F' F1; F I_F 1} "
               "x (-1)^(F'-m') sqrt((2F+1)(2F'+1)) (F' 1 F; -m' 0 m), with "
               "<I||T1(I)||I> from B&C (5.179) PDF p.206 / book p.174. "
               "(5.175)'s delta_{j1 j1'} makes the FORMULA diagonal in F1; the "
               "declared Delta F1 = 0, +-1 is deliberately WIDER than that "
               "(brief Step 3: declare wider where the reach is not pinned by "
               "S9), which costs gate A5 only its outside-the-rules half. "
               "g_N mu_N = 4.008 kHz/G. [HAM] S2.8")
def zeeman_nuclear_F(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "F1", "mF"):
        return 0.0
    _, I_F = _spins(ctx)
    G = float(ket["F1"])
    F1, m1 = float(ket["F"]), float(ket["mF"])
    F2, m2 = float(bra["F"]), float(bra["mF"])
    six = w6j(I_F, F2, G, F1, I_F, 1)
    if six == 0.0:
        return 0.0
    a = _ph(F2 + G + I_F + 1.0) * np.sqrt(I_F * (I_F + 1.0) * (2 * I_F + 1.0)) * six
    b = (_ph(F2 - m2) * np.sqrt((2 * F1 + 1.0) * (2 * F2 + 1.0))
         * w3j(F2, 1, F1, -m2, 0, m1))
    return -ctx.mu_N * a * b


# -------------------------------------------------------- the field terms

_FIELD = Rules(dJ=(-1, 0, 1), dOm=(0.0,), dF1=(-1, 0, 1), dF=(-1, 0, 1), dmF=(0,))


@_term_c2(name="stark_z", param=("d_mf", "E_z"), rules=_FIELD,
          hermitian=True, real=True,
          cite="Ng thesis Eq. C.5 p.320, H_Stark = -d_mf n^ . E, evaluated with "
               "the [HAM] S9.1 two-spectator geometry at k = 1, q = 0, p = 0 -- "
               "B&C (5.172) + (5.174) twice + (5.186), PDF pp.205-207 / book "
               "pp.173-175. Identical to the v1 elements_c.stark_z apart from "
               "the added I_Th spectator line, which is exactly 1 at I_Th = 0. "
               "PARITY-ODD: within a J it connects the two parity components of "
               "the Omega-doublet. Odd in n^, so it carries n_hat_sign "
               "([HAM] OPEN-11). Delta F1 = 0, +-1 from the I_Th 6j's triangle. "
               "[HAM] S2.7, S9.1")
def stark_z(bra, ket, ctx):
    if bra["Om"] != ket["Om"]:
        return 0.0
    return -n_hat_sign(ctx.conventions) * axial_geometry(bra, ket, ctx, k=1, q=0, p=0)


@_term_c2(name="zeeman_Gpar", param=("G_par", "B_z"), rules=_FIELD,
          hermitian=True, real=True,
          cite="Ng thesis Eq. C.6 p.321 with the sign CORRECTED to +G_par mu_B "
               "(J.n^)(n^.B) ([HAM] S2.8, OPEN-3), evaluated with the [HAM] S9.1 "
               "two-spectator geometry at k = 1, q = 0, p = 0 -- same four cited "
               "B&C equations as stark_z with -d_mf E_p -> +G_par mu_B B_p "
               "Omega. PARITY-EVEN and EVEN in n^ (quadratic), so unlike the "
               "Stark and PT-odd terms it does NOT carry n_hat_sign. Delta F1 = "
               "0, +-1 from the I_Th 6j's triangle. [HAM] S2.8, S9.1")
def zeeman_Gpar(bra, ket, ctx):
    if bra["Om"] != ket["Om"]:
        return 0.0
    sign = 1.0 if ctx.conventions.zeeman_sign == "plus_Gpar" else -1.0
    return sign * ctx.mu_B * float(ket["Om"]) * axial_geometry(
        bra, ket, ctx, k=1, q=0, p=0)
