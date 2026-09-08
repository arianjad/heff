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
from dataclasses import replace

import numpy as np

from .conventions import n_hat_sign
from .elements_c import _ph, _same
from .elements_c import (centrifugal as _v1_centrifugal,
                         hyperfine_A_par as _v1_hyperfine_A_par,
                         hyperfine_A_par_dJ1 as _v1_hyperfine_A_par_dJ1,
                         omega_doubling as _v1_omega_doubling,
                         pt_odd_edm as _v1_pt_odd_edm,
                         pt_odd_scalar_pseudoscalar as _v1_pt_odd_sps,
                         rotation as _v1_rotation,
                         spin_rotation_cI as _v1_spin_rotation_cI)
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
    # NAMING, so the lines below read against B&C: the locals F1/F2 are the KET
    # and BRA values of the TOTAL F (the ket_c2 field "F"), following B&C's
    # own use of F for the total; G1/G2 are the intermediate F1 = J + I_Th (the
    # field "F1"). The trailing digit is 1 = ket, 2 = bra throughout this
    # module; primes are on the bra ([HAM] S9 convention note).
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
    # Same naming as axial_geometry: G1/G2 are the intermediate F1 = J + I_Th
    # (ket/bra), F is the total F (diagonal here); 1 = ket, 2 = bra.
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
          rules=Rules(dJ=(0,), dOm=(0.0,), dF1=(0,), dF=(-1, 0, 1),
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
               "(5.175)'s delta_{j1 j1'} is exact with j1 = F1 here, so the "
               "declared rule is Delta F1 = 0. g_N mu_N = 4.008 kHz/G. "
               "[HAM] S2.8")
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


# --------------- the Th operators: S9.2's three scalars, then the S9.2.1 Zeeman

def _inner_view(row):
    """A KET_C2 row seen as the one-spin ket |J, Omega, I_Th, F1>.

    A plain dict, because every v1 element reads its ket by field name and
    nothing else: handing it F1 where it looks for F, and (through _inner_spin)
    I_Th where it looks for ctx.I, IS the whole substitution [HAM] S9.2 licenses.
    """
    return {"J": row["J"], "Om": row["Om"], "F": row["F1"], "mF": row["mF"]}


def _inner_spin(fn):
    """Field adapter: run a v1 one-spin element on the INNER pair (J, I_Th, F1).

    The mirror image of _delegate. _delegate keeps the v1 (J, Omega, F, m_F)
    reading and only adds delta_{F1 F1'}; this one REPLACES the v1 (I, F) pair by
    (I_Th, F1) and adds delta_{F F'} delta_{m_F m_F'}.

    That is exactly B&C Eq. (5.176), PDF p.205 / book p.173 ([HAM] S9.2): a
    scalar built from the rotational/electronic degrees of freedom and I_Th alone
    acts on the INNER part j1 = F1 of |((J I_Th) F1, I_F) F, m_F>, so its element
    is diagonal in F and m_F, independent of them and of I_F, and equal to the
    one-spin element evaluated inside |J, Omega, I_Th, F1>. One indirection, and
    NO NEW ALGEBRA is written for any term built through it -- [HAM] S9.2's
    substitution table is the derivation, and it is why S9.3, not S9.2, is where
    the v2 work was.
    """
    def wrapped(bra, ket, ctx):
        if not _same(bra, ket, "F", "mF"):
            return 0.0
        I_Th, _ = _spins(ctx)
        return fn(_inner_view(bra), _inner_view(ket), replace(ctx, I=I_Th))
    wrapped.__name__ = f"{fn.__name__}_Th"
    wrapped.__doc__ = (f"Inner-pair adapter over elements_c.{fn.__name__}: the "
                       "v1 closed form with (I, F) -> (I_Th, F1). [HAM] S9.2.")
    return wrapped


hyperfine_A_par_Th = _term_c2(
    name="hyperfine_A_par_Th", param=("A_par_Th",), rules=_DIAG,
    hermitian=True, real=True,
    cite="[HAM] S9.2 substitution table row 1: B&C Eq. (9.50) PDF p.636 / book "
         "p.604 with I -> I_Th and F -> F1, A_par^Th [F1(F1+1) - J(J+1) - "
         "I_Th(I_Th+1)] / [2 J(J+1)], times delta_{FF'} delta_{mm'}. Exact, not "
         "approximate: B&C Eq. (5.176) PDF p.205 / book p.173 makes a scalar "
         "built from the rotational degrees of freedom and I_Th alone diagonal "
         "in F and m_F and independent of I_F. The body IS elements_c."
         "hyperfine_A_par through the inner-pair adapter -- no new algebra. "
         "Coefficients checked against [TH] S4.1 at J = 1-4 (gate V18) and "
         "reproduced in [HAM] S9.2's own cross-check run. The A_par_Th Param is "
         "SIGNED; the sign selects the trusted calculation via params.thf_v2"
         "(..., a_par_th_sign=...), default 'negative' ([HAM] S9.4.3, "
         "OPEN-16). [HAM] S9.2"
)(_inner_spin(_v1_hyperfine_A_par))

hyperfine_A_par_Th_dJ1 = _term_c2(
    name="hyperfine_A_par_Th_dJ1", param=("A_par_Th",),
    rules=Rules(dJ=(-1, 1), dOm=(0.0,), dF1=(0,), dF=(0,), dmF=(0,)),
    hermitian=True, real=True,
    cite="[HAM] S9.2 substitution table row 2: B&C Eq. (9.51) PDF p.636 / book "
         "p.604 with I -> I_Th and F -> F1, J the LARGER of the two and the "
         "brace = A_par^Th / Omega; diagonal in F1, F, m_F by B&C Eq. (5.176) "
         "PDF p.205 / book p.173. The body IS elements_c.hyperfine_A_par_dJ1 "
         "through the inner-pair adapter -- no new algebra -- and that function "
         "already excludes Delta J = 0 itself, so there is no double counting "
         "against hyperfine_A_par_Th. NOT a small correction for Th: [TH] S4.2 "
         "puts the J = 1 <-> 2 element at ~2 GHz with a ~135 MHz second-order "
         "shift, five orders above the 19F analogue's ~2.6 kHz ([HAM] S2.5), so "
         "the J truncation has to be re-tested. The A_par_Th Param is SIGNED; "
         "the sign selects the trusted calculation via params.thf_v2(..., "
         "a_par_th_sign=...), default 'negative' ([HAM] S9.4.3, OPEN-16). "
         "[HAM] S9.2"
)(_inner_spin(_v1_hyperfine_A_par_dJ1))

spin_rotation_cI_Th = _term_c2(
    name="spin_rotation_cI_Th", param=("c_I_Th",), rules=_DIAG,
    hermitian=True, real=True,
    cite="[HAM] S9.2 substitution table row 3: B&C Eq. (8.7) PDF p.410 / book "
         "p.378, H_nsr = c_I T1(J).T1(I), with its coupled-basis element Eq. "
         "(8.20) PDF p.414 / book p.382 read with I -> I_Th and F -> F1, "
         "c_I^Th [F1(F1+1) - I_Th(I_Th+1) - J(J+1)] / 2; diagonal in F1, F, m_F "
         "by B&C Eq. (5.176) PDF p.205 / book p.173. The body IS elements_c."
         "spin_rotation_cI through the inner-pair adapter -- no new algebra. "
         "c_I_Th has NO value in any source: [TH] S4.6 declines to pick one and "
         "brackets it at ~1 kHz to ~1 MHz, so the Param is held at 0 (gap G3, "
         "OPEN-19). [HAM] S9.2"
)(_inner_spin(_v1_spin_rotation_cI))


@_term_c2(name="zeeman_nuclear_Th", param=("g_N_Th", "B_z"),
          rules=Rules(dJ=(0,), dOm=(0.0,), dF1=(-1, 0, 1), dF=(-1, 0, 1),
                      dmF=(0,)),
          hermitian=True, real=True,
          cite="H = -g_N^Th mu_N T1(I_Th).T1(B), the ONE Th operator B&C Eq. "
               "(5.176) does not cover: it is a LAB-frame rank-1 operator on "
               "the inner spin, not a scalar, so [HAM] S9.2's substitution rule "
               "does not apply and the element is a two-step recoupling written "
               "out in [HAM] S9.2.1 -- B&C Eq. (5.172) PDF p.205 / book p.173 "
               "in F and m_F; Eq. (5.174) same pages with I_F a spectator "
               "(identical to axial_geometry's line 2 at k = 1); Eq. (5.175) "
               "same pages with I_Th the SECOND constituent of F1 = J + I_Th "
               "and J the spectator (zeeman_nuclear_F's recoupler with "
               "(F1, I_F) -> (J, I_Th), i.e. elements_c.zeeman_nuclear's "
               "structure); and Eq. (5.179) PDF p.206 / book p.174 for "
               "<I_Th||T1(I_Th)||I_Th> = [I(I+1)(2I+1)]^(1/2). Same sign "
               "convention as zeeman_nuclear_F, -g_N mu_N. Delta F1 = 0, +-1 "
               "and Delta F = 0, +-1 from the two rank-1 6j triangles. Pinned "
               "element by element against a decoupled-basis rebuild in "
               "tests/test_elements_c2_th.py (controller ruling R13); the "
               "bra/ket-F1 phase swap is invisible to Hermiticity and is caught "
               "there. [HAM] S9.2.1, S2.8")
def zeeman_nuclear_Th(bra, ket, ctx):
    if not _same(bra, ket, "J", "Om", "mF"):
        return 0.0
    I_Th, I_F = _spins(ctx)
    J = float(ket["J"])
    G1, F1, m1 = float(ket["F1"]), float(ket["F"]), float(ket["mF"])
    G2, F2, m2 = float(bra["F1"]), float(bra["F"]), float(bra["mF"])
    six_F = w6j(G1, F1, I_F, F2, G2, 1)
    if six_F == 0.0:
        return 0.0
    six_Th = w6j(I_Th, G2, J, G1, I_Th, 1)
    if six_Th == 0.0:
        return 0.0
    line1 = _ph(F2 - m2) * w3j(F2, 1, F1, -m2, 0, m1)
    line2 = (_ph(F1 + G2 + 1.0 + I_F)
             * np.sqrt((2 * F2 + 1.0) * (2 * F1 + 1.0)) * six_F)
    line3 = (_ph(G2 + J + I_Th + 1.0)
             * np.sqrt((2 * G2 + 1.0) * (2 * G1 + 1.0)) * six_Th
             * np.sqrt(I_Th * (I_Th + 1.0) * (2 * I_Th + 1.0)))
    return -ctx.mu_N * line1 * line2 * line3


# ------------------------------------------------------ the Th quadrupole

# The explicit -1 that encodes conventions.quadrupole_convention =
# 'bc_q0_is_negative_efg'. B&C (9.52) prints the prefactor -(1/2) eQ
# <T2_q(grad E)>; their constant is defined by "q0 is the negative of the
# electric field gradient", i.e. eq_qQ = -2 eQ <T2_q(grad E)> ([HAM] S9.4.1,
# read off by comparing (9.52) at q = 0 with (9.53)). Read with the OPPOSITE
# convention (q0 = +EFG) the element would be -eq_qQ/4; B&C's convention turns
# that into +eq_qQ/4, and this factor is that turn, written out instead of
# folded into a +1/4 literal. Drop it and every quadrupole splitting flips
# sign: gate V20's uniform ratio of -1 against the textbook Casimir function
# becomes +1 ([HAM] S9.4.2, [TH] S3.4).
_Q0_IS_NEGATIVE_EFG = -1.0

# The smallest I that HAS a quadrupole moment: <I||T2(Q)||I> needs the triangle
# (I, 2, I). Named rather than inlined so gate V21's FAIL demonstration can
# monkeypatch it to 0.0 and run the REAL product path unguarded (which is where
# the 0/0 lives), instead of a hand copy of the formula.
_MIN_I_FOR_QUADRUPOLE = 1.0


def _q0_sign(ctx):
    """conventions.quadrupole_convention -> the sign factor above.

    Mirrors quadrupole_eQq2_Th's eqq2_norm guard: the constant is a convention
    field's consequence, so it is read from the Ctx rather than hardcoded, and
    an unrecognised fork raises instead of silently keeping B&C's. _ALLOWED
    currently lists one value, so the raise is a tripwire for the day a second
    q0 convention is added, not a branch reachable through a valid Conventions.
    """
    conv = ctx.conventions.quadrupole_convention
    if conv != "bc_q0_is_negative_efg":
        raise NotImplementedError(
            f"quadrupole_convention={conv!r} is not implemented: the only "
            "transcribed convention is B&C's 'bc_q0_is_negative_efg' (q0 is the "
            "NEGATIVE of the electric field gradient, [HAM] S9.4.1). A new fork "
            "needs its own sign derived from its own printed equations, not a "
            "guess.")
    return _Q0_IS_NEGATIVE_EFG


def _quadrupole_body(bra, ket, ctx, q):
    """B&C Eq. (9.52) at molecule-frame component q, in units of the constant.

    [HAM] S9.4.1, copied with the primes moved onto the bra so that q =
    Om_bra - Om_ket as everywhere else in S9, and with I -> I_Th, F -> F1 (S9.2,
    B&C Eq. (5.176)):

      <J',Om',F1',F',m'| H_Q |J,Om,F1,F,m>
        = (e q_q Q / 4) (-1)^(J + I_Th + F1' + J' - Om')                <- (9.53) phase
          x sqrt((2J'+1)(2J+1)) { J  I_Th  F1' ;  I_Th  J'  2 }
          x ( J'  2  J ;  -Om'  q  Om )
          x ( I_Th  2  I_Th ;  -I_Th  0  I_Th )^(-1)
          x delta_{F1 F1'} delta_{F F'} delta_{m m'}

    (primes = bra throughout, so B&C's own unprimed labels are the bra's here).
    Returns the geometry in units of eq_qQ; the caller's `param` supplies the
    constant. Delta J = 0, +-1, +-2 and q = Om_bra - Om_ket, both read straight
    off the (J' 2 J) 3j -- B&C say so themselves in the sentence under (9.53).

    THE I < 1 GUARD IS PHYSICS, NOT DEFENSIVE CODING ([HAM] S9.4.1, [TH] S1.4).
    A nucleus with I < 1 has no quadrupole moment at all: the moment is the
    rank-2 element <I||T2(Q)||I>, which needs the triangle (I, 2, I), so eQ = 0
    and H_Q = 0 identically. The (I 2 I; -I 0 I) in the DENOMINATOR vanishes on
    that same triangle, so the printed expression at I <= 1/2 is 0/0 -- singular,
    not zero -- and a literal transcription divides by zero. Reading the
    vanishing denominator as the reason for the vanishing element is backwards.
    So 227ThF+ (I_Th = 1/2) has no Th quadrupole, exactly as 19F has none
    ([HAM] S2.11), and 232ThF+ has none because I_Th = 0. Gate V21.
    """
    if not _same(bra, ket, "F1", "F", "mF"):
        return 0.0
    I_Th, _ = _spins(ctx)
    if I_Th < _MIN_I_FOR_QUADRUPOLE:
        return 0.0
    J1, Om1 = float(ket["J"]), float(ket["Om"])
    J2, Om2 = float(bra["J"]), float(bra["Om"])
    if abs(Om2 - Om1 - float(q)) > 1e-9:
        return 0.0
    G = float(ket["F1"])
    den = w3j(I_Th, 2, I_Th, -I_Th, 0, I_Th)
    return (_q0_sign(ctx) * -0.25
            * _ph(J1 + I_Th + G + J2 - Om2)
            * np.sqrt((2 * J2 + 1.0) * (2 * J1 + 1.0))
            * w6j(J1, I_Th, G, I_Th, J2, 2)
            * w3j(J2, 2, J1, -Om2, float(q), Om1) / den)


@_term_c2(name="quadrupole_eQq0_Th", param=("eQq0_Th",),
          rules=Rules(dJ=(-2, -1, 0, 1, 2), dOm=(0.0,), dF1=(0,), dF=(0,),
                      dmF=(0,)),
          hermitian=True, real=True,
          cite="B&C Eq. (9.53) PDF p.637 / book p.605, the q = 0 specialisation "
               "of (9.52), with I -> I_Th and F -> F1 ([HAM] S9.2, B&C Eq. "
               "(5.176) PDF p.205 / book p.173) -- see _quadrupole_body for the "
               "transcription. Delta J = 0, +-1, +-2 (B&C's own sentence under "
               "(9.53)), Delta Omega = 0, Delta F1 = Delta F = Delta m_F = 0. "
               "The explicit -1 of conventions.quadrupole_convention = "
               "'bc_q0_is_negative_efg' is in _quadrupole_body: at Omega = 0, "
               "J' = J this element is MINUS the textbook Casimir function, "
               "uniformly, which is B&C's stated q0 sign ([HAM] S9.4.2, [TH] "
               "S3.4, gate V20). Third source for the same rank-2 case-(c) "
               "skeleton: Skripnikov, Petrov, Titov & Flambaum arXiv:1408.5368 "
               "Eq. (5). eQq0_Th has no published ThF+ value at all -- the "
               "-2600 MHz Param is an uncalibrated sensitivity placeholder ([TH] S4.3, "
               "gap G2). [HAM] S9.4.1")
def quadrupole_eQq0_Th(bra, ket, ctx):
    return _quadrupole_body(bra, ket, ctx, 0.0)


@_term_c2(name="quadrupole_eQq2_Th", param=("eQq2_Th",),
          rules=Rules(dJ=(-2, -1, 0, 1, 2), dOm=(-2.0, 2.0), dF1=(0,), dF=(0,),
                      dmF=(0,)),
          hermitian=True, real=True,
          cite="B&C Eq. (9.52) PDF p.636 / book p.604 at q = Om_bra - Om_ket = "
               "-+2, with I -> I_Th and F -> F1 ([HAM] S9.2, B&C Eq. (5.176)). "
               "Inside the Omega = +-1 block of 3Delta1 this is the Omega = +1 "
               "<-> Omega = -1 element: parity-EVEN, Delta F1 = Delta F = "
               "Delta m_F = 0, Delta J = 0, +-1, +-2, sitting in the same matrix "
               "position as omega_doubling with a different (J, F1) law -- which "
               "is what gate V22 discriminates. At [TH] S4.4's ~200-400 MHz it "
               "would be 13-30x the Omega-doubling operator and would set the "
               "parity-doublet structure of 229ThF+, the behaviour Petrov 2018 "
               "SV reports in 177HfF+. "
               "NORMALISATION, OPEN-17: this computes in B&C (9.52)'s own "
               "q = +-2 normalisation, where the constant is eq2Q = -2 eQ "
               "<eta,Lam|T2_(+-2)(grad E)|eta,Lam'> -- [HAM] S9.4.4's [derived] "
               "extension of B&C's own eq0Q definition, since B&C print no "
               "eq2Q. That is conventions.eqq2_norm = 'bc_9p52_q2', the only "
               "implemented value. The bridge to Petrov 2018 Eq. (23) is NOT "
               "resolved: an unresolved sqrt(2) (whether Petrov's sqrt(2 pi/5) "
               "is intended or is a typo for sqrt(4 pi/5) = C^2_q) and an "
               "unresolved sign leave two candidate factors, -1/sqrt(3) and "
               "-1/sqrt(6), that the printed equations do not discriminate, so "
               "'petrov2018_eq23' raises rather than guessing a branch. What IS "
               "pinned is the RELATIVE bridge: a Petrov-style eQq0 entered into "
               "(9.53) as-is requires eQq2/sqrt(6) entered here. [HAM] S9.4.4")
def quadrupole_eQq2_Th(bra, ket, ctx):
    norm = ctx.conventions.eqq2_norm
    if norm != "bc_9p52_q2":
        raise NotImplementedError(
            f"eqq2_norm={norm!r} is not implemented: OPEN-17, the eQq2 "
            "normalisation bridge between B&C Eq. (9.52) at q = +-2 and Petrov "
            "2018 PRA 98 042502 Eq. (23), is NOT resolved. An unresolved "
            "sqrt(2) (is the sqrt(2 pi/5) of Petrov's Eqs. (19), (22), (23) "
            "intended, or a typo for sqrt(4 pi/5) = Racah's C^2_q?) and an "
            "unresolved sign leave two candidate factors, -1/sqrt(3) and "
            "-1/sqrt(6); the printed equations do not discriminate and "
            "[HAM] S9.4.4's probe of Petrov's own Eq. (24) prefactor was "
            "inconclusive. A guessed factor is a stop-work condition. Use "
            "eqq2_norm='bc_9p52_q2' and read eQq2_Th as B&C's eq2Q, or supply "
            "the answer to either question in [HAM] S9.4.4's OPEN-17.")
    q = float(bra["Om"]) - float(ket["Om"])
    # |q| = 2 ONLY: at Delta Omega = 0 the same body is the eQq0 element, and
    # returning it here would double-count it under the wrong constant (and put
    # non-zero elements outside this term's declared Delta Omega = +-2, which is
    # what gate A5 measures).
    if abs(q) != 2.0:
        return 0.0
    return _quadrupole_body(bra, ket, ctx, q)


# --------------------------------------------------- J-truncation reporting

def j_convergence(isotopologue, *, J_maxes=(2, 4, 6), mF=None, n_levels=8, knobs=None):
    """Report how the lowest `n_levels` levels move as J_max is truncated.

    A REPORTING helper (task-6 brief): it never raises on a magnitude, only
    on a structural mismatch (an m_F the requested block does not hold). Each
    J_max's own lowest `n_levels` eigenvalues are compared to the largest
    J_max's (the "converged" reference), and the largest absolute difference,
    in MHz, is what is reported -- this is what [HAM] S2.5/S9.2's Th
    Delta-J = +-1 hyperfine (a ~2 GHz off-diagonal element, [TH] S4.2) needs a
    convergence table for, and V23 (tests/test_assemble_c2.py) is a direct
    demonstration of this function, not a parallel computation.

    Field-free by default (knobs=None): all field knobs (E_z, B_z) are absent
    from thf_v2's own ParamSet, so `hamiltonian`'s knob lookup defaults them
    to 0 -- J-truncation is a Hamiltonian-structure question, not a Stark/
    Zeeman one. Pass `knobs` to report the same convergence at a field point.

    `mF` has no single value that is a level of all three isotopologues:
    F (and so m_F) is INTEGER for 229ThF+ and 227ThF+ (I_Th half-integer + 19F
    1/2 -> F1 half-integer -> F integer, [HAM] S9's coupling scheme) and
    HALF-INTEGER for 232ThF+ (the v1 F = J +- 1/2 basis) -- verified
    numerically, not assumed. The default is therefore per-isotopologue:
    0.0 for '229'/'227', 0.5 for '232'.

    '232' has no coupled Th spin (thf_spec('232').spins == ()), so it builds
    on the v1 KET_C basis through the default (v1) registry -- gate V16
    (tests/test_elements_c2_reduction.py) already proves every REGISTRY_C2
    term equals its v1 twin at I_Th = 0, so this is not a second code path,
    just the one thf_spec('232') already hands back. '229' and '227' build on
    KET_C2 through REGISTRY_C2.
    """
    from .assemble import build_term_matrices, hamiltonian
    from .params import thf_v2
    from .spec import block_by_mF, enumerate_kets, thf_spec
    from .terms import ctx_from

    if mF is None:
        mF = 0.0 if isotopologue in ("229", "227") else 0.5
    knobs = dict(knobs) if knobs else {}
    pset = thf_v2(isotopologue)
    two_spin = isotopologue in ("229", "227")

    lowest, dims = {}, {}
    for J_max in J_maxes:
        spec = thf_spec(isotopologue, J_max=J_max)
        kets = enumerate_kets(spec)
        blocks = block_by_mF(kets)
        if mF not in blocks.index:
            raise ValueError(
                f"mF={mF} is not a level of {isotopologue}ThF+ at J_max={J_max} "
                f"(available: {sorted(blocks.index)})")
        idx = blocks.index[mF]
        ctx = ctx_from(spec, pset)
        if two_spin:
            tm = build_term_matrices(kets[idx], ctx, case="c2", registry=REGISTRY_C2)
        else:
            tm = build_term_matrices(kets[idx], ctx)
        w = np.linalg.eigvalsh(hamiltonian(tm, pset, knobs))
        lowest[J_max] = np.sort(w)[:n_levels]
        dims[J_max] = len(idx)

    ref_J = max(J_maxes)
    ref = lowest[ref_J]
    by_J_max = {}
    for J_max in J_maxes:
        n = min(len(lowest[J_max]), len(ref))
        shift = float(np.max(np.abs(lowest[J_max][:n] - ref[:n]))) if n else float("nan")
        by_J_max[J_max] = {"dimension": dims[J_max], "max_shift_MHz": shift,
                           "n_levels_compared": n}
    return {"isotopologue": isotopologue, "mF": mF, "n_levels": n_levels,
            "reference_J_max": ref_J, "field_free": not knobs, "by_J_max": by_J_max}
