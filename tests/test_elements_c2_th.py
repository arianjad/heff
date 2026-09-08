"""Th hyperfine, Zeeman, and quadrupole gates at I_Th != 0.

V16 covers the I_Th=0 reduction; independent checks here cover its blind spots.
See [HAM] S9.2/S9.4 and [TH] S1.4/S3.4/S4.
"""
import numpy as np
import pytest

from heff.assemble import build_term_matrices, coefficients
from heff.conventions import Conventions, parity_operator
from heff.elements_c import _ph
from heff.elements_c2 import REGISTRY_C2, _quadrupole_body
from heff.params import MU_B, MU_N, thf_v2
from heff.spec import KET_C2, Spin, StateSpec, enumerate_kets, thf_spec
from heff.terms import Ctx, check_selection_rules, ctx_from, terms_for_case
from heff.wigner import w3j, w6j

TH_TERMS = ("hyperfine_A_par_Th", "hyperfine_A_par_Th_dJ1", "spin_rotation_cI_Th",
            "zeeman_nuclear_Th", "quadrupole_eQq0_Th", "quadrupole_eQq2_Th")


# ------------------------------------------------------------------ helpers

def row(J, Om, F1, F, mF):
    """One KET_C2 record, for the kernel-level gates that need a (J, I) pair
    the ThF+ basis does not contain (V20 runs at Omega = 0 and I = 1, 3/2)."""
    return np.array([(J, Om, F1, F, mF)], dtype=KET_C2)[0]


def ctx_at(I_Th, *, I_F=0.5, conventions=None):
    """A two-spin Ctx at an arbitrary inner spin, built directly rather than
    through a StateSpec -- V20's (J, I) pairs are textbook cases, not ThF+."""
    return Ctx(S=1.0, Lam=2.0, I=I_F, mu_B=MU_B, mu_N=MU_N,
               conventions=conventions or Conventions(), frame="rotating",
               spins=(Spin(label="X", I=I_Th, couple_to="J"),
                      Spin(label="19F", I=I_F, couple_to="F1")))


def setup_229(J_max=2):
    spec = thf_spec("229", J_max=J_max)
    return enumerate_kets(spec), ctx_from(spec, thf_v2("229")), spec


def dense(fn, kets, ctx):
    """Every (i, j) with no selection-rule mask -- a corrupted element must show
    up as a number, not as a masked zero."""
    d = len(kets)
    M = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            M[i, j] = float(np.real(fn(kets[i], kets[j], ctx)))
    return M


def diagonal_over_F1(name, J, kets, ctx):
    """The diagonal element of one term at each F1 of a given J, on the
    Omega = +1, F = F1 + 1/2, m_F = F component (the term is independent of F
    and m_F by B&C (5.176), which test_th_terms_are_independent_of_F pins)."""
    fn = REGISTRY_C2[name].fn
    sel = (kets["J"] == J) & (kets["Om"] == 1.0)
    out = {}
    for F1 in sorted({float(v) for v in kets["F1"][sel]}):
        hit = np.flatnonzero(sel & (kets["F1"] == F1) & (kets["F"] == F1 + 0.5))
        i = hit[np.argmax(kets["mF"][hit])]
        out[F1] = float(fn(kets[i], kets[i], ctx))
    return out


# --------------------------------------------------------- V18: B&C (9.50)

# [TH] S4.1, the diagonal (9.50) manifold in units of A_par^Th at I_Th = 5/2.
TH_HYPERFINE_MANIFOLD = {
    1: (-1.7500, -0.5000, +1.2500),
    2: (-1.1667, -0.9167, -0.5000, +0.0833, +0.8333),
    3: (-0.8333, -0.7083, -0.5000, -0.2083, +0.1667, +0.6250),
    4: (-0.6250, -0.5000, -0.3250, -0.1000, +0.1750, +0.5000),
}


def test_th_hyperfine_manifold_coefficients_match_the_digest():
    """V18. B&C Eq. (9.50) with I -> I_Th, F -> F1 ([HAM] S9.2), over every F1
    at J = 1-4, against [TH] S4.1 (and reproduced in [HAM] S9.2's cross-check
    run for J = 1).

    Uniquely catches mis-coupling I_Th to F instead of to J. That error leaves a
    plausible manifold -- same number of levels, same order of magnitude -- with
    the wrong spacings, and no symmetry gate sees it: it is Hermitian, parity-
    even and traceless just like the right answer.

    The digest tabulates the coefficients in units of A_par^Th (R14: the
    element function itself is now unsigned -- the SIGNED physical constant
    lives entirely in the A_par_Th Param, applied at assembly, not here -- so
    no multiplier needs to be divided out to compare against the digest).

    TOLERANCE. [TH] S4.1 prints four decimals, so the comparison is ABSOLUTE to
    half the last printed digit (abs = 5e-5), the same convention V19 already
    uses against the same digest. The brief's rel = 1e-4 would reject the J = 2,
    F1 = 7/2 entry, where the exact coefficient is 1/12 = 0.0833333 and the
    table rounds it to 0.0833 -- a rounding artefact of the table, not a
    disagreement, and 4e-4 relative on the smallest entry is still four orders
    tighter than any error this gate exists to catch.
    """
    kets, ctx, _ = setup_229(J_max=4)
    for J, expect in TH_HYPERFINE_MANIFOLD.items():
        got = diagonal_over_F1("hyperfine_A_par_Th", J, kets, ctx)
        assert len(got) == len(expect), f"J={J}: F1 count {sorted(got)}"
        for (F1, value), want in zip(sorted(got.items()), expect):
            assert value == pytest.approx(want, abs=5e-5), f"J={J}, F1={F1}"


def test_th_hyperfine_is_linear_in_the_signed_A_par_Th_param():
    """R14: A_par(Th)'s sign is no longer a conventions.py fork; it lives in
    the SIGNED A_par_Th Param, selected by params.thf_v2(..., a_par_th_sign=).
    The raw element functions are unsigned (previous test), so the physical
    sign only shows up once the term matrix is scaled by the Param's value
    (heff.assemble.coefficients) -- exactly like any other Param, and unlike
    n_hat_sign/zeeman_sign, which the element itself still multiplies in.
    """
    kets, ctx, spec = setup_229(J_max=2)
    tm = build_term_matrices(kets, ctx, case="c2", registry=REGISTRY_C2)
    pset_neg = thf_v2("229")
    pset_pos = thf_v2("229", a_par_th_sign="positive")
    assert pset_neg.value("A_par_Th") == pytest.approx(-1510)
    assert pset_pos.value("A_par_Th") == pytest.approx(1510)
    c_neg = coefficients(tm, pset_neg, {})
    c_pos = coefficients(tm, pset_pos, {})
    for k, name in enumerate(tm.names):
        M_neg = c_neg[k] * tm.mats[k]
        M_pos = c_pos[k] * tm.mats[k]
        if name in ("hyperfine_A_par_Th", "hyperfine_A_par_Th_dJ1"):
            assert np.max(np.abs(M_neg)) > 0.0, name
            assert np.max(np.abs(M_neg + M_pos)) < 1e-12, (
                f"{name} does not flip sign with a_par_th_sign")
        else:
            assert np.max(np.abs(M_neg - M_pos)) < 1e-12, (
                f"{name} changed between a_par_th_sign values; only A_par_Th "
                f"should")


def test_227_A_par_Th_diagonal_carries_the_placeholder_sign_at_J1():
    """R14: 227ThF+'s A_par_Th = +39821 MHz is a SIGNED physical value ([HAM]
    S9.6), not a magnitude paired with a_par_th_sign, so the actual assembled
    J = 1 diagonal element must show the sign of a POSITIVE A_par: from the
    closed form A_par[F1(F1+1) - J(J+1) - I(I+1)] / [2 J(J+1)] at J = 1,
    I_Th = 1/2, the F1 = 3/2 coefficient is +0.25 and the F1 = 1/2 coefficient
    is -0.5, so with A_par_Th > 0 the LARGER F1 (3/2) is positive and the
    smaller F1 (1/2) is negative.
    """
    spec = thf_spec("227", J_max=1)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("227"))
    ps = thf_v2("227")
    a_par = ps.value("A_par_Th")
    assert a_par == pytest.approx(39821, abs=1)
    got = diagonal_over_F1("hyperfine_A_par_Th", 1, kets, ctx)
    signed = {F1: coeff * a_par for F1, coeff in got.items()}
    assert signed[1.5] > 0.0, signed
    assert signed[0.5] < 0.0, signed


def test_th_terms_are_independent_of_F_and_m_F():
    """B&C (5.176) as a testable statement ([HAM] S9.2): a Th scalar's element
    depends on F1 and not on F, m_F or I_F. Without this, V18's choice of the
    F = F1 + 1/2 row would be arbitrary."""
    kets, ctx, _ = setup_229(J_max=2)
    for name in ("hyperfine_A_par_Th", "spin_rotation_cI_Th", "quadrupole_eQq0_Th"):
        fn = REGISTRY_C2[name].fn
        seen = {}
        for i in range(len(kets)):
            key = (float(kets["J"][i]), float(kets["Om"][i]), float(kets["F1"][i]))
            val = float(fn(kets[i], kets[i], ctx))
            if key in seen:
                assert val == pytest.approx(seen[key], abs=1e-12), f"{name} {key}"
            seen[key] = val


# --------------------------------------------------- V20: the Casimir ratio

def casimir(J, I, F):
    """The textbook quadrupole function, [TH] S3.4 / [HAM] S9.4.2."""
    C = F * (F + 1.0) - I * (I + 1.0) - J * (J + 1.0)
    return (0.75 * C * (C + 1.0) - I * (I + 1.0) * J * (J + 1.0)) / (
        2.0 * I * (2.0 * I - 1.0) * (2.0 * J - 1.0) * (2.0 * J + 3.0))


CASIMIR_CASES = ((1, 1.0), (2, 1.0), (2, 1.5), (3, 2.5))


def casimir_ratios():
    """(J, I, F, body, casimir) for every case of [HAM] S9.4.2's printed table,
    skipping the single 0/0 row it prints as `nan`."""
    out = []
    for J, I in CASIMIR_CASES:
        ctx = ctx_at(I)
        for twoF in range(int(round(2 * abs(J - I))), int(round(2 * (J + I))) + 1, 2):
            F = twoF / 2.0
            k = row(J, 0.0, F, F + 0.5, 0.5)
            out.append((J, I, F, _quadrupole_body(k, k, ctx, 0.0), casimir(J, I, F)))
    return out


def test_quadrupole_reproduces_the_casimir_function_with_a_ratio_of_minus_one():
    """V20. At Omega = 0, J' = J, B&C Eq. (9.53) equals MINUS the textbook
    Casimir function, uniformly, for every (J, I, F) of [HAM] S9.4.2's table
    ([TH] S3.4 reports the same ratio from an independent implementation).

    That -1 IS B&C's stated convention that q0 is the negative of the electric
    field gradient. Uniquely catches the convention being dropped: nothing else
    in the package sees it, because a global sign flip of one term is Hermitian,
    parity-even and invisible to every reduction gate (all of them run at
    I_Th = 0 where the quadrupole vanishes identically).
    """
    rows = casimir_ratios()
    assert len(rows) == 16, "S9.4.2 prints sixteen rows"
    checked = 0
    for J, I, F, body, cas in rows:
        if abs(cas) < 1e-12:
            assert abs(body) < 1e-12, f"({J},{I}) F={F}: 0/0 row is not 0/0"
            continue
        assert body / cas == pytest.approx(-1.0, abs=1e-12), f"({J},{I}) F={F}"
        checked += 1
    assert checked == 15, "fifteen non-degenerate rows, one 0/0 row"


# ------------------- V20b: (9.52) transcribed independently, off-diagonal J

def bc_9p52(bra, ket, ctx, q=0.0):
    """Independent B&C 9.52/9.53 transcription (PDF pp.636–637; book pp.604–605).

    B&C primes are kets, so this uses its unprimed bra and I→I_Th, F→F1
    substitutions ([HAM] S9.2); B&C q0 is negative EFG.
    """
    I = ctx.spins[0].I
    if I < 1.0:                          # no rank-2 moment; the printed 0/0
        return 0.0
    if not all(float(bra[f]) == float(ket[f]) for f in ("F1", "F", "mF")):
        return 0.0
    J_bra, Om_bra = float(bra["J"]), float(bra["Om"])
    J_ket, Om_ket = float(ket["J"]), float(ket["Om"])
    F1 = float(ket["F1"])
    return (0.25
            * _ph(J_ket + I + F1 + J_bra - Om_bra)
            * np.sqrt((2 * J_bra + 1.0) * (2 * J_ket + 1.0))
            * w6j(J_ket, I, F1, I, J_bra, 2)
            * w3j(J_bra, 2, J_ket, -Om_bra, q, Om_ket)
            / w3j(I, 2, I, -I, 0, I))


def _eQq0_matrices(J_max=3):
    """(kets, |eQq0_Th| geometry from the registry, the same from bc_9p52)."""
    spec = thf_spec("229", J_max=J_max)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("229"))
    got = dense(REGISTRY_C2["quadrupole_eQq0_Th"].fn, kets, ctx)
    want = dense(bc_9p52, kets, ctx)
    return kets, got, want


def _dJ_classes(kets):
    return np.abs(kets["J"][:, None] - kets["J"][None, :])


def test_V20b_quadrupole_dJ_elements_match_an_independent_transcription_of_bc_9p52():
    """V20b: compare all 360-ket eQq0 elements with independent B&C 9.52/9.53.

    Unlike V20's diagonal Casimir limit, this catches J-dependent phase and 3j
    order errors; the rebuild shares only wigner and ``_ph``.
    """
    kets, got, want = _eQq0_matrices()
    assert len(kets) == 360
    dev = float(np.max(np.abs(got - want)))
    assert dev < 1e-12, f"max |registry - independent (9.52)| = {dev:.3e}"

    dJ = _dJ_classes(kets)
    reach = {}
    for k in (0, 1, 2):
        cls = np.abs(got[dJ == k])
        reach[k] = float(np.max(cls))
        assert np.sum(cls > 1e-12) > 0, f"Delta J = {k} class is empty"
    # Off-diagonal-J reach is measured rather than assumed.
    assert reach[1] > reach[0], f"Delta J reach {reach}"


# ------------------------------------------- V21: no quadrupole below I = 1

def test_quadrupole_is_identically_zero_for_I_Th_at_most_one_half():
    """V21. Both quadrupole matrices are exactly 0.0 at I_Th = 1/2 and I_Th = 0
    -- FROM THE FORMULA's own guard, not from a parameter being zero.

    The physics is the numerator, not the denominator ([HAM] S9.4.1, [TH] S1.4):
    a nucleus with I < 1 has no quadrupole moment, because <I||T2(Q)||I> needs
    the triangle (I, 2, I). The (I 2 I; -I 0 I) in the denominator of (9.52)
    vanishes on that same triangle, so the printed expression at I <= 1/2 is
    0/0. Uniquely catches a code path that manufactures a 227ThF+ quadrupole:
    227Th's params carry NO eQq0_Th/eQq2_Th at all (structurally absent, not
    zero), so an unguarded formula would be caught by nothing downstream -- the
    assembler would default its knob to 0.0 and hide a nan.
    """
    for iso, I_Th in (("227", 0.5), ("232", 0.0)):
        spec = thf_spec("227", J_max=2)
        spec = StateSpec(case="c", electronic=spec.electronic, I=0.5,
                         J_range=spec.J_range, v=0, M="blocks", frame="rotating",
                         spins=(Spin(label=f"{iso}Th", I=I_Th, couple_to="J"),
                                Spin(label="19F", I=0.5, couple_to="F1")))
        kets = enumerate_kets(spec)
        ctx = ctx_from(spec, thf_v2(iso))
        for name in ("quadrupole_eQq0_Th", "quadrupole_eQq2_Th"):
            M = dense(REGISTRY_C2[name].fn, kets, ctx)
            assert np.all(M == 0.0), f"I_Th = {I_Th}: {name} has a non-zero element"


# ------------------------------------------------------------- V22: eQq2

# [TH] S4.4, <J, Om=-1, F1| H_Q |J, Om=+1, F1> in units of eq2Q, I_Th = 5/2.
EQQ2_J1 = (+0.1715, -0.1960, +0.0612)


def eqq2_J1_coefficients(kets, ctx):
    fn = REGISTRY_C2["quadrupole_eQq2_Th"].fn
    out = {}
    sel = kets["J"] == 1.0
    for F1 in (1.5, 2.5, 3.5):
        F, mF = F1 + 0.5, F1 + 0.5
        base = sel & (kets["F1"] == F1) & (kets["F"] == F) & (kets["mF"] == mF)
        i = np.flatnonzero(base & (kets["Om"] == -1.0))[0]
        j = np.flatnonzero(base & (kets["Om"] == +1.0))[0]
        out[F1] = float(fn(kets[i], kets[j], ctx))
    return out


def test_eQq2_J1_coefficients_match_the_digest():
    """V22. The Omega = +1 <-> Omega = -1 element at J = 1, over F1 = 3/2, 5/2,
    7/2, against [TH] S4.4: +0.1715, -0.1960, +0.0612 in units of eq2Q.

    NORMALISATION. Those are B&C (9.52)'s own q = +-2 normalisation, which is
    what the package computes and what conventions.eqq2_norm='bc_9p52_q2'
    names. [HAM] S9.4.4 could NOT pin the bridge to Petrov 2018 Eq. (23) -- an
    unresolved sqrt(2) and an unresolved sign (OPEN-17) -- so no rescaling is
    applied here and none is implied; the digest's numbers are reproduced as
    printed, with no derived factor.

    Uniquely catches confusing eQq2 with omega_doubling. The two occupy the
    SAME matrix position (Delta Omega = +-2, diagonal in J, F1, F, m_F) and
    differ only in their (J, F1) dependence, so a term registered with the
    wrong law is invisible to parity, hermiticity and A5 alike.
    """
    kets, ctx, _ = setup_229(J_max=1)
    got = eqq2_J1_coefficients(kets, ctx)
    for F1, want in zip((1.5, 2.5, 3.5), EQQ2_J1):
        assert got[F1] == pytest.approx(want, rel=1e-3), f"F1={F1}"


def test_eQq2_is_dOmega_two_diagonal_in_J_F1_F_mF_and_parity_even():
    """V22, the structural half: every non-zero eQq2 element has |Delta Omega|
    = 2 and Delta F1 = Delta F = Delta m_F = 0, and the term is parity-even.

    CORRECTION to the name, which follows [TH] S3.2's looser phrase "diagonal in
    J, F1, F, m_F": the term is NOT diagonal in J. B&C's own sentence under
    (9.53) and [TH] S3.5 both give Delta J = 0, +-1, +-2 from the 3j
    (J' 2 J; -Om' q Om), and [HAM] S9.4.1 reads the same two implications
    straight off it. The Delta J != 0 elements are asserted non-zero here so the
    correction is pinned by data rather than by a docstring.
    """
    kets, ctx, _ = setup_229(J_max=2)
    M = dense(REGISTRY_C2["quadrupole_eQq2_Th"].fn, kets, ctx)
    nz = np.argwhere(np.abs(M) > 1e-12)
    assert len(nz) > 0
    dJ = set()
    for i, j in nz:
        assert abs(kets["Om"][i] - kets["Om"][j]) == 2.0, "Delta Omega != 2"
        for f in ("F1", "F", "mF"):
            assert kets[f][i] == kets[f][j], f"off-diagonal in {f}"
        dJ.add(float(kets["J"][i] - kets["J"][j]))
    # >= not ==: {0, +-1} is all a J = 1-2 basis can show, but the rule (and
    # B&C's sentence under (9.53)) is Delta J = 0, +-1, +-2, so a wider basis
    # must be allowed to populate +-2 without failing this assertion.
    assert dJ >= {-1.0, 0.0, 1.0} and dJ <= {-2.0, -1.0, 0.0, 1.0, 2.0}, (
        f"Delta J reach should be 0, +-1 (and +-2 once J = 3 is in the basis), "
        f"got {sorted(dJ)}")
    P = parity_operator(kets, 1.0, ell=0, s=0)
    assert np.max(np.abs(P @ M @ P.T - M)) < 1e-12, "eQq2 is not parity-even"


def test_eQq2_refuses_the_petrov_normalisation_naming_OPEN_17():
    """[HAM] S9.4.4: Petrov 2018 Eq. (23) leaves the normalization unresolved.

    The candidate factors, -1/sqrt(3) and -1/sqrt(6), are indistinguishable
    from the printed equations, so the converter raises instead of guessing.
    """
    kets, _, _ = setup_229(J_max=1)
    ctx = ctx_at(2.5, conventions=Conventions(eqq2_norm="petrov2018_eq23"))
    i = int(np.flatnonzero((kets["J"] == 1.0) & (kets["Om"] == -1.0)
                           & (kets["F1"] == 1.5))[0])
    j = int(np.flatnonzero((kets["J"] == 1.0) & (kets["Om"] == +1.0)
                           & (kets["F1"] == 1.5) & (kets["F"] == kets["F"][i])
                           & (kets["mF"] == kets["mF"][i]))[0])
    with pytest.raises(NotImplementedError) as e:
        REGISTRY_C2["quadrupole_eQq2_Th"].fn(kets[i], kets[j], ctx)
    msg = str(e.value)
    assert "OPEN-17" in msg and "sqrt(2)" in msg and "sign" in msg


# ------------------------------------------- the Delta J = +-1 magnitude gate

def test_th_hyperfine_dJ1_is_hundreds_of_MHz_not_kHz():
    """[TH] S4.2: at |A_par^Th| = 1510 MHz the J = 1 <-> 2, F1 = 7/2 element is
    ~2 GHz and its second-order shift |ME|^2/(2 B J_upper) is ~135 MHz, five
    orders of magnitude above the 19F analogue's ~2.6 kHz ([HAM] S2.5).

    Not a precision gate -- a MAGNITUDE gate at rel = 0.1. What it catches is
    the term being registered against `A_par` (the 19F constant, 20.1 MHz)
    instead of `A_par_Th`, which is a one-character error that every angular
    gate in this file passes: the matrix is identical, only the knob is wrong.
    """
    assert REGISTRY_C2["hyperfine_A_par_Th_dJ1"].param == ("A_par_Th",)
    assert REGISTRY_C2["hyperfine_A_par_F_dJ1"].param == ("A_par",)
    kets, ctx, _ = setup_229(J_max=2)
    pset = thf_v2("229")
    sel = (kets["Om"] == 1.0) & (kets["F1"] == 3.5) & (kets["F"] == 4.0) \
        & (kets["mF"] == 4.0)
    i = int(np.flatnonzero(sel & (kets["J"] == 2.0))[0])
    j = int(np.flatnonzero(sel & (kets["J"] == 1.0))[0])
    geo = float(REGISTRY_C2["hyperfine_A_par_Th_dJ1"].fn(kets[i], kets[j], ctx))
    element = abs(pset.value("A_par_Th") * geo)
    shift = element ** 2 / (2.0 * pset.value("B0") * 2.0)
    assert element == pytest.approx(1980.0, rel=0.1), f"element = {element:.1f} MHz"
    assert shift == pytest.approx(135.0, rel=0.1), f"shift = {shift:.1f} MHz"

    f_geo = float(REGISTRY_C2["hyperfine_A_par_F_dJ1"].fn(kets[i], kets[j], ctx))
    f_element = abs(pset.value("A_par") * f_geo)
    f_shift = f_element ** 2 / (2.0 * pset.value("B0") * 2.0)
    assert shift / f_shift > 1e4, (
        f"Th shift {shift:.1f} MHz vs 19F {f_shift * 1e3:.2f} kHz")


# ---------------------------------------- R13: the Th nuclear Zeeman rebuild

def clebsch_gordan(j1, m1, j2, m2, J3, M3):
    """<j1 m1 j2 m2 | J3 M3> from the 3j -- B&C Eq. (5.83), used ONLY by the
    independent rebuild below, so the rebuild shares no machinery with the
    recoupled formula beyond w3j itself."""
    if abs(m1 + m2 - M3) > 1e-9:
        return 0.0
    return _ph(j1 - j2 + M3) * np.sqrt(2 * J3 + 1.0) * w3j(j1, j2, J3, m1, m2, -M3)


def coupled_and_product(J, I_Th, I_F):
    """The coupled |((J I_Th) F1, I_F) F, m_F> list, the fully decoupled
    |J m_J>|I_Th m1>|I_F m2> list, and the CG transformation between them."""
    coupled = []
    for twoG in range(int(round(2 * abs(J - I_Th))), int(round(2 * (J + I_Th))) + 1, 2):
        G = twoG / 2.0
        for twoF in range(int(round(2 * abs(G - I_F))), int(round(2 * (G + I_F))) + 1, 2):
            F = twoF / 2.0
            for twom in range(int(round(-2 * F)), int(round(2 * F)) + 1, 2):
                coupled.append((G, F, twom / 2.0))
    product = [(mJ, m1, m2) for mJ in np.arange(-J, J + 1)
               for m1 in np.arange(-I_Th, I_Th + 1)
               for m2 in np.arange(-I_F, I_F + 1)]
    U = np.zeros((len(coupled), len(product)))
    for a, (G, F, mF) in enumerate(coupled):
        for b, (mJ, m1, m2) in enumerate(product):
            U[a, b] = (clebsch_gordan(J, mJ, I_Th, m1, G, mJ + m1)
                       * clebsch_gordan(G, mJ + m1, I_F, m2, F, mF))
    return coupled, product, U


def test_R13_th_nuclear_zeeman_matches_a_decoupled_basis_rebuild():
    """The Th nuclear Zeeman is the ONE lab-frame Th operator, so B&C (5.176)
    and [HAM] S9.2's substitution rule do not cover it ([HAM] S9.2.1). Its
    element is a two-step recoupling -- (5.172), then (5.174) with I_F a
    spectator, then (5.175) with I_Th the second constituent of F1, then (5.179)
    -- and none of that is checked by any other gate here: it vanishes at
    I_Th = 0 (so V16 cannot see it) and it is the only Th term off-diagonal in
    F1 and F.

    So it is pinned the way [HAM] S9.3 pinned the (5.173) phase: rebuild
    T1_0(I_Th) = I_Th,z in the fully decoupled |J m_J>|I_Th m1>|I_F m2> basis,
    transform with two layers of Clebsch-Gordan coefficients, and compare EVERY
    element. The rebuild uses neither (5.174) nor (5.175).

    zeeman_nuclear_F is checked on the same fixture at one extra assert -- it
    agrees analytically with (5.175), and this numerical pin is
    free here.
    """
    J, I_Th, I_F = 1.0, 2.5, 0.5
    coupled, product, U = coupled_and_product(J, I_Th, I_F)
    assert len(coupled) == len(product) == 36
    assert np.max(np.abs(U @ U.T - np.eye(len(coupled)))) < 1e-12, "CG not unitary"
    ctx = ctx_at(I_Th, I_F=I_F)
    kets = np.array([(J, 1.0, G, F, mF) for G, F, mF in coupled], dtype=KET_C2)

    for name, which, fn in (
            ("zeeman_nuclear_Th", 1, REGISTRY_C2["zeeman_nuclear_Th"].fn),
            ("zeeman_nuclear_F", 2, REGISTRY_C2["zeeman_nuclear_F"].fn)):
        D = np.diag([p[which] for p in product])
        want = -ctx.mu_N * (U @ D @ U.T)
        got = dense(fn, kets, ctx)
        dev = float(np.max(np.abs(got - want)))
        assert dev < 1e-12, f"{name}: max |formula - rebuild| = {dev:.3e}"



# ------------------------------------------------- registry-wide obligations

def test_every_th_term_carries_a_citation():
    """Every element in this package names its primary source in the decorator
    (the elements_c/elements_c2 module contract). For the four terms that are v1
    closed forms read through the inner-pair adapter, the cite must also say so:
    no new algebra was written for them ([HAM] S9.2)."""
    for name in TH_TERMS:
        cite = REGISTRY_C2[name].cite
        assert len(cite) > 40 and "B&C" in cite and "[HAM]" in cite, name
    for name in ("hyperfine_A_par_Th", "hyperfine_A_par_Th_dJ1",
                 "spin_rotation_cI_Th"):
        assert "no new algebra" in REGISTRY_C2[name].cite, name
    assert "OPEN-17" in REGISTRY_C2["quadrupole_eQq2_Th"].cite
    assert "bc_q0_is_negative_efg" in REGISTRY_C2["quadrupole_eQq0_Th"].cite


def test_A5_holds_for_every_th_term():
    """Gate A5 on the real 229ThF+ basis (thf_spec('229', J_max=2), 192 kets):
    every term is non-zero SOMEWHERE inside its declared rules and EXACTLY zero
    outside them.

    Run over the whole of REGISTRY_C2, not just the six Th terms: the eleven
    The term-level A5 check at I_Th = 0 does not exercise these terms, and a rule
    set can be right at I_Th = 0 and wrong at I_Th = 5/2 (the I_Th 6j opens
    Delta F1 = +-1 that does not exist at I_Th = 0).
    """
    kets, ctx, _ = setup_229(J_max=2)
    for t in terms_for_case("c2", registry=REGISTRY_C2):
        r = check_selection_rules(t, kets, ctx)
        assert r["n_nonzero_outside"] == 0, (
            f"{t.name}: {r['n_nonzero_outside']} non-zero elements outside its "
            f"declared rules, max {r['max_outside']:.3e}")
        assert r["n_nonzero_inside"] > 0, f"{t.name}: dead operator"


def test_every_th_term_is_hermitian_on_the_229_basis():
    """Declared hermitian=True, so the matrix must be symmetric on the real
    basis. Free, and it is the check build_term_matrices would raise on."""
    kets, ctx, _ = setup_229(J_max=2)
    for name in TH_TERMS:
        M = dense(REGISTRY_C2[name].fn, kets, ctx)
        assert np.max(np.abs(M - M.T)) < 1e-12, name


def test_every_th_term_vanishes_at_I_Th_zero():
    """Why V16 (tests/test_elements_c2_reduction.py) remains unchanged:
    every Th term is identically zero on the 96-ket I_Th = 0 basis the master
    gate runs on -- the magnetic ones because their Casimir bracket and their
    <I||T1(I)||I> both vanish, the quadrupole ones by the I < 1 guard."""
    v1 = thf_spec(J_max=4)
    spec = StateSpec(case="c", electronic=v1.electronic, I=0.5,
                     J_range=v1.J_range, v=0, M="blocks", frame="rotating",
                     spins=(Spin(label="232Th", I=0.0, couple_to="J"),
                            Spin(label="19F", I=0.5, couple_to="F1")))
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v2("232"))
    assert len(kets) == 96
    for name in TH_TERMS:
        M = dense(REGISTRY_C2[name].fn, kets, ctx)
        assert np.all(M == 0.0), f"{name} is non-zero at I_Th = 0"
