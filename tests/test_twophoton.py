"""Task 7 gates: the rank-K effective two-photon (2 x E1) operator within X 3Delta1.

Sources, in the order the gates use them:
  [HAM] = docs/thf-plus-x3delta1-effective-hamiltonian.md S9.5 (S9.5.1 the rank
          decomposition and the closure form, S9.5.2 which (K, dOmega) channels
          exist, S9.5.3 the OPEN-21 ruling on K = 1, S9.5.4 the selection rules
          as data, S9.5.5 validity)
  [2g]  = docs/digest-literature-two-photon.md S0 item 6, S3.0, S3.3
  B&C   = Brown & Carrington Eqs. (5.141), (5.142), PDF p.198 / book p.166

Two things these gates pin that nothing else in the suite can:

  * The operator is a TRANSITION operator living in its own registry
    (REGISTRY_2G), so it can never be summed into an assembled Hamiltonian
    ([SPEC-v2] S3.2). test_no_two_photon_term_appears_in_a_hamiltonian_registry
    is the check.
  * The registered rank set is K in {0, 2}. [HAM] S9.5.3 rules that K = 1 is
    identically zero in exact closure (it is the antisymmetric part of the
    dyad, and the Cartesian components of d commute), so there is no
    alpha_K1_* parameter and no K = 1 operator. `dyad_weights` still RETURNS
    the K = 1 weights -- they are a property of the polarisation pair, not of
    the molecule, and the sin^2(theta) law Cossel measured lives in them (see
    test_dyad_weights_reproduce_the_known_polarisation_limits).
"""
import numpy as np
import pytest

from heff.conventions import parity_operator
from heff.elements_c2 import REGISTRY_C2, axial_geometry
from heff.params import thf_v2
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import REGISTRY, check_selection_rules, terms_for_case, ctx_from
from heff.twophoton import (REGISTRY_2G, dyad_weights, two_photon_geometry,
                            two_photon_line_strengths, two_photon_matrix)

SQ2 = np.sqrt(2.0)
# Jones vectors in the [HAM] S2 Condon-Shortley phase: sigma+ is
# e_{+1} = -(x + iy)/sqrt(2), the vector S9.5.1(3) anchors the whole
# convention on ("this gives c[+1] = +1 and nothing else").
SIGMA_P = np.array([-1.0, -1.0j, 0.0]) / SQ2
SIGMA_M = np.array([1.0, -1.0j, 0.0]) / SQ2
PI_Z = np.array([0.0, 0.0, 1.0])

CHANNELS = ((0, 0), (2, 0), (2, 2))       # (K, |dOmega|), [HAM] S9.5.2


def setup_229(J_max=2):
    spec = thf_spec("229", J_max=J_max)
    return enumerate_kets(spec), ctx_from(spec, thf_v2("229"))


@pytest.fixture(scope="module")
def basis2():
    """229ThF+ at J_max = 2 (192 kets), the whole basis (every m_F).

    Merged rather than m_F-blocked because a two-photon operator with P != 0
    connects DIFFERENT m_F blocks, and because conventions.parity_operator
    needs both signed-Omega partners of every ket present.
    """
    return setup_229(J_max=2)


@pytest.fixture(scope="module")
def mats(basis2):
    """{(K, dOmega, P): matrix} on the J_max = 2 basis, P = -3 .. +3.

    Built once (about 3 s) and shared by V24 and V25; P = +-3 is included
    because "never +-3" is the claim V25(b) has to be able to falsify.
    """
    kets, ctx = basis2
    return {(K, dOm, P): two_photon_matrix(kets, kets, ctx, K=K, dOmega=dOm, P=P)
            for (K, dOm) in CHANNELS for P in range(-3, 4)}


# ------------------------------------------------------------------ V24

def test_two_photon_operator_commutes_with_parity_at_zero_field(basis2, mats):
    """V24: [T, P] = 0 for every registered (K, dOmega) channel and every P.

    [2g] S0 item 6 / S3.3, `[derived]` there: P d P+ = -d, and at zero field
    parity commutes with H, so the parity partner |i'> = P|i> carries the same
    detuning and P T_eff P+ = (-1)^2 T_eff. [2g] measured max|[T,P]| <= 5.6e-17
    with max|T| ~ 0.6 on a resolved-intermediate probe; this gate re-runs the
    check on the rank-K operator itself rather than citing that number.

    The dOmega = +-2 channel is ONE channel: the q = +2 and q = -2 components
    share one alpha (alpha_K2_dOm2) and parity maps each into the other, so
    splitting them would make this identity false by construction -- which is
    exactly what test_..._fails_if_the_q_plus_2_component_changes_sign shows.
    """
    kets, ctx = basis2
    Par = parity_operator(kets, ctx.S, ell=0.0, s=0.0)
    worst, biggest = 0.0, 0.0
    for (K, dOm) in CHANNELS:
        for P in range(-K, K + 1):
            T = mats[(K, dOm, P)]
            assert np.max(np.abs(T)) > 0.1, f"K={K} dOm={dOm} P={P}: dead channel"
            biggest = max(biggest, float(np.max(np.abs(T))))
            worst = max(worst, float(np.max(np.abs(T @ Par - Par @ T))))
    assert worst < 1e-12, f"max|[T,P]| = {worst:.3e} (max|T| = {biggest:.3f})"
    assert biggest > 0.1                      # non-vacuous: the operator is real


def test_V24_fails_if_the_q_plus_2_component_changes_sign(basis2):
    """FAIL demo for V24, through the real code path.

    Every element is the true two_photon_geometry output; only the SIGN of the
    q = +2 half of the dOmega = +-2 channel is flipped. That is the corruption
    V24 uniquely catches: magnitudes stay right, hermiticity survives, and the
    only symptom is that the operator now connects e to f -- a fake two-photon
    parity flip ([2g] S3.3: "the notebook should say this rather than advertise
    a two-photon parity flip").
    """
    kets, ctx = basis2
    Par = parity_operator(kets, ctx.S, ell=0.0, s=0.0)
    n = len(kets)
    T = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            q = float(kets["Om"][i]) - float(kets["Om"][j])
            if abs(q) != 2.0:
                continue
            T[i, j] = np.sign(-q) * two_photon_geometry(kets[i], kets[j], ctx,
                                                        K=2, P=0)
    assert np.max(np.abs(T)) > 0.1                       # non-vacuous
    assert np.max(np.abs(T @ Par - Par @ T)) > 0.1, (
        "flipping the q = +2 sign should have broken [T, P] = 0")


# ------------------------------------------------------------------ V25

def test_channel_reach(basis2, mats):
    """V25, three parts, each with its own negative control ([HAM] S9.5.2/S9.5.4).

    (a) dOmega = +-2 is a K = 2 channel only. The molecule-frame 3j
        (J' K J; -Om' q Om) vanishes identically unless |q| <= K, so K = 0
        cannot reach |q| = 2 -- and there is no K = 1 operator to ask about
        ([HAM] S9.5.3, OPEN-21), which this gate asserts directly.
    (b) Delta m_F = P with |P| <= K <= 2, so the reachable set over the whole
        basis is {0, +-1, +-2} and NEVER +-3.
    (c) With sigma+- polarisations only (the JILA constraint, Ng thesis p.102
        fn.4) the dyad narrows the reach to {0, +-2}.
    """
    kets, ctx = basis2

    # (a) -- the dOmega = +-2 matrix at K = 0 is identically zero ...
    K0_dOm2 = two_photon_matrix(kets, kets, ctx, K=0, dOmega=2, P=0)
    assert np.max(np.abs(K0_dOm2)) == 0.0
    # ... and non-zero at K = 2 (negative control: the channel is not just
    # empty for want of connected states)
    assert np.max(np.abs(mats[(2, 2, 0)])) > 0.1
    # the leak this closure guards against is a hard error one level down
    with pytest.raises(ValueError):
        axial_geometry(kets[0], kets[0], ctx, k=0, q=2.0, p=0)
    # OPEN-21: no K = 1 operator exists at all
    assert not [n for n in REGISTRY_2G if "K1" in n], sorted(REGISTRY_2G)
    assert all(t.param[0] in ("alpha_K0_dOm0", "alpha_K2_dOm0", "alpha_K2_dOm2")
               for t in REGISTRY_2G.values())

    # (b)
    reach = {P for (K, dOm, P), M in mats.items() if np.max(np.abs(M)) > 1e-12}
    assert reach == {0, -1, 1, -2, 2}, sorted(reach)
    assert all(np.max(np.abs(mats[(K, dOm, P)])) == 0.0
               for (K, dOm) in CHANNELS for P in (-3, 3))

    # (c) sigma+- only; the negative control is the sigma/pi pair, which does
    # reach +-1
    sigma_reach, pi_reach = set(), set()
    for e1 in (SIGMA_P, SIGMA_M):
        for e2 in (SIGMA_P, SIGMA_M):
            sigma_reach |= {P for (K, P), w in dyad_weights(e1, e2).items()
                            if abs(w) > 1e-12}
        pi_reach |= {P for (K, P), w in dyad_weights(e1, PI_Z).items()
                     if abs(w) > 1e-12}
    assert sigma_reach == {0, -2, 2}, sorted(sigma_reach)
    assert pi_reach == {-1, 1}, sorted(pi_reach)


# ------------------------------------------------------------------ V26

def sum_rule(kets, ctx, i, K, dOm):
    """Sum over every final basis state and every P of |<f|alpha^K_P|i>|^2."""
    return sum(float(np.sum(two_photon_matrix(kets[i:i + 1], kets, ctx, K=K,
                                              dOmega=dOm, P=P) ** 2))
               for P in range(-K, K + 1))


def test_rank_K_sum_rule_and_reciprocity():
    """V26, the B7/B7b analogues at rank K.

    Sum rule: the total strength out of a state, summed over final states and
    over P, is independent of m_F -- sum_P T^K+_P T^K_P is a rotational scalar,
    so its expectation cannot depend on the orientation of the initial state.
    In heff's normalisation (conventions.two_photon_norm = 'bc_5p142_reduced',
    unit one-photon reduced elements) the closed form is exactly 1.0 per
    channel, by 3j/6j completeness -- but only while the basis holds every J'
    the operator can reach, |dJ| <= K = 2, hence J_max = 3 and J = 1 initial
    states here.

    Reciprocity: the geometry is real, so the Wigner-Eckart relation
    T^K+_P = (-1)^P T^K_{-P} collapses to M_P(a, b) = (-1)^P M_{-P}(b, a)^T --
    the rank-K form of gate B7b. THIS is the phase check: the sum rule above
    sums squares and no unimodular phase can move it (v1's own correction to
    B7's docstring), so a missing (-1)^(F'-m'_F) is invisible there and visible
    here. The FAIL demo is the next test.
    """
    kets3, ctx3 = setup_229(J_max=3)
    for K, dOm in CHANNELS:
        for (J, Om, F1, F) in ((1, 1.0, 3.5, 4.0), (1, 1.0, 2.5, 2.0)):
            vals = []
            for mF in np.arange(-F, F + 0.5, 1.0):
                i = int(np.flatnonzero((kets3["J"] == J) & (kets3["Om"] == Om)
                                       & (kets3["F1"] == F1) & (kets3["F"] == F)
                                       & (kets3["mF"] == mF))[0])
                vals.append(sum_rule(kets3, ctx3, i, K, dOm))
            assert np.allclose(vals, vals[0], rtol=1e-10), f"K={K} dOm={dOm}: {vals}"
            assert vals[0] == pytest.approx(1.0, rel=1e-10), (
                f"K={K} dOm={dOm} (J={J}, F1={F1}, F={F}): {vals[0]}")

    kets, ctx = setup_229(J_max=2)
    blocks = block_by_mF(kets)
    for P, mF_b in ((-1, 1.0), (-2, 2.0)):
        a, b = kets[blocks.index[0.0]], kets[blocks.index[mF_b]]
        for K, dOm in ((2, 0), (2, 2)):
            lhs = two_photon_matrix(a, b, ctx, K=K, dOmega=dOm, P=P)
            rhs = (-1.0) ** P * two_photon_matrix(b, a, ctx, K=K, dOmega=dOm,
                                                  P=-P).T
            assert np.max(np.abs(lhs)) > 0.0, f"K={K} dOm={dOm} P={P}: vacuous"
            assert np.allclose(lhs, rhs), f"K={K} dOm={dOm} P={P}: reciprocity broken"


def test_V26_fails_if_the_wigner_eckart_phase_is_dropped():
    """FAIL demo for V26's reciprocity half, exactly as v1's
    test_B7b_fails_if_the_wigner_eckart_phase_is_dropped does it: multiply the
    real element by its own (-1)^(F'-m'_F) (the BRA's F and m_F -- primes are
    the bra throughout S9), which cancels that phase because it is +-1.

    Built by calling heff.elements_c2.axial_geometry directly rather than by
    monkeypatching, so the demo does not depend on twophoton.py's import style.
    """
    kets, ctx = setup_229(J_max=2)
    blocks = block_by_mF(kets)
    a, b = kets[blocks.index[0.0]], kets[blocks.index[1.0]]

    def corrupted(ka, kb, P):
        out = np.zeros((len(ka), len(kb)))
        for i in range(len(ka)):
            for j in range(len(kb)):
                q = float(ka["Om"][i]) - float(kb["Om"][j])
                if abs(q) > 2.0:
                    continue
                out[i, j] = (axial_geometry(ka[i], kb[j], ctx, k=2, q=q, p=P)
                             * (-1.0) ** (ka["F"][i] - ka["mF"][i]))
        return out

    lhs = corrupted(a, b, -1)
    rhs = (-1.0) ** (-1) * corrupted(b, a, +1).T
    assert np.max(np.abs(lhs)) > 0.0                      # non-vacuous
    assert not np.allclose(lhs, rhs), (
        "dropping the (-1)^(F'-m') phase should have broken reciprocity")


# ------------------------------------------------------------------ V28

def test_amplitudes_are_summed_over_channels_then_squared():
    """V28: the two-photon form of v1's gate B8 and of Cossel's measured
    sigma-pathway cancellation ([2g] S3.0; his Fig. 6.18, p.224: "the two-sigma
    pathways ... cancel because of the signs of the Wigner 3j coefficients").

    Real code path: two channels (K = 0 and K = 2, both dOmega = 0, P = 0) on
    the real 229ThF+ basis, with alpha_K2_dOm0 tuned so that the two amplitudes
    into one chosen line are equal and opposite. Summing first gives zero;
    squaring each channel first gives twice one channel's strength.
    """
    kets, ctx = setup_229(J_max=2)
    a = kets[block_by_mF(kets).index[0.0]]
    ev = np.zeros(len(a)), np.eye(len(a))             # bare basis states

    def strengths(alphas):
        return two_photon_line_strengths(ev[0], ev[1], a, ev[0], ev[1], a, ctx,
                                         eps1=PI_Z, eps2=PI_Z, alphas=alphas)[1]

    one = {"alpha_K0_dOm0": 1.0, "alpha_K2_dOm0": 0.0, "alpha_K2_dOm2": 0.0}
    two = {"alpha_K0_dOm0": 0.0, "alpha_K2_dOm0": 1.0, "alpha_K2_dOm2": 0.0}
    s0, s2 = strengths(one), strengths(two)
    i, j = np.unravel_index(int(np.argmax(np.minimum(s0, s2))), s0.shape)
    assert s0[i, j] > 1e-6 and s2[i, j] > 1e-6        # both channels are live

    # amplitudes are real here (real geometry, real pi-pi dyad), so the ratio
    # of the strengths gives the alpha that makes them cancel, up to a sign
    t = np.sqrt(s0[i, j] / s2[i, j])
    both = [strengths({"alpha_K0_dOm0": 1.0, "alpha_K2_dOm0": s * t,
                       "alpha_K2_dOm2": 0.0})[i, j] for s in (-1.0, 1.0)]
    assert min(both) == pytest.approx(0.0, abs=1e-18), both
    assert max(both) == pytest.approx(4.0 * s0[i, j], rel=1e-9)
    # square-then-sum would have given this instead, and never zero
    assert s0[i, j] + t ** 2 * s2[i, j] == pytest.approx(2.0 * s0[i, j])


# ------------------------------------------------------------------- A5

def test_A5_holds_for_every_registered_two_photon_channel():
    """Gate A5 on the merged 229ThF+ basis at J_max = 3 (360 kets): every
    registered channel is non-zero SOMEWHERE inside its declared Rules and
    EXACTLY zero outside them.

    J_max = 3, not 2: |dJ| <= K = 2 is part of the declared rule set, and a
    J = 1, 2 basis holds no dJ = +-2 pair at all, so the gate would be blind to
    a wrong dJ rule there.
    """
    kets, ctx = setup_229(J_max=3)
    assert len(REGISTRY_2G) == 3
    for t in terms_for_case("c2", registry=REGISTRY_2G):
        r = check_selection_rules(t, kets, ctx)
        assert r["n_nonzero_outside"] == 0, (
            f"{t.name}: {r['n_nonzero_outside']} non-zero elements outside its "
            f"declared rules, max {r['max_outside']:.3e}")
        assert r["n_nonzero_inside"] > 0, f"{t.name}: dead operator"


# ------------------------------------------------------------- the dyad

def test_dyad_weights_reproduce_the_known_polarisation_limits():
    """The polarisation dyad of [HAM] S9.5.1(3), against limits derived from
    its own printed formula.

    THE READING IS RAMAN, ep_2 CONJUGATED: T_eff = (d.eps2*)|i><i|(d.eps1)/D,
    so slot 1 (bra-side) carries eps2* and slot 2 (ket-side) carries eps1, and
    c[p] = (-1)^p eps_{-p}. Conjugating eps2 flips the sign of its helicity
    label, so the beam-pair -> Delta m_F map is INVERTED relative to the ladder
    reading where both photons are absorbed. Derived from S9.5.1(3) and
    reproduced verbatim in its table (rows d):

        eps1 = sigma+, eps2 = sigma+  ->  legs (p_a, p_b) = (-1, +1)  ->  P = 0
        eps1 = sigma+, eps2 = sigma-  ->  legs (p_a, p_b) = (+1, +1)  ->  P = +2

    The task brief states the opposite pairing (sigma+sigma- -> 0 and
    sigma+sigma+ -> +2); that is the LADDER reading (rows e of the same table,
    and the [2g] S3.3 probe, whose rows are labelled by the LEG components
    (p_a, p_b), not by the two beams). The reachable SET {0, +-2} is the same
    either way -- only the labelling differs -- and [HAM] S9.5.1(3) is the
    source this package follows.

    THE PER-COMPONENT WEIGHTS BELOW ARE FRAME-SPECIFIC (this test's coplanar
    frame, eps1 = x, eps2 = (cos theta, sin theta, 0)); the frame-INDEPENDENT
    statement, true in any frame for two real linear polarisations at
    relative angle theta, is the rank-resolved sum

        sum_P |w^0_P|^2 = cos^2(theta)/3
        sum_P |w^1_P|^2 = sin^2(theta)/2
        sum_P |w^2_P|^2 = 1/2 + cos^2(theta)/6         (the three sum to 1)

    verified below in BOTH this coplanar frame and a second, z-anchored frame
    (eps1 = z, eps2 = (sin theta, 0, cos theta)), where the per-component
    weights differ (|w^2_{+-1}|^2 = |w^1_{+-1}|^2 = sin^2(theta)/4,
    |w^2_0|^2 = (2/3) cos^2(theta), |w^0_0|^2 = cos^2(theta)/3) but the three
    sums are the same functions of theta.

    The sin^2(theta) Cossel measured ([2g] S3.0; thesis Eq. 6.35, p.219: no
    transfer at theta = 0, maximum at theta = pi/2) is the K = 1 sum:
    sum_P |w^1_P|^2 = sin^2(theta)/2. That is consistent, not contradictory:
    Cossel's 1Sigma+(J=0) -> 3Pi_0+ -> 3Delta1(J=1) transfer has dJ = 1 from
    J = 0, and the 3j (J' K J; ...) admits ONLY K = 1 there, so his
    measurement is a pure-K = 1 observation -- and his single resolved
    intermediate is exactly the restricted manifold in which [HAM] S9.5.3 says
    K = 1 survives at O(1).
    """
    # the anchor the whole convention rests on: c[+1] = +1 for sigma+, and
    # nothing else ([HAM] S9.5.1(3))
    # eps1 = sigma+ gives c_b[+1] = +1 (the anchor, "and nothing else");
    # eps2 = sigma- conjugated gives c_a[+1] = -1; <1 +1 1 +1|2 +2> = 1.
    w = dyad_weights(SIGMA_P, SIGMA_M)
    assert w[(2, 2)] == pytest.approx(-1.0)
    assert abs(dyad_weights(SIGMA_P, SIGMA_P)[(0, 0)]) == pytest.approx(1 / np.sqrt(3))

    reach = {}
    for n1, e1 in (("s+", SIGMA_P), ("s-", SIGMA_M), ("pi", PI_Z)):
        for n2, e2 in (("s+", SIGMA_P), ("s-", SIGMA_M), ("pi", PI_Z)):
            got = sorted({P for (K, P), v in dyad_weights(e1, e2).items()
                          if abs(v) > 1e-12})
            reach[(n1, n2)] = got
    assert reach[("s+", "s+")] == [0], reach          # RAMAN reading, S9.5.1(3)
    assert reach[("s+", "s-")] == [2], reach
    assert reach[("s-", "s+")] == [-2], reach
    assert reach[("s-", "s-")] == [0], reach
    assert reach[("pi", "pi")] == [0], reach
    assert reach[("s+", "pi")] == [1], reach          # negative control: +-1
    assert reach[("pi", "s+")] == [-1], reach         # exists off the sigma pair

    for theta in (0.0, np.pi / 6, np.pi / 3, np.pi / 2):
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([np.cos(theta), np.sin(theta), 0.0])
        w = dyad_weights(e1, e2)
        assert abs(w[(1, 0)]) ** 2 == pytest.approx(np.sin(theta) ** 2 / 2)
        assert w[(1, 1)] == pytest.approx(0.0) and w[(1, -1)] == pytest.approx(0.0)
        assert w[(0, 0)] == pytest.approx(-np.cos(theta) / np.sqrt(3))
        assert w[(2, 0)] == pytest.approx(-np.cos(theta) / np.sqrt(6))
        assert abs(w[(2, 2)]) == pytest.approx(0.5)
        assert abs(w[(2, -2)]) == pytest.approx(0.5)
        # the rank-resolved sum is frame-independent (this coplanar frame)
        s0 = sum(abs(w[(0, P)]) ** 2 for P in (0,))
        s1 = sum(abs(w[(1, P)]) ** 2 for P in (-1, 0, 1))
        s2 = sum(abs(w[(2, P)]) ** 2 for P in (-2, -1, 0, 1, 2))
        assert s0 == pytest.approx(np.cos(theta) ** 2 / 3, abs=1e-12)
        assert s1 == pytest.approx(np.sin(theta) ** 2 / 2, abs=1e-12)
        assert s2 == pytest.approx(0.5 + np.cos(theta) ** 2 / 6, abs=1e-12)

    # the same sums, in a SECOND frame (z-anchored) -- the per-component
    # weights differ from the coplanar frame above, but sum_P |w^K_P|^2 does
    # not: it is a property of theta, not of the frame.
    for theta in (0.0, np.pi / 6, np.pi / 4, np.pi / 2):
        e1 = np.array([0.0, 0.0, 1.0])
        e2 = np.array([np.sin(theta), 0.0, np.cos(theta)])
        w = dyad_weights(e1, e2)
        s0 = sum(abs(w[(0, P)]) ** 2 for P in (0,))
        s1 = sum(abs(w[(1, P)]) ** 2 for P in (-1, 0, 1))
        s2 = sum(abs(w[(2, P)]) ** 2 for P in (-2, -1, 0, 1, 2))
        assert s0 == pytest.approx(np.cos(theta) ** 2 / 3, abs=1e-12)
        assert s1 == pytest.approx(np.sin(theta) ** 2 / 2, abs=1e-12)
        assert s2 == pytest.approx(0.5 + np.cos(theta) ** 2 / 6, abs=1e-12)
    # non-vacuous in the direction that matters: the K = 1 weight really does
    # vanish at theta = 0 and peak at pi/2
    assert abs(dyad_weights(np.array([1.0, 0, 0]), np.array([1.0, 0, 0]))[(1, 0)]) == 0.0


def test_dyad_weights_fails_if_the_two_slots_are_swapped():
    """FAIL demo for the dyad: swapping the two tensor slots (equivalently,
    transposing the Clebsch-Gordan) leaves K = 0 and K = 2 untouched and flips
    the sign of every K = 1 weight, since <1 p_b 1 p_a|K P> = (-1)^K
    <1 p_a 1 p_b|K P>. The completeness identity of [HAM] S9.5.1(3) --
    c_a[p_a] c_b[p_b] = sum_K <1 p_a 1 p_b|K P> w^K_P, i.e. row (a) of its
    printed check -- is what detects it.
    """
    from heff.twophoton import _cg, _leg

    e1, e2 = SIGMA_P, np.array([1.0, 0.0, 0.0])
    good = dyad_weights(e1, e2)
    ca, cb = _leg(np.conj(e2)), _leg(e1)
    bad = {(K, P): (-1.0) ** K * good[(K, P)] for (K, P) in good}   # slots swapped

    def rebuilt(w, pa, pb):
        return sum(_cg(1, pa, 1, pb, K, pa + pb) * w[(K, pa + pb)]
                   for K in (0, 1, 2) if abs(pa + pb) <= K)

    devs = [(abs(rebuilt(good, pa, pb) - ca[pa] * cb[pb]),
             abs(rebuilt(bad, pa, pb) - ca[pa] * cb[pb]))
            for pa in (-1, 0, 1) for pb in (-1, 0, 1)]
    assert max(d[0] for d in devs) < 1e-12, "the honest dyad must reconstruct"
    assert max(d[1] for d in devs) > 0.1, "a swapped dyad should not"


# --------------------------------------------------- registry containment

def test_K2_entries_are_declared_non_hermitian_and_the_builder_accepts_the_registry():
    """Fix round 1, finding 1: the registered fn evaluates the operator at
    P = Delta m_F (heff.terms' fn(bra, ket, ctx) contract), and by reciprocity
    (M_P(a, b) = (-1)^P M_{-P}(b, a), test_rank_K_sum_rule_and_reciprocity)
    that object is antisymmetric under bra<->ket at odd Delta m_F. K = 0 stays
    hermitian=True (Delta m_F = 0 only, so the antisymmetric part is moot);
    K = 2 must be hermitian=False or heff.build_term_matrices raises.
    """
    import heff

    assert REGISTRY_2G["two_photon_K0_dOm0"].hermitian is True
    assert REGISTRY_2G["two_photon_K2_dOm0"].hermitian is False
    assert REGISTRY_2G["two_photon_K2_dOm2"].hermitian is False

    kets, ctx = setup_229(J_max=2)
    d = len(kets)
    devs = {}
    for name, t in REGISTRY_2G.items():
        M = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                if not t.rules.allows(kets[i], kets[j]):
                    continue
                M[i, j] = t.fn(kets[i], kets[j], ctx)
        devs[name] = float(np.max(np.abs(M - M.T)))
    assert devs["two_photon_K0_dOm0"] == 0.0, devs
    assert devs["two_photon_K2_dOm0"] > 0.1, devs
    assert devs["two_photon_K2_dOm2"] > 0.1, devs

    heff.build_term_matrices(kets, ctx, case="c2", registry=REGISTRY_2G)


def test_no_two_photon_term_appears_in_a_hamiltonian_registry():
    """[SPEC-v2] S3.2: REGISTRY_2G is a THIRD registry precisely so a
    transition operator can never be summed into an assembled Hamiltonian.
    Importing heff.twophoton must leave both Hamiltonian registries untouched.
    """
    import heff

    assert set(REGISTRY_2G) & set(REGISTRY) == set()
    assert set(REGISTRY_2G) & set(REGISTRY_C2) == set()
    assert not [n for n in (*REGISTRY, *REGISTRY_C2) if "photon" in n or "alpha" in n]
    assert not [t for t in (*REGISTRY.values(), *REGISTRY_C2.values())
                if any(p.startswith("alpha_K") for p in t.param)]
    # OPEN-21 again, from the parameter side: no alpha_K1_* knob exists
    assert not [k for k in thf_v2("229").params if k.startswith("alpha_K1")]
    for name in ("REGISTRY_2G", "two_photon_geometry", "two_photon_matrix",
                 "dyad_weights", "two_photon_line_strengths"):
        assert name in heff.__all__ and hasattr(heff, name), name
    for t in REGISTRY_2G.values():
        assert "[HAM] S9.5" in t.cite and "5.14" in t.cite, t.name
        assert "bc_5p142_reduced" in t.cite, t.name
