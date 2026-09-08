"""Gate V27 -- the closure test.

Uniquely catches: a wrong phase or normalisation in the rank-K reduction of
B&C Eq. (5.142), the one step that turns a resolved intermediate sum into a
rank-K reduced element ([HAM] S9.5.1(1)). Every other two-photon gate survives
that error: parity (V24), channel reach (V25), the rank-2 sum rule (V26) and
sum-then-square (V28) never evaluate the closure relation at all -- the 6j
{1 1 K; j' j j''} and the (-1)^(K+j+j') (2K+1)^(1/2) factor of (5.142) appear
NOWHERE in heff.twophoton, which is the whole content of the closure form.
THE INTERMEDIATE MANIFOLD IS SYNTHETIC ON PURPOSE -- the real ThF+ ladder is
contested (docs/digest-literature-two-photon.md S2.2, gap 3: nothing measured
between 3150 and 10472 cm^-1, no 0+/0- labels, no excited-state hyperfine), and
this test checks algebra, not a molecule.

WHAT THE FIXTURE IS. A model whose dipole operator connects the X block
(Omega = +-1) to intermediate electronic states at Omega_i in {0, +-1} and to
nothing else, so the intermediate projector acts as the identity on d|g> and
closure is EXACT by construction. Omega_i = 0 is what opens |dOmega| = 2
(each E1 leg has |q| <= 1, so dOmega = q1 + q2 = +-2 needs q1 = q2 = +-1);
Omega_i = +-1 gives dOmega = 0. One common detuning Delta = 137, which is the
condition [HAM] S9.5.1(2)'s spectator reduction needs (unresolved intermediate
hyperfine). Molecule-frame reduced dipoles 0.7, 1.3, 0.4, -2.1 -- none of them
0 or +-1, so a dropped factor cannot hide. The two legs carry INDEPENDENT
reduced dipoles (D_UP for the bra-side leg, D_DN for the ket-side leg), which
is what makes the test sensitive to the slot assignment of S9.5.1(3).

THE CONTRACTION, DERIVED ONCE ([HAM] S9.5, RAMAN reading -- photon 2 emitted,
hence conjugated; the same reading heff.twophoton.dyad_weights implements).

  T_eff = sum_i (d.eps2*) |i><i| (d.eps1) / Delta
        = (1/Delta) sum_{p_a p_b} c_a[p_a] c_b[p_b] sum_i <f|T^1_{p_a}|i><i|T^1_{p_b}|g>

  with c_a[p] = (-1)^p (eps2*)_{-p} on the BRA-side leg and c_b[p] =
  (-1)^p (eps1)_{-p} on the KET-side leg (S9.5.1(3)), written out here from
  that line rather than imported, so V27 cross-checks `dyad_weights` too.

  RESOLVED SIDE. Each leg is elements_c2.axial_geometry at k = 1 with the
  molecule-frame component q forced by the Omega's, times the reduced dipole
  of that (Omega, Omega_i) pair. The intermediate sum runs over
  (Omega_i, J'', F1'', F'', m''_F) with the FULL hyperfine set at each
  (Omega_i, J'') and J'' = 0..3, complete for an X block at J = 1, 2.

  RANK-K SIDE.  A = sum_K sum_P w^K_P alpha^K(Om_f, Om_g) G_K(f, g, P),
  w^K_P = dyad_weights(eps1, eps2), G_K = two_photon_geometry (K = 0, 2, the
  registered channels) and axial_geometry at k = 1 for the unregistered K = 1.
  The geometry is NOT rescaled: heff's normalisation is
  <eta'||alpha^K||eta> == 1 per channel (conventions.two_photon_norm =
  'bc_5p142_reduced'), so the closure-valued alpha is computed here, on the
  resolved side, from the same reduced dipoles -- B&C (5.141) inverted in the
  MOLECULE frame, which is (5.142) read backwards one level in:

      alpha^K_q = (1/Delta) (-1)^q (2K+1)^(1/2)
                  sum_{q1 + q2 = q} (1 1 K; q1 q2 -q) D_UP[Om_f, Om_i] D_DN[Om_i, Om_g]

  with Om_i = Om_g + q2. No spectator 6j enters this independent coefficient.

K = 1, [HAM] S9.5.3 / ruling R10. This fixture's manifold is RESTRICTED in
Omega (no Omega_i = +-2), which is the physical ThF+/HfF+ situation, so K = 1
survives at O(1): the K = 1 piece is 6-9 % of the amplitude here (and exactly
zero for sigma+sigma-, where the K = 1 dyad weights vanish). It is reported,
never asserted against. What kills it is COMPLETENESS OF THE MANIFOLD, not a
common Delta and not the J'' range: test_K1_vanishes_when_the_manifold_is_complete
adds Omega_i = +-2 with a dipole product matched to the Omega_i = 0 path and
alpha^1 drops to 0.0 exactly, while the J'' range governs the closure identity
for EVERY K (test_closure_fails_with_a_truncated_intermediate_J).
"""
import numpy as np
import pytest
import test_twophoton as V

import heff.twophoton as tp
from heff.elements_c2 import axial_geometry
from heff.params import thf_v2
from heff.spec import KET_C2, enumerate_kets, thf_spec
from heff.terms import ctx_from
from heff.twophoton import dyad_weights, two_photon_matrix
from heff.wigner import w3j, w6j

SQ2 = np.sqrt(2.0)
DELTA = 137.0
#: <eta_X, Om_f| d_{q} |eta_i, Om_i>, the BRA-side (slot 1) leg. q = Om_f - Om_i.
D_UP = {(1.0, 0.0): 0.7, (-1.0, 0.0): 0.7, (1.0, 1.0): 1.3, (-1.0, -1.0): 1.3}
#: <eta_i, Om_i| d_{q} |eta_X, Om_g>, the KET-side (slot 2) leg. q = Om_i - Om_g.
D_DN = {(0.0, 1.0): 0.4, (0.0, -1.0): 0.4, (1.0, 1.0): -2.1, (-1.0, -1.0): -2.1}
OM_I = (0.0, 1.0, -1.0)
J_INT = (0, 3)                      # complete for an X block at J = 1, 2

#: (label, eps1, eps2) -- sigma+sigma-, two linear at 45 deg, and two controls.
POLS = (("sigma+ sigma-", V.SIGMA_P, V.SIGMA_M),
        ("linear 45 deg", np.array([1.0, 0.0, 0.0]),
         np.array([1.0, 1.0, 0.0]) / SQ2),
        ("sigma+ pi", V.SIGMA_P, V.PI_Z),
        ("sigma+ sigma+", V.SIGMA_P, V.SIGMA_P))


def _leg(eps):
    """c[p] = (-1)^p eps_{-p}, [HAM] S9.5.1(3), written out from that line."""
    ex, ey, ez = (complex(c) for c in np.asarray(eps).reshape(3))
    s = {1: -(ex + 1j * ey) / SQ2, 0: ez, -1: (ex - 1j * ey) / SQ2}
    return {p: (-1.0) ** p * s[-p] for p in (-1, 0, 1)}


def _inter_kets(Om_i, ctx, J_range):
    """|((J I_Th) F1, I_F) F, m_F> at one Omega_i -- spec._enumerate_kets_two_spin
    with the (-Om0, +Om0) loop replaced by the one signed Omega the caller asks
    for, since a StateSpec cannot carry Omega = 0 as a two-sided pair."""
    t1, t2 = int(round(2 * ctx.spins[0].I)), int(round(2 * ctx.I))
    rows = [(float(J), float(Om_i), twoF1 / 2.0, twoF / 2.0, twom / 2.0)
            for J in range(J_range[0], J_range[1] + 1) if abs(Om_i) <= J
            for twoF1 in range(abs(2 * J - t1), 2 * J + t1 + 1, 2)
            for twoF in range(abs(twoF1 - t2), twoF1 + t2 + 1, 2)
            for twom in range(-twoF, twoF + 1, 2)]
    return np.array(rows, dtype=KET_C2)


def _leg_matrix(bras, kets, ctx, p, dip):
    """One E1 leg: axial_geometry(k=1) x the reduced dipole of its Omega pair.

    The mask is k = 1 selection rules (dm_F = p, |dJ| <= 1, |dF| <= 1,
    |dF1| <= 1) plus "this Omega pair has a dipole"; anything it wrongly
    dropped would show up as a closure failure, so it is self-validating.
    """
    d = np.zeros((len(bras), len(kets)))
    for (Ob, Ok), val in dip.items():
        d[np.ix_(bras["Om"] == Ob, kets["Om"] == Ok)] = val
    ok = ((np.abs(bras["mF"][:, None] - kets["mF"][None, :] - p) < 1e-9)
          & (np.abs(bras["J"][:, None] - kets["J"][None, :]) <= 1.0)
          & (np.abs(bras["F"][:, None] - kets["F"][None, :]) <= 1.0)
          & (np.abs(bras["F1"][:, None] - kets["F1"][None, :]) <= 1.0)
          & (d != 0.0))
    out = np.zeros((len(bras), len(kets)))
    for i, j in zip(*np.nonzero(ok)):
        out[i, j] = d[i, j] * axial_geometry(
            bras[i], kets[j], ctx, k=1,
            q=float(bras["Om"][i]) - float(kets["Om"][j]), p=p)
    return out


def closure_alpha(K, Om_f, Om_g, d_up=D_UP, d_dn=D_DN):
    """alpha^K_q = (1/Delta) <eta_f|T^K_q(d, d)|eta_g>, q = Om_f - Om_g.

    B&C (5.141) at k1 = k2 = 1 in the MOLECULE frame -- the closure form of
    [HAM] S9.5.1(1) with the intermediate sum done explicitly over the
    fixture's Omega_i. Free of 6j's by construction.
    """
    q = Om_f - Om_g
    if abs(q) > K:
        return 0.0
    tot = 0.0
    for q1 in (-1.0, 0.0, 1.0):
        q2 = q - q1
        if abs(q2) > 1.0:
            continue
        up, dn = d_up.get((Om_f, Om_g + q2)), d_dn.get((Om_g + q2, Om_g))
        if up is not None and dn is not None:
            tot += w3j(1, 1, K, q1, q2, -q) * up * dn
    return (-1.0) ** q * np.sqrt(2 * K + 1) * tot / DELTA


def build_legs(X, ctx, J_range=J_INT, om_i=OM_I, d_up=D_UP, d_dn=D_DN):
    """{(Omega_i, p): matrix} for both legs over the whole intermediate manifold."""
    ups, dns = {}, {}
    for Om_i in om_i:
        I = _inter_kets(Om_i, ctx, J_range)
        if len(I) == 0:
            continue
        for p in (-1, 0, 1):
            ups[(Om_i, p)] = _leg_matrix(X, I, ctx, p, d_up)
            dns[(Om_i, p)] = _leg_matrix(I, X, ctx, p, d_dn)
    return ups, dns


def resolved(X, ups, dns, eps1, eps2):
    """The resolved amplitude matrix, sum_i (d.eps2*)|i><i|(d.eps1)/Delta."""
    ca, cb = _leg(np.conj(np.asarray(eps2, dtype=complex))), _leg(eps1)
    A = np.zeros((len(X), len(X)), dtype=complex)
    for pa in (-1, 0, 1):
        for pb in (-1, 0, 1):
            if abs(ca[pa]) == 0.0 or abs(cb[pb]) == 0.0:
                continue
            A += ca[pa] * cb[pb] * sum(
                ups[key] @ dns[(key[0], pb)] for key in ups if key[1] == pa)
    return A / DELTA


def rank_geometry(X, ctx):
    """{(K, P): matrix}. K = 0, 2 through the SHIPPED two_photon_geometry (the
    registered channels), K = 1 through axial_geometry at k = 1 -- there is no
    K = 1 operator to call ([HAM] S9.5.3, OPEN-21)."""
    G = {}
    for K in (0, 1, 2):
        for P in range(-K, K + 1):
            M = np.zeros((len(X), len(X)))
            for i, j in zip(*np.nonzero(
                    np.abs(X["mF"][:, None] - X["mF"][None, :] - P) < 1e-9)):
                q = float(X["Om"][i]) - float(X["Om"][j])
                if abs(q) > K:
                    continue
                M[i, j] = (tp.two_photon_geometry(X[i], X[j], ctx, K=K, P=P)
                           if K in (0, 2) else
                           axial_geometry(X[i], X[j], ctx, k=1, q=q, p=P))
            G[(K, P)] = M
    return G


def rank_amplitude(X, G, eps1, eps2, Ks=(0, 1, 2), alpha=closure_alpha):
    """sum_{K in Ks} sum_P w^K_P alpha^K(Om_f, Om_g) G_K(f, g, P)."""
    w = dyad_weights(eps1, eps2)
    A = np.zeros((len(X), len(X)), dtype=complex)
    for K in Ks:
        al = np.zeros((len(X), len(X)))
        for a in (1.0, -1.0):
            for b in (1.0, -1.0):
                al[np.ix_(X["Om"] == a, X["Om"] == b)] = alpha(K, a, b)
        for P in range(-K, K + 1):
            A += w[(K, P)] * al * G[(K, P)]
    return A


@pytest.fixture(scope="module")
def fix():
    """(X kets, ctx, leg matrices, rank-K geometry, resolved amplitudes).

    Cache the independent resolved and rank-K calculations for the module.
    """
    spec = thf_spec("229", J_max=2)
    X, ctx = enumerate_kets(spec), ctx_from(spec, thf_v2("229"))
    ups, dns = build_legs(X, ctx)
    return {"X": X, "ctx": ctx, "ups": ups, "dns": dns,
            "G": rank_geometry(X, ctx),
            "A": {n: resolved(X, ups, dns, e1, e2) for n, e1, e2 in POLS}}


# ------------------------------------------------------------------ V27

def test_V27_resolved_intermediate_sum_equals_the_rank_K_form(fix):
    """V27: the resolved sum over the synthetic intermediate manifold equals the
    rank-K closure form, to atol = 1e-12, for four polarisation pairs.

    This is B&C (5.142) as an IDENTITY between two computations that share no
    algebra: the left side runs elements_c2.axial_geometry twice at k = 1 and
    sums over every intermediate (Omega_i, J'', F1'', F'', m''_F); the right
    side runs it once at k = K and multiplies by the closure alpha. Their
    agreement is the {1 1 K; j' j j''} 6j and the (-1)^(K+j+j') (2K+1)^(1/2)
    phase of (5.142), at all three spectator levels at once.
    """
    for name, e1, e2 in POLS:
        A = fix["A"][name]
        dev = np.max(np.abs(A - rank_amplitude(fix["X"], fix["G"], e1, e2)))
        assert np.max(np.abs(A)) > 1e-3, f"{name}: vacuous, the amplitude is dead"
        assert dev < 1e-12, f"{name}: max|resolved - rank-K| = {dev:.3e}"


def test_V27_registered_channels_alone_carry_the_sigma_pair(fix):
    """The K = 0 and K = 2 projections against the REGISTERED channels only.

    For eps1 = sigma+, eps2 = sigma- the K = 1 dyad weights vanish identically
    (eps2* is parallel to eps1, so eps2* x eps1 = 0), so the two registered
    channels alone must reproduce the whole resolved amplitude. For the other
    three pairs they cannot, and the shortfall IS the K = 1 piece -- printed,
    not asserted against ([HAM] S9.5.3, ruling R10: this fixture's manifold is
    restricted in Omega, so K = 1 survives at O(1)).
    """
    w1 = dyad_weights(V.SIGMA_P, V.SIGMA_M)
    assert max(abs(w1[(1, P)]) for P in (-1, 0, 1)) == 0.0
    report = {}
    for name, e1, e2 in POLS:
        A = fix["A"][name]
        K02 = rank_amplitude(fix["X"], fix["G"], e1, e2, Ks=(0, 2))
        K1 = rank_amplitude(fix["X"], fix["G"], e1, e2, Ks=(1,))
        report[name] = (np.max(np.abs(A - K02)), np.max(np.abs(K1)),
                        np.max(np.abs(A)))
        assert np.max(np.abs(A - K02 - K1)) < 1e-12, name
    print("\n  K=1 residual (restricted Omega manifold, [HAM] S9.5.3):")
    for name, (dev, k1, amp) in report.items():
        print(f"    {name:14s} |A - (K0+K2)| = {dev:.3e}  |K1| = {k1:.3e}"
              f"  ({100 * k1 / amp:.1f} % of max|A| = {amp:.3e})")
    assert report["sigma+ sigma-"][0] < 1e-12
    assert all(r[1] > 1e-5 for n, r in report.items() if n != "sigma+ sigma-")


def test_closure_alpha_is_one_scalar_per_channel(fix):
    """[HAM] S9.5.2's channel structure, read off the closure alphas.

    alpha^K is the SAME number for Omega = +1 and Omega = -1 at K = 0 and
    K = 2, and the same for q = +2 and q = -2 -- which CORROBORATES (it does not
    derive) heff.twophoton carrying one alpha per (K, |dOmega|) channel, and
    makes
    [T, parity] = 0 (gate V24) true rather than assumed. At K = 1 the two
    Omega-diagonal values are OPPOSITE, an independent corroboration that a
    K = 1 term has no place in the parity-even registered set.
    """
    for K in (0, 2):
        assert closure_alpha(K, 1.0, 1.0) == pytest.approx(
            closure_alpha(K, -1.0, -1.0))
        assert closure_alpha(K, 1.0, -1.0) == pytest.approx(
            closure_alpha(K, -1.0, 1.0))
    assert closure_alpha(1, 1.0, 1.0) == pytest.approx(
        -closure_alpha(1, -1.0, -1.0))
    assert abs(closure_alpha(1, 1.0, 1.0)) > 1e-4              # non-vacuous
    # |dOmega| = 2 exists at K = 2 only ([HAM] S9.5.2)
    assert closure_alpha(0, 1.0, -1.0) == 0.0
    assert abs(closure_alpha(2, 1.0, -1.0)) > 1e-4
    print("\n  closure alphas (Delta = %.1f):" % DELTA)
    for K in (0, 1, 2):
        print(f"    K={K}  " + "  ".join(
            f"({a:+.0f}<-{b:+.0f}) {closure_alpha(K, a, b):+.8f}"
            for a, b in ((1.0, 1.0), (-1.0, -1.0), (1.0, -1.0))))


def test_closure_fails_with_a_truncated_intermediate_J(fix):
    """Non-vacuity of V27, and the validity condition it rests on.

    The closure identity needs the intermediate projector to act as the
    identity on d|g>, so the J'' range must cover |J -+ 1| for every J in the X
    block. Truncating it to J'' = 1, 2 (the X block runs to J = 2, so J'' = 3
    is missing) breaks the identity for EVERY K by 1e-3 .. 8e-3 on an amplitude
    of order 6e-3 -- i.e. V27 can fail, and this is what it is measuring.
    """
    X, ctx = fix["X"], fix["ctx"]
    ups, dns = build_legs(X, ctx, J_range=(1, 2))
    devs = {n: np.max(np.abs(resolved(X, ups, dns, e1, e2)
                             - rank_amplitude(X, fix["G"], e1, e2)))
            for n, e1, e2 in POLS}
    print("\n  truncated J'' = 1..2, closure deviation:")
    for n, d in devs.items():
        print(f"    {n:14s} {d:.3e}")
    assert min(devs.values()) > 1e-3, devs


def test_K1_vanishes_when_the_manifold_is_complete(fix):
    """[HAM] S9.5.3, ruling R10, with both outcomes reachable.

    S9.5.3 is explicit that two suppressions of K = 1 must not be conflated:
    COMPLETENESS of the intermediate manifold (exact zero) and COMMONNESS of
    the detuning (which alone buys nothing once the manifold is restricted).
    This fixture already has a common Delta and, in the test above, a complete
    J'' range -- and K = 1 is still O(1), because the manifold is restricted in
    Omega. Completing it in Omega is what kills K = 1: add Omega_i = +-2 with
    D_UP[+1,+2] D_DN[+2,+1] = D_UP[+1,0] D_DN[0,+1], which makes
    P_X d_a P d_b P_X symmetric in a <-> b, and alpha^1 goes to exactly 0 while
    the closure identity keeps holding at 1e-18.
    """
    d_up = {**D_UP, (1.0, 2.0): 0.35, (-1.0, -2.0): 0.35}
    d_dn = {**D_DN, (2.0, 1.0): 0.8, (-2.0, -1.0): 0.8}     # 0.35*0.8 == 0.7*0.4
    X, ctx = fix["X"], fix["ctx"]
    ups, dns = build_legs(X, ctx, om_i=(0.0, 1.0, -1.0, 2.0, -2.0),
                          d_up=d_up, d_dn=d_dn)

    def alpha(K, a, b):
        return closure_alpha(K, a, b, d_up=d_up, d_dn=d_dn)

    assert alpha(1, 1.0, 1.0) == 0.0 and alpha(1, -1.0, -1.0) == 0.0
    assert abs(alpha(0, 1.0, 1.0)) > 1e-3 and abs(alpha(2, 1.0, 1.0)) > 1e-3
    for name, e1, e2 in POLS:
        A = resolved(X, ups, dns, e1, e2)
        dev = np.max(np.abs(A - rank_amplitude(X, fix["G"], e1, e2, Ks=(0, 2),
                                               alpha=alpha)))
        assert np.max(np.abs(A)) > 1e-3, name
        assert dev < 1e-12, f"{name}: {dev:.3e}"
