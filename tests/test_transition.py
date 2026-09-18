import numpy as np
import pytest

from heff.assemble import hamiltonian
from heff.spec import block_by_mF
from heff.spectra import dipole_matrix  # noqa: F401  (import check)
from heff.transition import (DipoleOperator, TwoPhotonOperator, diagonalize, padded_thf,
                             select, sweep_transition_matrix, transition_matrix)
from heff.twophoton import POLARIZATIONS, two_photon_line_strengths
from heff.wigner import w3j, w6j

SP, SM, PI, X, Y = (POLARIZATIONS[k] for k in ("sigma+", "sigma-", "pi", "x", "y"))
# Raman reading (default): photon 1 absorbed, photon 2 emitted, Delta m_F = p1 - p2.
PAIRS = {("sigma+", "sigma-"): 2, ("sigma-", "sigma+"): -2, ("pi", "pi"): 0,
         ("sigma+", "sigma+"): 0, ("sigma+", "pi"): 1, ("pi", "sigma+"): -1}


@pytest.fixture(scope="module")
def sys232():
    kets, ctx, tm, pset = padded_thf("232", J_max=3)
    op = TwoPhotonOperator()
    return kets, ctx, tm, pset, op, op.channels(kets, kets, ctx)


def _eig(s, E_z=0.0, B_z=1.0, **kw):
    kets, ctx, tm, pset, op, ch = s
    return diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=B_z, **kw)


def test_gate1_delta_mF_equals_p1_minus_p2(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2, 3))
    for (e1, e2), dm in PAIRS.items():
        t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1=e1, eps2=e2)
        i, j = np.nonzero(t.strength > 1e-12)
        assert len(i) > 0
        got = {t.labels_b[b]["mF"] - t.labels_a[a]["mF"] for a, b in zip(i, j)}
        assert got == {float(dm)}, (e1, e2, got)


def test_gate2_parity_conserved_at_zero_E(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=1.0)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2, 3))
    t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="pi", eps2="pi")
    i, j = np.nonzero(t.strength > 1e-12)
    assert all(t.labels_a[a]["parity"] == t.labels_b[b]["parity"] for a, b in zip(i, j))
    assert all(t.labels_a[a]["ef"] == t.labels_b[b]["ef"] for a, b in zip(i, j)
               if t.labels_a[a]["J"] == t.labels_b[b]["J"])
    eigE = _eig(sys232, E_z=20.0, B_z=1.0)
    tE = transition_matrix(op, eigE, select(eigE, J=1), select(eigE, J=(1, 2, 3)), ctx,
                           channels=ch, eps1="pi", eps2="pi")
    i, j = np.nonzero(tE.strength > 1e-6)
    assert any(tE.labels_a[a]["parity"] != tE.labels_b[b]["parity"] for a, b in zip(i, j))


def test_gate3_sum_rule_in_eigenbasis(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=1e-3)
    rows, cols = select(eig, J=1), np.arange(len(kets))
    for dOm in (0, 2):
        tot = np.zeros(len(rows))
        for P in range(-2, 3):
            M = ch[(2, dOm, P)]
            A = eig.evecs[:, rows].conj().T @ M @ eig.evecs[:, cols]
            tot += np.sum(np.abs(A) ** 2, axis=1)
        assert np.allclose(tot, 1.0, rtol=1e-6), (dOm, tot)


def test_gate4_random_phases_leave_strength_unchanged(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=2)
    t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="x", eps2="y")
    rng = np.random.default_rng(0)
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi, len(kets)))
    from heff.transition import Eigensystem
    eig2 = Eigensystem(eig.kets, eig.evals, eig.evecs * ph[None, :], eig.labels, eig.field)
    t2 = transition_matrix(op, eig2, rows, cols, ctx, channels=ch, eps1="x", eps2="y")
    assert np.allclose(t.strength, t2.strength)
    assert not np.allclose(t.amp, t2.amp)


def test_gate5_blockwise_equals_full(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    a, b = _eig(sys232, B_z=2.0), _eig(sys232, B_z=2.0, blockwise=False)
    ra, ca = select(a, J=1), select(a, J=(1, 2))
    rb, cb = select(b, J=1), select(b, J=(1, 2))
    ta = transition_matrix(op, a, ra, ca, ctx, channels=ch, eps1="sigma+", eps2="pi")
    tb = transition_matrix(op, b, rb, cb, ctx, channels=ch, eps1="sigma+", eps2="pi")
    assert np.allclose(a.evals, b.evals, atol=1e-8)
    assert np.allclose(ta.strength, tb.strength, atol=1e-10)


def test_gate6_x_is_the_phased_sigma_sum(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2))
    tx = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="x", eps2="x")
    parts = {}
    for e1 in ("sigma+", "sigma-"):
        for e2 in ("sigma+", "sigma-"):
            parts[(e1, e2)] = transition_matrix(op, eig, rows, cols, ctx, channels=ch,
                                                eps1=e1, eps2=e2).amp
    # x = (sigma- - sigma+)/sqrt2 in the Condon-Shortley vectors above
    coh = 0.5 * (parts[("sigma-", "sigma-")] - parts[("sigma-", "sigma+")]
                 - parts[("sigma+", "sigma-")] + parts[("sigma+", "sigma+")])
    assert np.allclose(tx.amp, coh, atol=1e-12)
    incoh = 0.25 * sum(np.abs(v) ** 2 for v in parts.values())
    assert not np.allclose(tx.strength, incoh)


def test_gate7_matches_two_photon_line_strengths(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=0.0)
    ma, mb = select(eig, mF=0.5), select(eig, mF=2.5)
    # mF 0.5 -> 2.5 is Delta m_F = +2: (sigma+ absorbed, sigma- emitted) in the Raman reading
    t = transition_matrix(op, eig, ma, mb, ctx, channels=ch, eps1="sigma+", eps2="sigma-")
    blk = block_by_mF(kets)
    ia, ib = blk.index[0.5], blk.index[2.5]
    H = hamiltonian(tm, pset, {"E_z": 0.0, "B_z": 0.0})
    wa, va = np.linalg.eigh(H[np.ix_(ia, ia)])
    wb, vb = np.linalg.eigh(H[np.ix_(ib, ib)])
    # kets_a is the bra (final) side of two_photon_matrix: pass the mF=2.5 block first
    _, S = two_photon_line_strengths(wb, vb, kets[ib], wa, va, kets[ia], ctx,
                                     eps1=SP, eps2=SM, alphas=op.alphas)
    assert np.allclose(np.sort(t.strength.ravel()), np.sort(S.ravel()), atol=1e-10)


def test_gate8_padding_matches_one_spin(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    I = 0.5
    def geo1(bra, ket, k, q, p):
        J1, F1, m1, O1 = (float(ket[x]) for x in ("J", "F", "mF", "Om"))
        J2, F2, m2, O2 = (float(bra[x]) for x in ("J", "F", "mF", "Om"))
        if abs(m2 - m1 - p) > 1e-9 or abs(O2 - O1 - q) > 1e-9:
            return 0.0
        return ((-1) ** (F2 - m2) * w3j(F2, k, F1, -m2, p, m1)
                * (-1) ** (F1 + J2 + k + I) * np.sqrt((2 * F2 + 1) * (2 * F1 + 1)) * w6j(J1, F1, I, F2, J2, k)
                * (-1) ** (J2 - O2) * np.sqrt((2 * J2 + 1) * (2 * J1 + 1)) * w3j(J2, k, J1, -O2, q, O1))
    for (K, dOm, P), A in ch.items():
        B = np.zeros_like(A)
        for i in range(len(kets)):
            for j in range(len(kets)):
                q = float(kets["Om"][i]) - float(kets["Om"][j])
                if abs(q) == dOm:
                    B[i, j] = geo1(kets[i], kets[j], K, q, P)
        assert np.max(np.abs(A - B)) < 1e-14, (K, dOm, P)
    from heff.assemble import build_term_matrices
    from heff.params import thf_v1
    from heff.spec import enumerate_kets, thf_spec
    from heff.terms import REGISTRY, ctx_from
    v1 = thf_spec("232", J_max=3); k1 = enumerate_kets(v1)
    tm1 = build_term_matrices(k1, ctx_from(v1, thf_v1()), case="c", registry=REGISTRY)
    for E, B in ((0.0, 0.0), (10.0, 5.0)):
        H1 = hamiltonian(tm1, thf_v1(), {"E_z": E, "B_z": B})
        H2 = hamiltonian(tm, pset, {"E_z": E, "B_z": B})
        assert np.max(np.abs(H1 - H2)) < 1e-12


def test_sweep_is_continuous(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    B = np.linspace(1e-3, 20.0, 41)
    amps, freqs, la, lb = sweep_transition_matrix(
        op, kets, tm, pset, ctx, E_z=0.0, B_values=B,
        rows_at=lambda e: select(e, J=1), cols_at=lambda e: select(e, J=2),
        eps1="sigma+", eps2="sigma+")
    S = np.abs(amps) ** 2
    assert np.max(np.abs(np.diff(S, axis=0))) < 0.05 * max(S.max(), 1e-12)


def test_dipole_operator_selection_rule(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    d = DipoleOperator()
    t = transition_matrix(d, eig, select(eig, J=1), select(eig, J=2), ctx, eps="sigma+")
    i, j = np.nonzero(t.strength > 1e-12)
    assert {t.labels_b[b]["mF"] - t.labels_a[a]["mF"] for a, b in zip(i, j)} == {1.0}
