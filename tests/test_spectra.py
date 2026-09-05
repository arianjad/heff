"""Gates B7 (TDM sum rule) and B8 (sum-then-square).

B8 is the one-line ordering choice with a factor-of-anything consequence: the
thesis sums amplitudes and then squares (p.161), so intensity-borrowing paths
interfere. A square-then-sum implementation passes every other test here.
"""
import numpy as np
import pytest

from heff import elements_c  # noqa: F401
from heff.params import thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.spectra import dipole_matrix, line_strengths
from heff.terms import ctx_from


@pytest.fixture
def setup():
    spec = thf_spec()
    kets = enumerate_kets(spec)
    return spec, kets, block_by_mF(kets), ctx_from(spec, thf_v1())


def test_dipole_matrix_connects_only_the_right_mF_blocks(setup):
    spec, kets, blocks, ctx = setup
    a, b = kets[blocks.index[0.5]], kets[blocks.index[1.5]]
    # the 3j (F' 1 F; -m' p m) forces m_bra = m_ket + p, and bra comes from `a`
    assert np.max(np.abs(dipole_matrix(a, b, ctx, 0))) == 0.0        # p = 0 needs m_a = m_b
    assert np.max(np.abs(dipole_matrix(a, b, ctx, -1))) > 0.0        # 0.5 = 1.5 + (-1)
    assert np.max(np.abs(dipole_matrix(a, b, ctx, +1))) == 0.0
    assert np.max(np.abs(dipole_matrix(b, a, ctx, +1))) > 0.0        # and the reverse


def test_B7_total_strength_out_of_a_state_is_independent_of_mF(setup):
    """Gate B7: sum over polarisation and over ALL final basis states of
    |<f|d_p|i>|^2 is the same for every m_F of a given (J, F, Omega).

    Exact, by 3j orthogonality. Uniquely catches a missing or misplaced
    (-1)^(F-m_F) phase and a wrong m-sum -- both of which leave individual
    elements plausible and only break the sum. Analogue of Molecule-Structure
    test_sio_stark_tdm.py:130. F-independence is NOT asserted; only
    m_F-independence, which is the part that follows from orthogonality alone.
    """
    spec, kets, blocks, ctx = setup
    totals = {}
    for J, F, Om in ((1, 1.5, 1.0), (2, 2.5, 1.0), (3, 2.5, -1.0)):
        vals = []
        for mF in np.arange(-F, F + 0.5, 1.0):
            i = int(np.flatnonzero((kets["J"] == J) & (kets["F"] == F)
                                   & (kets["Om"] == Om) & (kets["mF"] == mF))[0])
            tot = 0.0
            for p in (-1, 0, 1):
                row = dipole_matrix(kets[i:i + 1], kets, ctx, p)
                tot += float(np.sum(np.abs(row) ** 2))
            vals.append(tot)
        totals[(J, F, Om)] = vals
        assert np.allclose(vals, vals[0], rtol=1e-10), f"J={J} F={F}: {vals}"
    assert all(v[0] > 0 for v in totals.values())


def test_B7_fails_if_a_polarisation_is_dropped(setup):
    """FAIL demo for gate B7, through the real code path (controller ruling 1):
    the brief's original version built a random matrix and called no heff
    code, which is not a real constraint on this implementation at all.

    Here every dipole_matrix element is exactly the real dipole_geometry
    output -- nothing is corrupted or faked -- but the m-sum over polarisation
    is deliberately incomplete (p = -1 dropped). That is "a wrong m-sum" in
    B7's own docstring: each surviving element is still a perfectly plausible
    number, and only the total, which orthogonality promises is m_F-independent
    when ALL of p in {-1,0,1} are included, comes out m_F-dependent.
    """
    spec, kets, blocks, ctx = setup
    totals = {}
    for J, F, Om in ((1, 1.5, 1.0), (2, 2.5, 1.0), (3, 2.5, -1.0)):
        vals = []
        for mF in np.arange(-F, F + 0.5, 1.0):
            i = int(np.flatnonzero((kets["J"] == J) & (kets["F"] == F)
                                   & (kets["Om"] == Om) & (kets["mF"] == mF))[0])
            tot = 0.0
            for p in (0, 1):                       # p = -1 dropped on purpose
                row = dipole_matrix(kets[i:i + 1], kets, ctx, p)
                tot += float(np.sum(np.abs(row) ** 2))
            vals.append(tot)
        totals[(J, F, Om)] = vals
    assert not all(np.allclose(v, v[0], rtol=1e-10) for v in totals.values()), (
        f"dropping p=-1 should have broken m_F-independence, got {totals}")


def test_B8_amplitudes_are_summed_then_squared():
    """Gate B8: two paths that cancel must give zero strength.

    Hand-built, no physics: the initial eigenvector has equal amplitude on two
    basis kets whose dipole matrix elements to the same final ket are equal and
    OPPOSITE. Sum-then-square gives 0; square-then-sum gives 1/2 * 2 = 1.
    """
    import heff.spectra as S

    kets_a = np.zeros(2, dtype=[("J", "f8"), ("Om", "f8"), ("F", "f8"), ("mF", "f8")])
    kets_b = np.zeros(1, dtype=kets_a.dtype)
    D = np.array([[1.0], [-1.0]])
    va = np.array([[1.0], [1.0]]) / np.sqrt(2)      # (d_a, n_states_a) = (2, 1)
    vb = np.array([[1.0]])
    freqs, strengths = S._strengths_from_matrices([0.0], va, [1.0], vb, {0: D},
                                                  weights={0: 1.0})
    assert strengths.shape == (1, 1)
    assert strengths[0, 0] == pytest.approx(0.0, abs=1e-24)
    # and the same two paths with the SAME sign do not cancel
    D2 = np.array([[1.0], [1.0]])
    _, s2 = S._strengths_from_matrices([0.0], va, [1.0], vb, {0: D2}, weights={0: 1.0})
    assert s2[0, 0] == pytest.approx(2.0)


def test_line_strengths_frequencies_are_upper_minus_lower(setup):
    spec, kets, blocks, ctx = setup
    from heff.assemble import build_term_matrices, hamiltonian

    pset = thf_v1()
    out = []
    for m in (0.5, 1.5):
        sub = kets[blocks.index[m]]
        tm = build_term_matrices(sub, ctx)
        w, v = np.linalg.eigh(hamiltonian(tm, pset, {"E_z": 0.0, "B_z": 0.0}))
        out.append((w, v, sub))
    (wa, va, ka), (wb, vb, kb) = out
    freqs, strengths = line_strengths(wa, va, ka, wb, vb, kb, ctx)
    assert freqs.shape == strengths.shape == (len(wa), len(wb))
    assert freqs[0, 0] == pytest.approx(wb[0] - wa[0])
    assert np.all(strengths >= 0.0)
    assert np.max(strengths) > 0.0
