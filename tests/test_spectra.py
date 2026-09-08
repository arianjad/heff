"""Gates B7/B7b (TDM sum rule / spherical-tensor reciprocity) and B8
(sum-then-square).

B8 is the one-line ordering choice with a factor-of-anything consequence: the
thesis sums amplitudes and then squares (p.161), so intensity-borrowing paths
interfere. A square-then-sum implementation passes every other test here.

B7 is a closed-form MAGNITUDE identity (sums |amplitude|^2), so it cannot see a
phase error; B7b is the companion identity that lives entirely in phases and
signs.
"""
import numpy as np
import pytest

from heff.params import thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.spectra import dipole_matrix, label_lines, line_strengths
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
    |<f|d_p|i>|^2 equals exactly 1.0 (dimensionless geometry units -- the
    closed form <n_hat.n_hat> = 1, in units of d_mf^2) for every m_F of a
    given (J, F, Omega), for J <= 3 in this J = 1-4 basis.

    Exact, by 3j/6j completeness -- but only while the basis actually holds
    every J' a dipole transition can reach (J' = J, J +- 1). J = 4 would need
    J' = 5, which this basis truncates, so J = 4 gives 7/15 instead of 1 and
    is deliberately excluded from the loop below; it is not asserted anywhere
    in this file.

    This gate sums |amplitude|^2, so no unimodular phase -- including a missing
    or misplaced (-1)^(F-m_F) -- changes the sum. It detects a wrong m-sum or
    dropped polarisation channel; B7b checks phase by spherical-tensor
    reciprocity. Analogue of Molecule-Structure test_sio_stark_tdm.py:130.
    The test asserts m_F-independence, not F-independence.
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
        assert vals[0] == pytest.approx(1.0, rel=1e-10), (
            f"J={J} F={F}: sum rule = {vals[0]}, expected the closed form 1.0")
    assert all(v[0] > 0 for v in totals.values())


def test_B7b_transition_dipole_obeys_spherical_tensor_reciprocity(setup):
    """Gate B7b checks the phase relation that B7 cannot detect.

    dipole_geometry never returns a phase (report's self-review: always
    float64), so the general Wigner-Eckart Hermiticity relation
    T^1_p+ = (-1)^p T^1_{-p} collapses, for THIS real operator, to a bare sign:
    <a|d_p|b> = (-1)^p <b|d_{-p}|a>, i.e. D(a, b, p) == (-1)^p * D(b, a, -p).T.
    Checked for p = +-1 across the two mF blocks that are actually connected
    (both directions), plus the p = 0 case, which only ever connects a block
    to itself (m_bra = m_ket + 0), so "the two blocks" collapse to one and the
    identity reduces to plain symmetry, D(a, a, 0) == D(a, a, 0).T.
    """
    spec, kets, blocks, ctx = setup
    a, b = kets[blocks.index[0.5]], kets[blocks.index[1.5]]

    for p in (-1, 1):
        lhs = dipole_matrix(a, b, ctx, p)
        rhs = (-1.0) ** p * dipole_matrix(b, a, ctx, -p).T
        assert np.allclose(lhs, rhs), f"p={p}: reciprocity broken"
    # non-vacuous: 0.5 = 1.5 + (-1) is the allowed direction, so this is > 0
    assert np.max(np.abs(dipole_matrix(a, b, ctx, -1))) > 0.0

    D0 = dipole_matrix(a, a, ctx, 0)
    assert np.max(np.abs(D0)) > 0.0
    assert np.allclose(D0, D0.T)


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


def test_line_strengths_accepts_complex_polarisation_weights(setup):
    """A complex weight (for example, a sigma+/sigma- phase)
    must not raise, and must give the same strength as an incoherent sum for
    a pair of DEFINITE-m_F blocks -- there is no interference to get wrong
    there (finding 4: the 3j selection rule already forces exactly one p
    non-zero for any state pair across two definite-m_F blocks), so this is
    an exact check, not an approximation.
    """
    spec, kets, blocks, ctx = setup
    a, b = kets[blocks.index[0.5]], kets[blocks.index[1.5]]
    eva, evb = np.eye(len(a)), np.eye(len(b))  # identity evecs: strengths = |D_p|^2 directly
    weights = {-1: 1 / np.sqrt(2), +1: 1j / np.sqrt(2)}

    Dm1 = dipole_matrix(a, b, ctx, -1)
    Dp1 = dipole_matrix(a, b, ctx, +1)
    assert np.max(np.abs(Dp1)) == 0.0  # only p = -1 connects 0.5 -> 1.5

    _, strengths = line_strengths(np.zeros(len(a)), eva, a, np.zeros(len(b)), evb, b,
                                  ctx, weights=weights)
    assert np.isrealobj(strengths)
    expected = abs(weights[-1]) ** 2 * Dm1 ** 2 + abs(weights[+1]) ** 2 * Dp1 ** 2
    assert np.allclose(strengths, expected)
    assert np.max(strengths) > 0.0


def test_line_labels_are_present_and_the_strongest_J1_to_J2_lines_flip_parity(setup):
    """`label_lines` records a superposition parity and
    an e/f name for every state, and the six strongest J = 1 -> J = 2 lines at
    zero field -- a Delta-m_F = 0 spectrum within one signed-m_F block -- all
    connect states of OPPOSITE parity (a real selection-rule consequence of
    the E1 operator being parity-odd, not a filter: every J = 1 -> J = 2 pair
    is included in `strengths`, labels only annotate it after the fact).
    """
    spec, kets, blocks, ctx = setup
    from heff.assemble import build_term_matrices, hamiltonian

    sub = kets[blocks.index[0.5]]
    tm = build_term_matrices(sub, ctx)
    w, v = np.linalg.eigh(hamiltonian(tm, thf_v1(), {"E_z": 0.0, "B_z": 0.0}))
    freqs, strengths = line_strengths(w, v, sub, w, v, sub, ctx)

    labels = label_lines(sub, v, ctx.S, rule=ctx.conventions.ef_rule, ell=0.0, s=0.0)
    assert len(labels) == len(sub)
    assert all(l["ef"] in ("e", "f") and l["parity"] in (-1, 1) for l in labels)

    J = np.array([l["J"] for l in labels])
    mask = np.outer(J == 1.0, J == 2.0)
    masked = np.where(mask, strengths, -1.0)  # exclude everything but J=1 -> J=2
    top6 = np.argsort(masked.ravel())[::-1][:6]
    for i, j in zip(*np.unravel_index(top6, masked.shape)):
        assert masked[i, j] > 0.0
        assert labels[i]["parity"] == -labels[j]["parity"], (
            f"line {i}->{j} (freq {freqs[i, j]:.3f} MHz) did not flip parity")
