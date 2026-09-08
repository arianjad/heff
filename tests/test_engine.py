"""The 1D/2D sweep entry point that C2V-Molecules does not have: one that
returns EIGENVECTORS. C2V's batch_diagonalize_scan returns (evals, expects) and
discards vectors (C2V digest Q3, explicit); v1 needs vectors for TDMs, parity
labels and expectation values, so this is not optional.
"""
import numpy as np
import pytest

from heff.assemble import build_term_matrices, hamiltonian
from heff.engine import eigh_batch, sweep
from heff.params import thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import ctx_from


@pytest.fixture
def tm():
    spec = thf_spec()
    kets = enumerate_kets(spec)
    idx = block_by_mF(kets).index[1.5]
    return build_term_matrices(kets[idx], ctx_from(spec, thf_v1()))


def test_eigh_batch_matches_per_point_eigh_and_respects_chunking():
    """Chunked eigh_batch (chunk=4 forces >= 3 chunks on n=23) agrees with the
    unbatched per-point np.linalg.eigh reference and stays energy-ascending."""
    rng = np.random.default_rng(42)
    H = rng.standard_normal((23, 5, 5))
    H = H + np.transpose(H, (0, 2, 1))
    w, v = eigh_batch(H, chunk=4)
    for n in range(23):
        wn, vn = np.linalg.eigh(H[n])
        assert np.allclose(w[n], wn)
        assert np.allclose(np.abs(v[n]), np.abs(vn))
    assert np.allclose(w, np.sort(w, axis=1))


def test_eigh_batch_default_chunk_bytes_forces_chunking_and_agrees_with_unchunked():
    """Ruling 1: chunk sizing is driven by a chunk_bytes memory budget, not a
    hardcoded point count. A budget too small for even one point's worth of
    slack still yields one point per chunk (n=10, so >= 3 chunks) and must
    still agree with a single unchunked call."""
    rng = np.random.default_rng(99)
    H = rng.standard_normal((10, 6, 6))
    H = H + np.transpose(H, (0, 2, 1))
    w_budget, v_budget = eigh_batch(H, chunk_bytes=1)
    w_full, v_full = eigh_batch(H, chunk=10)
    assert np.allclose(w_budget, w_full)
    assert np.allclose(np.abs(v_budget), np.abs(v_full))


def test_eigh_batch_promotes_int_input_instead_of_silently_zeroing_it():
    """Ruling 4: `v = np.empty(H.shape, dtype=H.dtype)` silently returns
    int-truncated (all-zero, since eigenvector components are non-integer)
    eigenvectors for an int-typed H. np.result_type(H.dtype, np.float64)
    promotes int/bool to float64 (complex input stays complex) and the
    eigenvectors come out correct."""
    H = np.array([[[2, 1], [1, 2]]], dtype=np.int64)
    w, v = eigh_batch(H)
    assert v.dtype.kind == "f"
    wn, vn = np.linalg.eigh(H[0].astype(float))
    assert np.allclose(w[0], wn)
    assert np.allclose(np.abs(v[0]), np.abs(vn))


def test_eigenvectors_are_columns_and_reconstruct_H():
    """v[n][:, k] is eigenvector k: v @ diag(w) @ v.T reconstructs H exactly."""
    rng = np.random.default_rng(8)
    H = rng.standard_normal((3, 4, 4))
    H = H + np.transpose(H, (0, 2, 1))
    w, v = eigh_batch(H)
    for n in range(3):
        assert np.allclose(v[n] @ np.diag(w[n]) @ v[n].T, H[n])


def test_sweep_agrees_with_point_by_point_assembly(tm):
    """The chunked batched sweep (chunk=3) matches per-point hamiltonian() + eigvalsh."""
    pset = thf_v1()
    E = np.linspace(0.0, 60.0, 7)
    res = sweep(tm, pset, {"E_z": E, "B_z": np.zeros_like(E)}, chunk=3)
    d = len(tm.kets)  # tm.mats is a tuple of per-term (d, d) arrays, not a stacked array
    assert res.evals.shape == (7, d)
    assert res.evecs.shape == (7, d, d)
    for n, e in enumerate(E):
        ref = np.linalg.eigvalsh(hamiltonian(tm, pset, {"E_z": e, "B_z": 0.0}))
        assert np.allclose(res.evals[n], ref, atol=1e-9)


def test_sweep_result_carries_the_label_convention_and_the_active_terms(tm):
    """C2V's hardest-won lesson (wart #9): the label convention must travel
    INSIDE the artifact, or a dataset repaired under one convention silently
    invalidates a catalogue built on the other."""
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(3), "B_z": np.linspace(0.5, 1.0, 3)})
    assert res.order == "energy" and res.gauge == "none"
    assert "zeeman_Gpar" in res.active_terms and "stark_z" not in res.active_terms
    assert res.manifest["conventions"]["zeeman_sign"] == "plus_Gpar"


def test_sweep_result_carries_the_assignment_strategy(tm):
    """Ruling 2: the assignment strategy used for tracking travels on the artifact too."""
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(2), "B_z": np.zeros(2)}, assignment="hungarian")
    assert res.assignment == "hungarian"


def test_sweep_result_carries_the_zero_field_reference(tm):
    """Ruling 5: `reference` (default 0) travels on the artifact so a reader
    knows which grid index adiabatic_zero_field anchored to."""
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(2), "B_z": np.zeros(2)})
    assert res.reference == 0
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(2), "B_z": np.zeros(2)}, reference=1)
    assert res.reference == 1


def test_sweep_with_pinned_gauge_is_reproducible(tm):
    """gauge='pinned' gives byte-reproducible eigenvectors across repeated calls."""
    kw = dict(gauge="pinned")
    a = sweep(tm, thf_v1(), {"E_z": np.linspace(0, 60, 5), "B_z": np.zeros(5)}, **kw)
    b = sweep(tm, thf_v1(), {"E_z": np.linspace(0, 60, 5), "B_z": np.zeros(5)}, **kw)
    assert np.array_equal(a.evecs, b.evecs)
    from heff import track
    assert np.all(track.pin_col(np.transpose(a.evecs, (0, 2, 1))) == 1)
