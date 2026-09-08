"""Sweep and batched-eigensolver contracts.

Sweeps return eigenvectors as well as eigenvalues so downstream code can
compute transition dipoles, parity labels, and expectation values.
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
    """`chunk_bytes` controls chunk sizing instead of a fixed point count.

    A budget below one point's working space still yields one point per chunk
    and agrees with one unchunked call.
    """
    rng = np.random.default_rng(99)
    H = rng.standard_normal((10, 6, 6))
    H = H + np.transpose(H, (0, 2, 1))
    w_budget, v_budget = eigh_batch(H, chunk_bytes=1)
    w_full, v_full = eigh_batch(H, chunk=10)
    assert np.allclose(w_budget, w_full)
    assert np.allclose(np.abs(v_budget), np.abs(v_full))


def test_eigh_batch_promotes_int_input_instead_of_silently_zeroing_it():
    """Integer and Boolean Hamiltonians promote to float64 for eigenvectors.

    Complex input remains complex; using the integer input dtype would truncate
    non-integer eigenvector components.
    """
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
    """The sweep manifest carries the label convention used to produce it."""
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(3), "B_z": np.linspace(0.5, 1.0, 3)})
    assert res.order == "energy" and res.gauge == "none"
    assert "zeeman_Gpar" in res.active_terms and "stark_z" not in res.active_terms
    assert res.manifest["conventions"]["zeeman_sign"] == "plus_Gpar"


def test_sweep_result_carries_the_assignment_strategy(tm):
    """The sweep result records the tracking assignment strategy."""
    res = sweep(tm, thf_v1(), {"E_z": np.zeros(2), "B_z": np.zeros(2)}, assignment="hungarian")
    assert res.assignment == "hungarian"


def test_sweep_result_carries_the_zero_field_reference(tm):
    """The sweep result records the grid index that anchors zero-field tracking."""
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
