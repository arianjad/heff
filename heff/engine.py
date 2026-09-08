"""Batched CPU diagonalisation with eigenvalues and eigenvectors retained.

`chunk_bytes` estimates the chunk size from matrix dimension and dtype;
`chunk` overrides it with an explicit point count. Result metadata records
state ordering, gauge, assignment strategy, and the tracking reference.
"""
from dataclasses import dataclass

import numpy as np

from .assemble import hamiltonian_batch, sweep_coefficients
from .track import apply_gauge, order_states

# Backend hook retained for the planned PyTorch implementation.
_EIGH = np.linalg.eigh

_DEFAULT_CHUNK_BYTES = 200_000_000


@dataclass(frozen=True)
class SweepResult:
    """Sweep arrays, label conventions, active terms, and matrix provenance.

    `reference` records the grid index used as the zero-field anchor by
    `order='adiabatic_zero_field'`; the sweep need not start at zero field.
    """
    knobs: dict
    evals: np.ndarray
    evecs: np.ndarray
    kets: np.ndarray
    order: str
    gauge: str
    assignment: str
    reference: int
    active_terms: tuple
    manifest: dict


def _chunk_size(n, d, dtype, chunk_bytes):
    """Estimate points per chunk from H and eigenvector storage, minimum one."""
    itemsize = np.dtype(dtype).itemsize
    per_point = 2 * d * d * itemsize
    return max(1, min(n, chunk_bytes // per_point))


def eigh_batch(H, *, chunk=None, chunk_bytes=_DEFAULT_CHUNK_BYTES):
    """Chunked batched eigh. Returns (w (n, d) ascending, v (n, d, d) columns).

    `chunk=None` (default) derives a chunk size from `chunk_bytes`; pass an
    explicit `chunk` to force a point count.
    """
    H = np.asarray(H)
    n, d, _ = H.shape
    if chunk is None:
        chunk = _chunk_size(n, d, H.dtype, chunk_bytes)
    w = np.empty((n, d), dtype=float)
    # Eigenvectors need floating-point storage even for integer/bool inputs.
    v = np.empty(H.shape, dtype=np.result_type(H.dtype, np.float64))
    for c0 in range(0, n, chunk):
        c1 = min(n, c0 + chunk)
        w[c0:c1], v[c0:c1] = _EIGH(H[c0:c1])
    return w, v


def sweep(tm, pset, knob_arrays, *, chunk=None, chunk_bytes=_DEFAULT_CHUNK_BYTES,
          order="energy", gauge="none", assignment="adaptive", reference=0):
    """Diagonalise one block over a collection of knob values.

    knob_arrays values are broadcast against each other, so a 1D sweep is
    {'E_z': E, 'B_z': zeros}, and a 2D grid is flattened after broadcasting.

    `reference` (default 0) is the grid index `order=
    'adiabatic_zero_field'` treats as the zero-field point; it is ignored by
    the other two order policies. Pass it whenever the sweep's own knob_arrays
    do not start at the zero-field point, or tracking silently anchors to
    whatever point 0 happens to be.
    """
    c = sweep_coefficients(tm, pset, knob_arrays)
    w, v = eigh_batch(hamiltonian_batch(tm, c), chunk=chunk, chunk_bytes=chunk_bytes)
    perm = order_states(w, v, order=order, strategy=assignment, reference=reference)
    w = np.take_along_axis(w, perm, axis=1)
    v = np.take_along_axis(v, perm[:, None, :], axis=2)
    if gauge != "none":
        v = np.transpose(apply_gauge(np.transpose(v, (0, 2, 1)), gauge), (0, 2, 1))
    # A term is active for the WHOLE sweep if any point gives it a nonzero
    # coefficient -- same rule hamiltonian_batch itself uses to decide what to
    # sum (assemble.py); a first-point-only check would silently mislabel a
    # term that starts at zero and turns on later in the grid.
    active_terms = tuple(n for n, nz in zip(tm.names, np.any(c != 0.0, axis=0)) if nz)
    return SweepResult(
        knobs={k: np.asarray(a, dtype=float) for k, a in knob_arrays.items()},
        evals=w, evecs=v, kets=tm.kets, order=order, gauge=gauge,
        assignment=assignment, reference=reference, active_terms=active_terms,
        manifest=dict(tm.manifest))
