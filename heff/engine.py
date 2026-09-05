"""Batched diagonalisation over a sweep, returning eigenvectors.

numpy/scipy only. At dim <= 16 per m_F block over a few thousand points a GPU is
pure overhead (spec S3.5); the torch path is a later lift from C2V
atm_core/config.py:96 _get_eigh, and 'syevj' stays out of it (30 % slower there,
and it once shipped transposed eigenvectors -- C2V wart #6).

Controller ruling 1 (deferred from Task 7, spec S3.3): chunking is on a memory
budget, not a hardcoded point count. `chunk_bytes` estimates n_points-per-chunk
from d, dtype and the fact that eigenvectors are always retained here (ruling
3 -- unlike C2V, both sweep entry points return them); `chunk=` still overrides
the point count directly for callers/tests that want an exact chunk size.
"""
from dataclasses import dataclass

import numpy as np

from .assemble import hamiltonian_batch, sweep_coefficients
from .track import apply_gauge, order_states

# ponytail: named hook for a future backend swap (spec S3.5 lifts C2V's
# config.py:96 _get_eigh when a polyatomic basis needs GPU eigh). v1 stays
# numpy/scipy only -- no torch import anywhere in this module.
_EIGH = np.linalg.eigh

_DEFAULT_CHUNK_BYTES = 200_000_000  # a few hundred MB (ruling 1)


@dataclass(frozen=True)
class SweepResult:
    """Every result carries its label convention and its active term list.

    C2V's hardest lesson (wart #9): energy-order and character labels disagree
    at ~3.7 % of magic nodes, and repairing a dataset invalidates a catalogue
    built on the other convention. So order=, gauge=, assignment= and the
    active terms travel inside the artifact (spec S3.5; ruling 2). `manifest`
    already carries the conventions stamp (assemble.build_term_matrices).
    """
    knobs: dict
    evals: np.ndarray
    evecs: np.ndarray
    kets: np.ndarray
    order: str
    gauge: str
    assignment: str
    active_terms: tuple
    manifest: dict


def _chunk_size(n, d, dtype, chunk_bytes):
    """Points per chunk that keep one chunk's H-slice + eigenvector output
    under `chunk_bytes` (ruling 1). Eigenvectors are always materialised
    (ruling 3), so both count toward the budget."""
    itemsize = np.dtype(dtype).itemsize
    per_point = 2 * d * d * itemsize
    return max(1, min(n, chunk_bytes // per_point))


def eigh_batch(H, *, chunk=None, chunk_bytes=_DEFAULT_CHUNK_BYTES):
    """Chunked batched eigh. Returns (w (n, d) ascending, v (n, d, d) columns).

    `chunk=None` (default) derives a chunk size from `chunk_bytes`; pass an
    explicit `chunk` to force a point count (what the gate tests below do).
    """
    H = np.asarray(H)
    n, d, _ = H.shape
    if chunk is None:
        chunk = _chunk_size(n, d, H.dtype, chunk_bytes)
    w = np.empty((n, d), dtype=float)
    v = np.empty(H.shape, dtype=H.dtype)
    for c0 in range(0, n, chunk):
        c1 = min(n, c0 + chunk)
        w[c0:c1], v[c0:c1] = _EIGH(H[c0:c1])
    return w, v


def sweep(tm, pset, knob_arrays, *, chunk=None, chunk_bytes=_DEFAULT_CHUNK_BYTES,
          order="energy", gauge="none", assignment="adaptive"):
    """Diagonalise one block over a collection of knob values.

    knob_arrays values are broadcast against each other, so a 1D sweep is
    {'E_z': E, 'B_z': zeros}, and a 2D grid is the ravel of a meshgrid (ruling
    3: both return eigenvalues AND eigenvectors -- this is the one thing C2V's
    equivalent entry point does not do).
    """
    c = sweep_coefficients(tm, pset, knob_arrays)
    w, v = eigh_batch(hamiltonian_batch(tm, c), chunk=chunk, chunk_bytes=chunk_bytes)
    perm = order_states(w, v, order=order, strategy=assignment)
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
        assignment=assignment, active_terms=active_terms, manifest=dict(tm.manifest))
