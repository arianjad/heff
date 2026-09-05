"""Term-matrix cache and Hamiltonian assembly.

The cache IS the architecture. Molecule-Structure measured 0.88 s to build a
dim-64 Hamiltonian with ~14 terms -- roughly 57 k sympy-Wigner evaluations. The
ThF+ slice at dim 96 with 11 terms is comparable, PER PARAMETER SET; ten
thousand trial sets would be hours of arithmetic that should be a dot product.
Here it is one build plus 10^4 resums (spec S3.3).

The sweep is a single tensor contraction: stack coefficients as c[n_sets,
n_terms] and term matrices as M[n_terms, d, d], then
H_batch = np.tensordot(c, M, axes=1).

Raising is reserved for STRUCTURE -- a term matrix whose shape disagrees with
the block, a missing term, a non-Hermitian matrix from a term declared
Hermitian. Nothing raises on a hash, ever; the manifest is a record: v1 has no
disk cache backend, so there is no fingerprint-mismatch code path here at all
(the key is designed now for that future one-function backend swap) --
task-7-report.md records this as the reconciliation for controller ruling 4.
"""
import json
import time
from dataclasses import dataclass

import numpy as np

from .terms import REGISTRY, terms_for_case


@dataclass(frozen=True)
class TermMatrices:
    """Parameter-free matrices for one block, plus the manifest that travels with them."""
    names: tuple
    params: tuple
    mats: np.ndarray
    kets: np.ndarray
    manifest: dict


def build_term_matrices(kets, ctx, *, case="c", registry=REGISTRY):
    """Evaluate every applicable term once over one block.

    The declared selection rules are used as a sparsity mask so only allowed
    (i, j) are evaluated -- the dominant build cost (spec S3.2).
    """
    terms = terms_for_case(case, registry=registry)
    if not terms:
        raise ValueError(f"no terms registered for case {case!r}; import heff.elements_c")
    d = len(kets)
    t0 = time.perf_counter()
    mats = []
    for t in terms:
        M = np.zeros((d, d), dtype=float if t.real else complex)
        for i in range(d):
            for j in range(d):
                if not t.rules.allows(kets[i], kets[j]):
                    continue
                M[i, j] = t.fn(kets[i], kets[j], ctx)
        if t.hermitian and not np.allclose(M, M.conj().T, atol=1e-10):
            raise ValueError(f"term {t.name!r} is declared hermitian but its matrix is not")
        mats.append(M)
    dtype = np.result_type(*[M.dtype for M in mats])
    stack = np.array([M.astype(dtype) for M in mats])
    manifest = {
        "case": case,
        "dimension": d,
        "n_terms": len(terms),
        "terms": [t.name for t in terms],
        "cites": {t.name: t.cite for t in terms},
        "conventions": ctx.conventions.stamp(),
        "wigner_backend": "sympy+lru_cache",
        "build_seconds": time.perf_counter() - t0,
    }
    # The cache KEY is designed now even though v1 has no disk backend: canonical
    # JSON of (case, dimension, sorted term names, conventions version), so the
    # backend is a one-function swap later (spec S3.3). It is a RECORD -- nothing
    # gates on it, a mismatch would warn and rebuild.
    manifest["key"] = json.dumps(
        {"case": case, "dimension": d, "terms": sorted(manifest["terms"]),
         "conventions": manifest["conventions"]}, sort_keys=True)
    return TermMatrices(names=tuple(t.name for t in terms),
                        params=tuple(t.param for t in terms),
                        mats=stack, kets=kets, manifest=manifest)


def _knob(symbol, pset, knobs):
    """A knob's scalar value: `knobs` first, then the ParamSet, then 0.0.

    Defaulting to 0 is what makes the PT-odd block opt-in and what lets the
    thesis's own H_K = 0 deperturbation study work -- and it is why `active()`
    exists, so a term that is off is reported rather than silently absent.
    """
    if knobs is not None and symbol in knobs:
        return float(knobs[symbol])
    return pset.value(symbol, default=0.0)


def coefficients(tm, pset, knobs):
    """The coefficient of every term: the product of its knob symbols."""
    out = np.empty(len(tm.names))
    for k, symbols in enumerate(tm.params):
        val = 1.0
        for s in symbols:
            val *= _knob(s, pset, knobs)
        out[k] = val
    return out


def active(tm, pset, knobs):
    """Names of the terms with a non-zero coefficient, in stack order."""
    c = coefficients(tm, pset, knobs)
    return tuple(n for n, v in zip(tm.names, c) if v != 0.0)


def hamiltonian(tm, pset, knobs):
    """H = sum_k c_k M_k for one parameter set and one field point."""
    c = coefficients(tm, pset, knobs)
    return np.tensordot(c, tm.mats, axes=1)


def sweep_coefficients(tm, pset, knob_arrays):
    """c[n_sets, n_terms] by broadcasting, not by re-reading records (spec S3.4)."""
    arrays = {k: np.asarray(v, dtype=float) for k, v in knob_arrays.items()}
    if arrays:
        shapes = np.broadcast_shapes(*[a.shape for a in arrays.values()])
        arrays = {k: np.broadcast_to(v, shapes).ravel() for k, v in arrays.items()}
        n_sets = int(np.prod(shapes))
    else:
        n_sets = 1
    out = np.empty((n_sets, len(tm.names)))
    for k, symbols in enumerate(tm.params):
        val = np.ones(n_sets)
        for s in symbols:
            val = val * (arrays[s] if s in arrays else pset.value(s, default=0.0))
        out[:, k] = val
    return out


def hamiltonian_batch(tm, c):
    """(n_sets, d, d) in one BLAS call. Chunking is the caller's job (heff.engine)."""
    c = np.atleast_2d(np.asarray(c, dtype=float))
    # ponytail: dense assembled H; if dim > ~5e3 the sum itself needs chunking.
    return np.tensordot(c, tm.mats, axes=1)


def vertex(tm, pset, knobs, knob):
    """Exact dH/d(knob) at the given point.

    Free, because dH/dc_k = M_k is already in the catalogue: a knob's vertex is
    the sum over the terms containing it of (product of the OTHER knobs) x M_k.
    This is what makes g-factors and dipoles exact derivatives rather than
    finite differences (spec S3.6), and it is the analytic fit Jacobian for
    Hamiltonian parameters at zero extra cost.
    """
    out = np.zeros(tm.mats.shape[1:], dtype=tm.mats.dtype)
    for k, symbols in enumerate(tm.params):
        if knob not in symbols:
            continue
        val = 1.0
        for s in symbols:
            if s == knob:
                continue
            val *= _knob(s, pset, knobs)
        out = out + val * tm.mats[k]
    return out
