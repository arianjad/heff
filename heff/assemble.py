"""Build reusable term matrices and assemble H = sum_k c_k M_k."""
import hashlib
import json
import time
from dataclasses import dataclass

import numpy as np

from .terms import REGISTRY, terms_for_case


@dataclass(frozen=True)
class TermMatrices:
    """Parameter-free matrices and provenance for one block."""
    names: tuple
    params: tuple
    mats: tuple
    kets: np.ndarray
    manifest: dict


def build_term_matrices(kets, ctx, *, case="c", registry=REGISTRY, term_names=None):
    """Evaluate each applicable term once, masking by its declared rules."""
    terms = terms_for_case(case, names=term_names, registry=registry)
    if not terms:
        fix = ("pass registry=heff.elements_c2.REGISTRY_C2" if case == "c2"
               else "import heff.elements_c")
        raise ValueError(f"no terms registered for case {case!r}; {fix}")
    d = len(kets)
    single_mF = len({float(v) for v in np.asarray(kets["mF"], dtype=float)}) == 1
    if single_mF:
        for t in terms:
            if any(float(x) != 0.0 for x in t.rules.dmF):
                raise ValueError(
                    f"term {t.name!r} declares dmF={tuple(t.rules.dmF)}, which couples "
                    f"Delta m_F != 0, but this block holds the single m_F = "
                    f"{float(kets['mF'][0])}; build it on the merged full basis "
                    f"(Blocking.merge(...) / StateSpec(M='all')) instead")
    t0 = time.perf_counter()
    mats = []
    for t in terms:
        M = np.zeros((d, d), dtype=float if t.real else complex)
        for i in range(d):
            for j in range(d):
                if not t.rules.allows(kets[i], kets[j]):
                    continue
                M[i, j] = t.fn(kets[i], kets[j], ctx)
        if t.hermitian and not np.allclose(M, M.conj().T, atol=1e-10, rtol=0):
            raise ValueError(f"term {t.name!r} is declared hermitian but its matrix is not")
        mats.append(M)
    stack = tuple(mats)
    try:
        import sympy
        wigner_version = sympy.__version__
    except ImportError:
        wigner_version = "unavailable"
    manifest = {
        "case": case,
        "dimension": d,
        "n_terms": len(terms),
        "terms": [t.name for t in terms],
        "cites": {t.name: t.cite for t in terms},
        "conventions": ctx.conventions.stamp(),
        "frame": ctx.frame,
        "wigner_backend": "sympy+lru_cache",
        "wigner_version": wigner_version,
        "build_seconds": time.perf_counter() - t0,
    }
    ket_hash = hashlib.sha256(np.ascontiguousarray(kets).tobytes()).hexdigest()
    spec_desc = {"case": case, "dimension": d, "terms": sorted(manifest["terms"]),
                 "conventions": manifest["conventions"], "ket_hash": ket_hash}
    manifest["spec_hash"] = hashlib.sha256(
        json.dumps(spec_desc, sort_keys=True).encode()).hexdigest()
    return TermMatrices(names=tuple(t.name for t in terms),
                        params=tuple(t.param for t in terms),
                        mats=stack, kets=kets, manifest=manifest)


def _known_knob_symbols(tm):
    """The union of every registered term's knob symbols for this TermMatrices."""
    return frozenset(s for symbols in tm.params for s in symbols)


def _validate_knobs(tm, knobs):
    """Reject knob/override symbols that no registered term uses."""
    if not knobs:
        return
    known = _known_knob_symbols(tm)
    unknown = sorted(set(knobs) - known)
    if unknown:
        raise ValueError(
            f"unknown knob symbol(s) {unknown} for this TermMatrices; "
            f"known symbols: {sorted(known)}")


def _knob(symbol, pset, knobs):
    """Return a knob from overrides, ParamSet, or zero."""
    if knobs is not None and symbol in knobs:
        return float(knobs[symbol])
    return pset.value(symbol, default=0.0)


def coefficients(tm, pset, knobs):
    """The coefficient of every term: the product of its knob symbols."""
    _validate_knobs(tm, knobs)
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
    """Assemble H = sum_k c_k M_k with active-term dtype promotion."""
    return hamiltonian_batch(tm, coefficients(tm, pset, knobs))[0]


def sweep_coefficients(tm, pset, knob_arrays):
    """c[n_sets, n_terms] by broadcasting, not by re-reading records (spec S3.4)."""
    _validate_knobs(tm, knob_arrays)
    arrays = {k: np.asarray(v, dtype=float) for k, v in knob_arrays.items()}
    if arrays:
        try:
            shapes = np.broadcast_shapes(*[a.shape for a in arrays.values()])
        except ValueError as e:
            detail = ", ".join(f"{k!r}: shape {v.shape}" for k, v in arrays.items())
            raise ValueError(
                f"knob arrays do not broadcast to a common shape ({detail}): {e}") from e
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
    """Return ``(n_sets, d, d)`` in one contraction; callers choose chunking."""
    c = np.atleast_2d(np.asarray(c))
    if c.shape[-1] != len(tm.names):
        raise ValueError(
            f"coefficient array has c.shape[-1] = {c.shape[-1]} but this "
            f"TermMatrices has len(tm.names) = {len(tm.names)} terms")
    # ponytail: dense assembled H; if dim > ~5e3 the sum itself needs chunking.
    d = tm.mats[0].shape[0]
    active_idx = [k for k in range(c.shape[1]) if np.any(c[:, k] != 0.0)]
    if not active_idx:
        return np.zeros((c.shape[0], d, d), dtype=float)
    stack = np.array([tm.mats[k] for k in active_idx])
    return np.tensordot(c[:, active_idx], stack, axes=1)


def vertex(tm, pset, knobs, knob):
    """Return exact dH/d(knob) from cached term matrices."""
    _validate_knobs(tm, knobs)
    known = _known_knob_symbols(tm)
    if knob not in known:
        raise ValueError(
            f"unknown knob symbol {knob!r} for this TermMatrices; "
            f"known symbols: {sorted(known)}")
    contributions = []
    for k, symbols in enumerate(tm.params):
        if knob not in symbols:
            continue
        val = 1.0
        for s in symbols:
            if s == knob:
                continue
            val *= _knob(s, pset, knobs)
        if val != 0.0:
            contributions.append((val, tm.mats[k]))
    d = tm.mats[0].shape[0]
    if not contributions:
        return np.zeros((d, d), dtype=float)
    dtype = complex if any(np.iscomplexobj(M) for _, M in contributions) else float
    out = np.zeros((d, d), dtype=dtype)
    for val, M in contributions:
        out = out + val * M
    return out
