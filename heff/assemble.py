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

The cache key is block-aware: it folds in a hash of the block's own kets (not
just case/dimension/term-names/conventions), because two blocks can share
every one of those and still hold different matrices -- e.g. the +m_F = 3/2
and -m_F = 3/2 blocks of one spec (fix round 1, finding 1).
"""
import hashlib
import json
import time
from dataclasses import dataclass

import numpy as np

from .terms import REGISTRY, terms_for_case


@dataclass(frozen=True)
class TermMatrices:
    """Parameter-free matrices for one block, plus the manifest that travels with them.

    `mats` is a tuple of per-term (d, d) arrays, EACH IN ITS OWN DTYPE (real or
    complex) -- promotion to a common dtype happens at assembly, not here, and
    only for terms actually active in that assembly (spec S3.3; fix round 1,
    finding 4).
    """
    names: tuple
    params: tuple
    mats: tuple
    kets: np.ndarray
    manifest: dict


def build_term_matrices(kets, ctx, *, case="c", registry=REGISTRY):
    """Evaluate every applicable term once over one block.

    The declared selection rules are used as a sparsity mask so only allowed
    (i, j) are evaluated -- the dominant build cost (spec S3.2).
    """
    terms = terms_for_case(case, registry=registry)
    if not terms:
        # Name the actual fix per case: the v2 two-spin terms live in their own
        # registry and are NOT reachable by importing a module (importing
        # heff.elements_c2 registers them into REGISTRY_C2, not the global one).
        fix = ("pass registry=heff.elements_c2.REGISTRY_C2" if case == "c2"
               else "import heff.elements_c")
        raise ValueError(f"no terms registered for case {case!r}; {fix}")
    d = len(kets)
    # Refuse to block a Delta-m_F != 0 term into a single-m_F block (spec S3.1
    # "what could go wrong" (iv)): the rules mask would zero every element the
    # term actually has and the caller would get a silent all-zero matrix
    # instead of a transverse/rotating-field coupling. Structural, so it raises.
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
    # Each term keeps its own dtype -- no cross-term dtype promotion at build
    # (fix round 1, finding 4; spec S3.3 "complex vs real is decided at
    # assembly, not at build").
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
    # The cache KEY is designed now even though v1 has no disk backend: canonical
    # JSON of (case, dimension, sorted term names, conventions version, a hash
    # of the block's own kets), so the backend is a one-function swap later
    # (spec S3.3). It is a RECORD -- nothing gates on it, a mismatch would warn
    # and rebuild. The ket hash is what makes the key BLOCK-aware: two blocks
    # with identical case/dimension/terms/conventions (e.g. +-m_F = 3/2) still
    # get different keys because their kets differ (fix round 1, finding 1).
    ket_hash = hashlib.sha256(np.ascontiguousarray(kets).tobytes()).hexdigest()
    spec_desc = {"case": case, "dimension": d, "terms": sorted(manifest["terms"]),
                 "conventions": manifest["conventions"], "ket_hash": ket_hash}
    manifest["spec_hash"] = hashlib.sha256(
        json.dumps(spec_desc, sort_keys=True).encode()).hexdigest()
    manifest["key"] = json.dumps(spec_desc, sort_keys=True)
    return TermMatrices(names=tuple(t.name for t in terms),
                        params=tuple(t.param for t in terms),
                        mats=stack, kets=kets, manifest=manifest)


def _known_knob_symbols(tm):
    """The union of every registered term's knob symbols for this TermMatrices."""
    return frozenset(s for symbols in tm.params for s in symbols)


def _validate_knobs(tm, knobs):
    """Raise on a knob/override symbol that no registered term uses (fix round
    1, finding 3) -- `sweep_coefficients(tm, pset, {"NOPE": ...})` used to
    return a valid array with the typo silently dropped."""
    if not knobs:
        return
    known = _known_knob_symbols(tm)
    unknown = sorted(set(knobs) - known)
    if unknown:
        raise ValueError(
            f"unknown knob symbol(s) {unknown} for this TermMatrices; "
            f"known symbols: {sorted(known)}")


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
    """H = sum_k c_k M_k for one parameter set and one field point.

    Promotion to complex happens HERE, not at build: each term matrix keeps
    its own dtype (build_term_matrices), and the assembled H is promoted to
    complex only if an ACTIVE (nonzero-coefficient) term is complex -- an
    all-real active set stays float64 (spec S3.3; fix round 1, finding 4).
    """
    c = coefficients(tm, pset, knobs)
    d = tm.mats[0].shape[0]
    active_idx = [k for k in range(len(c)) if c[k] != 0.0]
    if not active_idx:
        return np.zeros((d, d), dtype=float)
    stack = np.array([tm.mats[k] for k in active_idx])
    return np.tensordot(np.asarray(c)[active_idx], stack, axes=1)


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
    """(n_sets, d, d) in one BLAS call. Chunking is the caller's job (heff.engine).

    Same active-only promotion rule as `hamiltonian` (fix round 1, finding 4):
    a term is active for the batch if any set gives it a nonzero coefficient.
    No hard float cast on `c` -- the caller's own dtype (always real for a
    physical field/parameter sweep) is preserved.
    """
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
    """Exact dH/d(knob) at the given point.

    Free, because dH/dc_k = M_k is already in the catalogue: a knob's vertex is
    the sum over the terms containing it of (product of the OTHER knobs) x M_k.
    This is what makes g-factors and dipoles exact derivatives rather than
    finite differences (spec S3.6), and it is the analytic fit Jacobian for
    Hamiltonian parameters at zero extra cost.
    """
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
