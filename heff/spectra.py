"""E1 transition dipoles and line strengths within one electronic state.

Amplitudes sum before squaring, so paths interfere. Geometry is dimensionless
(units of d_mf); line positions are MHz. Adapted from
``Molecule-Structure/Jupyter Notebooks/RaX/gen_spectra.py`` at ``9eec91a``.
"""
import numpy as np

from .conventions import ef_label, parity_operator
from .elements_c import dipole_geometry


def _v1_geometry(bra, ket, ctx, p):
    """The v1 default, wrapped to `geometry`'s (bra, ket, ctx, p) shape."""
    return dipole_geometry(bra, ket, ctx.I, p)


def dipole_matrix(kets_a, kets_b, ctx, p, *, geometry=None):
    """Return dimensionless <a|d_p|b> geometry; optional ``geometry`` has v2 shape.

    Bind ``elements_c2.axial_geometry`` at ``k=1, q=0`` for v2; it reduces to
    the v1 expression at I_Th=0 ([HAM] S9.1).
    """
    fn = _v1_geometry if geometry is None else geometry
    out = np.zeros((len(kets_a), len(kets_b)))
    for i in range(len(kets_a)):
        for j in range(len(kets_b)):
            out[i, j] = fn(kets_a[i], kets_b[j], ctx, p=p)
    return out


def _strengths_from_matrices(evals_a, evecs_a, evals_b, evecs_b, mats, *, weights=None):
    """Sum the amplitudes over polarisation, then square.

    evecs_a: (d_a, n_a) with eigenvector k in column k; mats: {p: (d_a, d_b)}.
    weights are cast with `complex(...)` (a real weight, e.g. from the
    `weights=None` default, is just a zero-imaginary complex) and a channel is
    skipped when `abs(weight) == 0`, not merely `weight == 0.0`, so a weight
    given as `0j` is dropped too. A strength is |amplitude|^2, which is real by
    construction no matter how complex the weights or eigenvectors are; the
    explicit `.astype(float)` below only documents that, it changes no value.
    """
    amp = None
    for p, D in mats.items():
        c = 1.0 if weights is None else complex(weights.get(p, 0.0))
        if abs(c) == 0.0:
            continue
        term = c * (np.conj(evecs_a).T @ D @ evecs_b)
        amp = term if amp is None else amp + term
    if amp is None:
        amp = np.zeros((np.shape(evecs_a)[1], np.shape(evecs_b)[1]))
    freqs = np.asarray(evals_b)[None, :] - np.asarray(evals_a)[:, None]
    strengths = (np.abs(amp) ** 2).astype(float)
    return freqs, strengths


def line_strengths(evals_a, evecs_a, kets_a, evals_b, evecs_b, kets_b, ctx, *,
                   polarizations=(-1, 0, 1), weights=None, geometry=None):
    """Line positions (MHz) and strengths (units of d_mf^2) between two
    separately diagonalised blocks.

    weights: {p: complex} polarisation amplitudes (default 1 for each p in
    `polarizations`). Because amplitudes are summed before squaring, a coherent
    superposition of polarisations interferes -- which is the point.

    `geometry` is forwarded to `dipole_matrix` unchanged (see its docstring);
    `None` keeps the v1 formula, a v2 caller passes
    `functools.partial(heff.elements_c2.axial_geometry, k=1, q=0.0)`.

    CAVEAT on the default: `polarizations=(-1, 0, 1)` with unit real weights
    adds all three channels COHERENTLY before squaring. That default is only
    physically meaningful when both `kets_a`/`kets_b` blocks have definite m_F
    (e.g. one output of `block_by_mF` on each side) -- then the 3j selection
    rule `m_bra = m_ket + p` already forces exactly one p to be non-zero for
    any given pair of states, so the "coherent sum" is really a sum of one
    non-zero term and two exact zeros, and no cross-polarisation interference
    actually happens. Feeding two blocks that mix multiple m_F (e.g. the full
    basis on both sides) with the default weights would coherently add
    physically distinct polarisation channels that do not interfere in
    reality; pass explicit `weights` (and generally a single `polarizations`
    value per physical channel) in that case.
    """
    if weights is None:
        weights = {p: 1.0 for p in polarizations}
    mats = {p: dipole_matrix(kets_a, kets_b, ctx, p, geometry=geometry) for p in polarizations}
    return _strengths_from_matrices(evals_a, evecs_a, evals_b, evecs_b, mats,
                                    weights=weights)


def label_lines(kets, evecs, S, *, rule, ell, s):
    """One label per eigenvector column: its dominant (J, F) and the parity of
    the +-Omega superposition it sits in, as an 'e'/'f' name under `rule`.

    Returns a list of dicts, one per column of `evecs` (same order):
    `{"J": float, "F": float, "parity": +-1, "ef": "e" or "f"}`. A line can
    then be reported as e.g. "(J=1, F=3/2, +) -> (J=2, F=5/2, -)" by pairing a
    row of `label_lines(kets_a, evecs_a, ...)` with a column of
    `label_lines(kets_b, evecs_b, ...)`.

    When `kets.dtype` carries an `F1` field (the two-spin `KET_C2` basis) the
    dict also gets `"F1": float` (the SAME dominant component's F1, so it is
    read off the identical `dom` index as J and F, not a separately-weighted
    group) and `"purity": float` -- |amplitude|^2 of that dominant component,
    i.e. how much of the eigenvector actually sits on the labelled (J, F1, F)
    ket.

    READ `purity` AGAINST 0.5, NOT 1.0, at zero field. A parity eigenstate of
    this molecule is a +-Omega DOUBLET, (|+Om> +- |-Om>)/sqrt(2), and the two
    halves are separate kets of the basis, so a perfectly clean, perfectly
    assigned state has purity 0.5 exactly. Values near 0.5 mean a well-defined
    (J, F1, F); values well below it mean the label is a dominant-component
    assignment over a mixed state. A `KET_C` basis (no `F1` field) gets neither
    the `F1` nor the `purity` key.

    Labels are recorded on the output for a caller to report; `line_strengths`
    itself never calls this and a line's strength/frequency never depends on
    it -- labelling only names states after the transition calculation.

    The superposition parity comes from `heff.conventions.parity_operator`,
    which is built from the same +-Omega ket pairing used by
    `heff.conventions.superposition_parity`: `P` is applied to each
    eigenvector and its parity is read off as the (rounded) expectation value
    `<v|P|v>`, which is the exact eigenvalue (+1 or -1) whenever `v` really is
    a parity eigenstate (guaranteed at zero field, spec [HAM] V8) and a
    nearest-parity label otherwise. `S`, `ell`, `s` have no default here, same
    as in `heff.conventions` -- for ThF+ X 3Delta1 pass `S=ctx.S` (Ctx already
    carries the electronic S, spec S3.2), `ell=0.0`, `s=0.0`.

    Each eigenvector is labelled by the (J, F) of its largest-|amplitude|
    basis component; a state that is not dominated by one (J, F) still gets a
    label, just a less meaningful one -- this is a display convenience, not a
    physics claim.
    """
    P = parity_operator(kets, S, ell=ell, s=s)
    has_F1 = "F1" in kets.dtype.names
    out = []
    for k in range(evecs.shape[1]):
        v = evecs[:, k]
        dom = int(np.argmax(np.abs(v)))
        J, F = float(kets["J"][dom]), float(kets["F"][dom])
        parity = int(round(float(np.real(np.vdot(v, P @ v)))))
        row = {"J": J, "F": F, "parity": parity,
              "ef": ef_label(J, parity, rule=rule, S=S, ell=ell)}
        if has_F1:
            row["F1"] = float(kets["F1"][dom])
            row["purity"] = float(np.abs(v[dom]) ** 2)
        out.append(row)
    return out
