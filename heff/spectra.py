"""E1 transition dipoles and line strengths within one electronic state.

Follows the thesis procedure exactly (thesis digest S5): diagonalise the two
manifolds separately, build the transition-dipole matrix in the common
primitive basis, sandwich it between the two eigenvector sets, and take line
positions from eigenvalue differences.

AMPLITUDES ARE SUMMED AND THEN SQUARED (thesis p.161), so intensity-borrowing
paths interfere. That is a one-line ordering choice with a factor-of-anything
consequence, so it is encoded once, here, and gated by B8.

Within X the E1 operator is the same molecule-frame dipole as the Stark term,
evaluated at p = +-1, 0 -- so this module calls heff.elements_c.dipole_geometry
rather than re-deriving anything. A transition connects different m_F, hence
different blocks: pass the two blocks' kets and eigenvectors (spec S3.6, task
brief controller ruling 3) -- dipole_matrix indexes purely by POSITION in the
kets_a/kets_b arrays the caller passes, never by any assumed ordering in the
full basis, so a caller's Blocking.index selection is what places kets.

`dipole_matrix` returns dimensionless geometry (units of d_mf); `line_strengths`
therefore returns strengths in units of d_mf^2 -- multiply by d_mf^2 (params.py
thf_v1()['d_mf'].canonical, MHz/(V/cm)) squared, or by
params.DEBYE_TO_MHZ_PER_V_CM^2 if starting from Debye, to get a strength with
units of (MHz/(V/cm))^2. Line positions are in MHz because evals_a/evals_b (the
SweepResult.evals or a bare np.linalg.eigh output) are in MHz throughout heff.

Port in spirit, not verbatim, of C:/Users/Arian/Code/Molecule-Structure
Jupyter Notebooks/RaX/gen_spectra.py (Molecule-Structure HEAD 9eec91a):
`_strength_matrix`'s "one primitive, sum |<g|T^1_p(d)|e>|^2" idea and the
matrix-sandwich-then-square shape carry over; the coherent sum-then-square
ACROSS p here is the heff-specific generalisation the design spec calls for
(spec S3.6), which is stronger than gen_spectra.py's per-polarisation
incoherent |T_p|^2 sum (fine there because within one signed-m_F block only
one polarisation ever connects a given bra/ket pair anyway).
"""
import numpy as np

from .elements_c import dipole_geometry


def dipole_matrix(kets_a, kets_b, ctx, p):
    """<a|d_p|b> in the primitive basis, shape (len(a), len(b)).

    Dimensionless geometry; multiply by d_mf for a dipole in MHz/(V/cm).
    """
    out = np.zeros((len(kets_a), len(kets_b)))
    for i in range(len(kets_a)):
        for j in range(len(kets_b)):
            out[i, j] = dipole_geometry(kets_a[i], kets_b[j], ctx.I, p)
    return out


def _strengths_from_matrices(evals_a, evecs_a, evals_b, evecs_b, mats, *, weights=None):
    """Sum the amplitudes over polarisation, THEN square (gate B8).

    evecs_a: (d_a, n_a) with eigenvector k in column k; mats: {p: (d_a, d_b)}.
    """
    amp = None
    for p, D in mats.items():
        c = 1.0 if weights is None else float(weights.get(p, 0.0))
        if c == 0.0:
            continue
        term = c * (np.conj(evecs_a).T @ D @ evecs_b)
        amp = term if amp is None else amp + term
    if amp is None:
        amp = np.zeros((np.shape(evecs_a)[1], np.shape(evecs_b)[1]))
    freqs = np.asarray(evals_b)[None, :] - np.asarray(evals_a)[:, None]
    return freqs, np.abs(amp) ** 2


def line_strengths(evals_a, evecs_a, kets_a, evals_b, evecs_b, kets_b, ctx, *,
                   polarizations=(-1, 0, 1), weights=None):
    """Line positions (MHz) and strengths (units of d_mf^2) between two
    separately diagonalised blocks.

    weights: {p: complex} polarisation amplitudes (default 1 for each p in
    `polarizations`). Because amplitudes are summed before squaring, a coherent
    superposition of polarisations interferes -- which is the point.
    """
    if weights is None:
        weights = {p: 1.0 for p in polarizations}
    mats = {p: dipole_matrix(kets_a, kets_b, ctx, p) for p in polarizations}
    return _strengths_from_matrices(evals_a, evecs_a, evals_b, evecs_b, mats,
                                    weights=weights)
