"""Expectation values and exact Hellmann-Feynman/sum-over-states derivatives."""
import numpy as np

from .assemble import hamiltonian, vertex


# Adapted from ``C2V-Molecules/atm_core/qgt.py`` at ``59a067a``.
# ``d2`` is the full second derivative; zero coupling remains zero at a degeneracy.
def multi_curvature(H0, verts, *, eigh=np.linalg.eigh):
    """Exact derivatives of H(x)=H0+sum_k x_k verts[k] at x=0.

    Exact away from coupled degeneracies. Returns W, V, d, and ordered-pair d2.
    """
    W, V = eigh(H0)
    n = len(W)
    dW = W[:, None] - W[None, :]
    di = np.arange(n)
    dW[di, di] = np.inf                            # self-term excluded by INDEX (a9351b9)
    inv = 1.0 / dW
    A = {k: V.conj().T @ M @ V for k, M in verts.items()}

    def curv(Ai, Aj):
        num = 2.0 * np.real(Ai * np.conj(Aj))
        contrib = num * inv
        contrib[num == 0] = 0.0                    # zero coupling -> 0 even at zero gap
        return contrib.sum(axis=1)

    names = list(verts)
    d = {k: np.real(np.diag(A[k])).copy() for k in names}
    d2 = {(a, b): curv(A[a], A[b]) for i, a in enumerate(names) for b in names[i:]}
    return dict(W=W, V=V, d=d, d2=d2)


def expectation(evecs, op, *, axes="nik,ij,njk->nk"):
    """Return diagonal expectations for a batch of column-eigenvector sets."""
    return np.einsum(axes, np.conj(evecs), op, evecs)


def offdiag(evecs, op, i, j):
    """Return <psi_i|O|psi_j> per point."""
    return np.einsum("ni,ij,nj->n", np.conj(evecs)[..., i], op, evecs[..., j])


def g_factors(tm, pset, knobs, *, ctx):
    """Return exact g, d_eff [MHz/(V/cm)], and energy curvatures for one m_F block.

    Uses E = -g mu_B B m_F ([HAM] S2.8); m_F=0 is undefined.
    """
    mFs = np.unique(np.asarray(tm.kets["mF"], dtype=float))
    if len(mFs) != 1:
        raise ValueError(f"g_factors needs one signed m_F block, got {mFs}")
    m = float(mFs[0])
    if m == 0.0:
        raise ValueError("g is undefined at m_F = 0 under E = -g mu_B B m_F")
    H0 = hamiltonian(tm, pset, knobs)
    verts = {"E_z": vertex(tm, pset, knobs, "E_z"),
             "B_z": vertex(tm, pset, knobs, "B_z")}
    mc = multi_curvature(H0, verts)
    return {"W": mc["W"], "V": mc["V"],
            "g": -mc["d"]["B_z"] / (ctx.mu_B * m),
            "d_eff": -mc["d"]["E_z"],
            "d2": mc["d2"],
            "mF": m,
            "conventions": pset.conventions.stamp()}


def pair_differential(upper, lower, conventions):
    """Return a signed pair difference and its ``delta_g`` convention label."""
    diff = np.asarray(upper) - np.asarray(lower)
    if conventions.dg_def == "delta":
        return diff / 2.0, "delta = (g^u - g^l)/2  [Ng thesis App. C.4]"
    return diff, "Delta = g^u - g^l  [Petrov arXiv:2503.02840 Eq. 18; Ng 2022 paper]"
