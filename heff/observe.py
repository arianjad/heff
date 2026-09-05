"""Observables: expectation values, exact derivatives, pair differentials.

Derivatives are exact, never finite differences. Molecule-Structure computes
g-factors and dipoles by finite difference today (Energy_Levels.py:408, :422),
which costs a step-size parameter and a re-matching call per point; the
Hellmann-Feynman / sum-over-states kernel below removes both (spec S3.6).

Convention leakage is the trap: whether Delta-g or delta-g comes out depends on
the ParamSet's dg_def, so pair_differential returns the LABEL alongside the
number, and g_factors stamps the whole conventions block on its result. Two
people comparing numbers otherwise differ by exactly 2 and neither knows why
([HAM] OPEN-12).

No physical quantity has a default here. mu_B arrives on the frozen Ctx, the
same route every matrix element takes (heff.terms.Ctx), so a g-factor cannot be
computed in the wrong unit system by omitting an argument.
"""
import numpy as np

from .assemble import hamiltonian, vertex


# ---------------------------------------------------------------------------
# Lifted VERBATIM from C2V-Molecules atm_core/qgt.py:180 multi_curvature
# Source: C:/Users/Arian/Code/C2V-Molecules @ 59a067a
# Gated there by test_multi_axis.py and test_trap_cross.py, both
# Hamiltonian-agnostic. The self-term-excluded-by-INDEX rule and the
# zero-coupling-at-zero-gap guard are both hard-won -- do not "simplify" either.
# `d2` is the full second derivative d2E/dx_a dx_b, not the second-order
# perturbation coefficient (which is half of it).
# ---------------------------------------------------------------------------
def multi_curvature(H0, verts, *, eigh=np.linalg.eigh):
    """Exact 1st + 2nd derivatives of every eigenvalue of H(x) = H0 + sum_k x_k*verts[k]
    at x = 0, for ANY named set of linear field knobs -- the n-vertex generalization of
    energy_curvature (fixed He/Hz/Ha slots; gate-checked equal on that subset,
    test_multi_axis.py).  Same HF diagonal + SOS kernel, one eigh total; exact at any
    non-degenerate point, singular only at a true degeneracy of a COUPLED state.

    verts: {name: (n,n) Hermitian}.  Returns dict(W, V, d={name: (n,)},
    d2={(a, b): (n,)}) with d2 keyed on ordered pairs, a before b in verts order."""
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
    """<psi_k|O|psi_k> for a batch of eigenvector sets.

    One einsum with a caller-supplied subscript pattern, after C2V-Molecules
    atm_core/observables.py (29 lines). Default pattern expects evecs
    (n_points, d_basis, d_states) with eigenvector k in COLUMN k -- the shape
    heff.engine.SweepResult.evecs carries.
    """
    return np.einsum(axes, np.conj(evecs), op, evecs)


def offdiag(evecs, op, i, j):
    """<psi_i|O|psi_j> per point -- the named-pair off-diagonal mode.

    Required because an imaginary-Hermitian operator has identically zero
    diagonal expectation value on real eigenvectors; Molecule-Structure
    documents exactly this at Energy_Levels.py:481-483 for its NSD-PV operator
    and directs the reader to the off-diagonal element (spec S3.6).
    """
    return np.einsum("ni,ij,nj->n", np.conj(evecs)[..., i], op, evecs[..., j])


def g_factors(tm, pset, knobs, *, ctx):
    """Exact g-factors and induced dipoles for one signed-m_F block.

    g is defined by E = -g mu_B B m_F (Ng, Petrov; [HAM] S2.8), so
    g = -(dE/dB_z) / (mu_B m_F). d_eff = -dE/dE_z in MHz/(V/cm). `d2` comes
    back too -- multi_curvature computes it in the same pass, and d2[(E_z, E_z)]
    is the polarizability.

    mu_B is read from `ctx`, the frozen element context (heff.terms.Ctx), not
    from a module default: a physical constant that can be omitted is a
    physical constant that will be wrong somewhere.

    Requires a single non-zero signed m_F, which is what the default blocking
    gives; m_F = 0 has no g in this definition and raises.
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
    """Signed difference of a pair observable, plus the label of the convention.

    The two headline ThF+ observables are differences between NAMED partners --
    the upper and lower Omega-doublet components, or +-m_F -- so this is
    required, not a convenience (spec S3.6). The caller says which partner is
    which and the sign survives; C2V's `f >= 0` fold once turned a zero-g Newton
    objective into a V-minimum with no sign change (C2V digest Q4), and that is
    the mistake this signature exists to make impossible. Works on scalars and
    arrays.
    """
    diff = np.asarray(upper) - np.asarray(lower)
    if conventions.dg_def == "delta":
        return diff / 2.0, "delta = (g^u - g^l)/2  [Ng thesis App. C.4]"
    return diff, "Delta = g^u - g^l  [Petrov arXiv:2503.02840 Eq. 18; Ng 2022 paper]"
