"""Assignment, sign gauges and ordering policies.

The kernels are lifted from C2V-Molecules atm_core/matching/_utils.py
(C:/Users/Arian/Code/C2V-Molecules @ 59a067a): `_assign` (:265), `global_sign_col`
(:456), `pin_col` (:467) and `_PIN_TIE = 1e-6` (:453). The C2V grid machinery
(snake propagation, the 967-line repair layer) is deliberately NOT lifted --
v1 runs energy order on small blocks, and repair.py exists to certify adiabatic
labels on a 2D grid, which nothing here needs yet (spec S3.5, risk 5).

Three orthogonal axes, none forbidding another (spec S3.5):
  order ~ {energy, adiabatic_step, adiabatic_zero_field}
  gauge ~ {none, pinned, global}
  assignment strategy ~ {argmax, hungarian, adaptive}
`adiabatic_zero_field` is the thesis's stated rule (p.267: order by the
"adiabatically correlated free field state") and NEITHER source repo implements
it -- Molecule-Structure chains step-to-step overlap, C2V repairs an energy-order
seed.
"""
import numpy as np

# Lifted verbatim from C2V-Molecules atm_core/matching/_utils.py:453
# Form B dominant-component near-max tie tolerance (solver-stable tiebreak).
_PIN_TIE = 1e-6

_ORDERS = ("energy", "adiabatic_step", "adiabatic_zero_field")


def assign(matrix, *, mode="cost", strategy="adaptive"):
    """Adaptive Hungarian assignment (lifted from C2V _assign, :265).

    mode='cost' -> argmin on the matrix; mode='overlap' -> argmax on |matrix|.
    strategy='adaptive' tries argmin/argmax first and falls back to Hungarian
    the moment the result is not a permutation. scipy is imported here, not at
    module scope, to keep `import heff` under a second (gate A8).
    """
    work = -np.abs(np.asarray(matrix)) if mode == "overlap" else np.asarray(matrix)
    if strategy == "max_overlap":
        return np.argmin(work, axis=1).astype(int)
    if strategy in ("adaptive", "hungarian"):
        if strategy == "adaptive":
            idx = np.argmin(work, axis=1).astype(int)
            if len(set(idx.tolist())) == len(idx):
                return idx
        from scipy.optimize import linear_sum_assignment
        _, idx = linear_sum_assignment(work)
        return np.asarray(idx, dtype=int)
    raise ValueError(f"strategy must be adaptive|hungarian|max_overlap, got {strategy!r}")


def global_sign_col(ref_amp, ref_sign=1):
    """Per-point whole-set flip so a reference amplitude has sign ref_sign.

    Lifted from C2V matching/_utils.py:456. The per-point scalar commutes with
    any relabelling, so this gauge is matcher-invariant.
    """
    rs = 1 if ref_sign >= 0 else -1
    return np.where(np.asarray(ref_amp) * rs < 0, np.int8(-1), np.int8(1))


def pin_col(vecs, ref_sign=1):
    """Sign of each vector's dominant component, with the near-max tie rule.

    Lifted verbatim from C2V matching/_utils.py:467. `vecs` has the component
    axis LAST. Plain argmax(|.|) picks DIFFERENT components on CPU vs CUDA when
    two are near-tied with opposite signs; the lowest-index-near-max rule is what
    makes the pin machine-reproducible, and it also removes the np.sign(0) edge
    that can zero an eigenvector (Molecule-Structure Energy_Levels.py:312).
    """
    rs = 1 if ref_sign >= 0 else -1
    am = np.abs(vecs)
    mx = am.max(axis=-1, keepdims=True)
    dom = np.argmax(am >= mx * (1 - _PIN_TIE), axis=-1)
    amp = np.take_along_axis(vecs, dom[..., None], axis=-1)[..., 0]
    return np.where(amp * rs < 0, np.int8(-1), np.int8(1))


def apply_gauge(evecs, mode, *, ref_sign=1):
    """Return a NEW array with a sign gauge applied. Component axis LAST."""
    v = np.array(evecs, copy=True)
    if mode == "none":
        return v
    if mode == "pinned":
        return v * pin_col(v, ref_sign)[..., None]
    if mode == "global":
        return v * global_sign_col(v[..., 0], ref_sign)[..., None]
    raise ValueError(f"gauge must be none|pinned|global, got {mode!r}")


def order_states(evals, evecs, *, order, strategy="adaptive"):
    """Permutation (n_points, d) mapping output slot -> raw eigen index.

    'energy'                 identity (eigh already sorts ascending): exact by
                             construction, point-local, and the one policy that
                             cannot silently mislabel. Default.
    'adiabatic_step'         assign each point against the previous point.
    'adiabatic_zero_field'   assign every point against point 0 -- the thesis
                             rule (p.267).
    """
    evals = np.asarray(evals)
    evecs = np.asarray(evecs)
    n, d = evals.shape
    if order == "energy":
        return np.tile(np.arange(d), (n, 1))
    if order not in _ORDERS:
        raise ValueError(f"order must be one of {_ORDERS}, got {order!r}")
    perm = np.zeros((n, d), dtype=int)
    perm[0] = np.arange(d)
    ref = evecs[0]
    for i in range(1, n):
        prev = evecs[i - 1][:, perm[i - 1]] if order == "adiabatic_step" else ref
        overlap = prev.T @ evecs[i]                 # (slot, raw)
        perm[i] = assign(overlap, mode="overlap", strategy=strategy)
    return perm
