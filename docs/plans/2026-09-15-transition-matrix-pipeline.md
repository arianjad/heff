# Transition-matrix pipeline — implementation plan

**Goal:** a generator that turns (operator, initial manifold, final manifold, field point) into a labelled complex transition matrix, first exercised on two-photon transitions in 232ThF+ X 3Delta1, J = 1..5, for every photon-polarization pair, at fixed B and along a B sweep.

**Approach:** an operator is a dict of field-free channel matrices plus per-polarization complex weights; `transition.py` owns the eigenvector sandwich (sum over channels, then square), per-signed-m_F diagonalization with a pinned gauge, manifold filters and B sweeps. The existing `twophoton.py` geometry is audited and kept; it gains the ladder reading as default. Plots and a script produce the 232 results.

**Design:** `docs/plans/2026-09-15-transition-matrix-pipeline-design.md` (agreed 2026-09-15). Audit: `docs/lit/2026-09-15-twophoton-operator-audit.md`.

**Constraints (every task inherits):**
- Env: `conda run -n structure python ...`; tests: `conda run -n structure python -m pytest -q tests/<file>`.
- Units: energies MHz, E_z V/cm, B_z G, geometry dimensionless, alphas placeholder 1.0.
- Sigma+ Jones vector is Condon-Shortley `(-1, -i, 0)/sqrt2`; sigma- is `(1, -i, 0)/sqrt2`; pi is `z`.
- Both photons absorbed (ladder): Delta m_F = p1 + p2. This replaces the Raman default everywhere.
- Transition matrices use K = 2 channels only. K = 0 is the identity and never enters a transition matrix.
- Do not touch `heff/elements_c2.py`, `heff/conventions.py`, `heff/models/thf_plus.toml`, `heff/params.py`.
- Docstrings describe the code as it is after the task; delete stale sentences, do not annotate them.
- Commit per task by explicit pathspec, message under 72 chars on line 1, trailer:
  `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01P3zj3Mfmkbtc96ieebWZdn`.
- Delegation boundary for every subagent: write only the files listed in its task; no `git push`; no process kills; checkpoint file `docs/plans/checkpoints/2026-09-15-task<N>.md` written at the first green check with what passed; stop when the task's Check passes and the commit is made.
- Task 1 first. Tasks 2 and 4 may run in parallel after Task 1. Task 3 after Task 2.

## File map

- Modify: `heff/twophoton.py` — `dyad_weights(eps1, eps2, *, reading="ladder")`; public `leg_weights`; named `POLARIZATIONS`; `jones(name_or_vector)`; docstring sweep incl. K = 0 identity and audit defects 1-3.
- Modify: `tests/test_twophoton.py` — reading-dependent assertions flipped to ladder; phantom test citation removed; normalization sentence corrected; one new test that `reading="raman"` reproduces the old table.
- Modify: `notebooks/build_isotopologues.py` S11-13 text and `TWOP_PAIRS` mapping — ladder labels; then rebuild and re-execute the notebook (Task 4).
- Create: `heff/transition.py` — `Eigensystem`, `diagonalize`, `padded_thf`, `select`, `TwoPhotonOperator`, `DipoleOperator`, `TransitionMatrix`, `transition_matrix`, `sweep_transition_matrix`.
- Create: `tests/test_transition.py` — gates 1-8 of the design.
- Create: `heff/plot_transition.py` — `heatmap_grid`, `curves_vs_B`, hierarchical ticks.
- Create: `scripts/plot_thf_twophoton.py` — the 232 run; writes `results/thf-twophoton-2026-09-15/`.
- Create: `results/thf-twophoton-2026-09-15/README.md` (by the script).
- Modify: `docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5.1(3) and S9.5.4 — ladder reading, K = 0 identity, link to audit.
- Modify: `docs/open-questions.md` OPEN-21/23 — one sentence each on what the pipeline uses.
- Modify: `docs/architecture.md` — one paragraph on `transition.py`.

## Tasks

### Task 1: ladder reading and named polarizations in `twophoton.py`

Files: `heff/twophoton.py`, `tests/test_twophoton.py`.
Interfaces produced: `dyad_weights(eps1, eps2, *, reading="ladder") -> {(K, P): complex}`; `leg_weights(eps) -> {p: complex}`; `POLARIZATIONS: dict[str, np.ndarray]`; `jones(x) -> np.ndarray`.

Steps:

1. In `heff/twophoton.py` replace `_leg` and `dyad_weights` with:

```python
def leg_weights(eps):
    """Return c[p] = (-1)^p eps_{-p}, so that d.eps = sum_p c[p] d_p (B&C 5.111)."""
    s = _spherical(eps)
    return {p: (-1.0) ** p * s[-p] for p in (-1, 0, 1)}


_SQ2 = np.sqrt(2.0)
#: Cartesian Jones vectors relative to lab z (E and B axis), Condon-Shortley.
POLARIZATIONS = {
    "sigma+": np.array([-1.0, -1.0j, 0.0]) / _SQ2,
    "sigma-": np.array([1.0, -1.0j, 0.0]) / _SQ2,
    "pi": np.array([0.0, 0.0, 1.0]),
    "x": np.array([1.0, 0.0, 0.0]),
    "y": np.array([0.0, 1.0, 0.0]),
}


def jones(x):
    """Return a unit Cartesian Jones vector from a name in POLARIZATIONS or a 3-vector."""
    v = POLARIZATIONS[x] if isinstance(x, str) else np.asarray(x, dtype=complex).reshape(3)
    n = np.linalg.norm(v)
    if n == 0.0:
        raise ValueError("polarization vector must be non-zero")
    return v / n


def dyad_weights(eps1, eps2, *, reading="ladder"):
    """Return w[(K, P)] with T_eff = sum_{K,P} w[(K,P)] alpha^K_P (B&C 5.141).

    reading="ladder": both photons absorbed, T = (d.eps2)(d.eps1)/Delta, so
    Delta m_F = p1 + p2 and (sigma+, sigma+) reaches Delta m_F = +2.
    reading="raman": the second photon is emitted, T = (d.eps2*)(d.eps1)/Delta,
    so (sigma+, sigma+) reaches Delta m_F = 0.
    Slot a (left operator factor) carries eps2, slot b carries eps1; the
    coupling is B&C 5.141 with k1 = k2 = 1. With the K = 1 part absent the
    result is symmetric under eps1 <-> eps2.
    """
    if reading not in ("ladder", "raman"):
        raise ValueError(f"reading must be 'ladder' or 'raman', got {reading!r}")
    e2 = np.asarray(jones(eps2), dtype=complex)
    if reading == "raman":
        e2 = np.conj(e2)
    c_a = leg_weights(e2)
    c_b = leg_weights(jones(eps1))
    return {(K, P): sum(_cg(1, p_a, 1, P - p_a, K, P) * c_a[p_a] * c_b[P - p_a]
                        for p_a in (-1, 0, 1) if abs(P - p_a) <= 1)
            for K in (0, 1, 2) for P in range(-K, K + 1)}
```

   Keep `_leg = leg_weights` as an alias for one release so `tests/test_twophoton.py` imports still resolve, then update that import to `leg_weights`.

2. `two_photon_line_strengths(..., eps1, eps2, alphas, reading="ladder")`: pass `reading` through to `dyad_weights`; update its docstring: drop the "RAMAN reading" paragraph, say "both photons absorbed by default; `reading='raman'` conjugates eps2".

3. Docstring sweep in `heff/twophoton.py`:
   - Module docstring: state channels (K, |dOmega|) = (0,0), (2,0), (2,2); K = 0 is exactly the identity on any basis (probe 2026-09-15, gate `test_K0_is_identity`); transitions live in K = 2; K = 1 absent (OPEN-21); ladder reading default.
   - `two_photon_K0_dOm0.cite`: append "Exactly the identity matrix: a state-independent light shift, never a transition."
   - `two_photon_K2_dOm0.cite` and `two_photon_K2_dOm2.cite`: add the Cartesian interpretation: "T^2_0 ~ (2 a_zz - a_xx - a_yy); T^2_{+-2} ~ (a_xx - a_yy). For a case (c) Omega = +-1 state the operator statement is the reduced element per q; the Cartesian form is the interpretation."
   - `_CITE`: delete the sentence "not unit ONE-photon reduced elements, which is a different (and unmade) claim" and replace by "the rank-K reduced element is unity per (K, dOmega) channel".

4. `tests/test_twophoton.py`:
   - Line 68: delete the sentence citing `test_..._fails_if_the_q_plus_2_component_changes_sign` and add that test:

```python
def test_V24_fails_if_the_q_plus_2_component_changes_sign(basis2):
    """Negative control for the shared alpha of q = +2 and q = -2."""
    kets, ctx = basis2
    P_op = parity_operator(kets, S=ctx.S, ell=0.0, s=0.0)
    for P in range(-2, 3):
        T = two_photon_matrix(kets, kets, ctx, K=2, dOmega=2, P=P)
        good = np.max(np.abs(T @ P_op - P_op @ T))
        flip = np.where((kets["Om"][:, None] - kets["Om"][None, :]) > 0, -1.0, 1.0)
        bad = np.max(np.abs((T * flip) @ P_op - P_op @ (T * flip)))
        assert good < 1e-12 and bad > 0.5, (P, good, bad)
```

   - Line 148: "unit one-photon reduced elements" -> "unit rank-K reduced element per channel".
   - In `test_dyad_weights_reproduce_the_known_polarisation_limits`: the reach table becomes ladder: `("s+","s+") -> [2]`, `("s+","s-") -> [0]`, `("s-","s+") -> [0]`, `("s-","s-") -> [-2]`, `("pi","pi") -> [0]`, `("s+","pi") -> [1]`, `("pi","s+") -> [1]`. Anchor: `dyad_weights(SIGMA_P, SIGMA_P)[(2, 2)] == approx(1.0)`. Delete the "THE READING IS RAMAN" paragraph. Keep the theta sums unchanged (real linear vectors are reading-independent).
   - Add:

```python
def test_raman_reading_reproduces_the_previous_default():
    reach = lambda e1, e2: sorted({P for (K, P), v in dyad_weights(e1, e2, reading="raman").items() if abs(v) > 1e-12})
    assert reach(SIGMA_P, SIGMA_P) == [0] and reach(SIGMA_P, SIGMA_M) == [2]
    assert dyad_weights(SIGMA_P, SIGMA_M, reading="raman")[(2, 2)] == pytest.approx(-1.0)


def test_K0_is_identity(basis2):
    kets, ctx = basis2
    M0 = two_photon_matrix(kets, kets, ctx, K=0, dOmega=0, P=0)
    assert np.array_equal(M0, np.eye(len(kets)))
```

   - Any other assertion in `tests/test_twophoton.py` or `tests/test_twophoton_closure.py` that encodes the Raman map (grep `sig-,sig+`, `Delta m_F = 0` next to `sigma+, sigma+`, `np.conj(eps2)`) is flipped to ladder or given `reading="raman"` explicitly, whichever the test's own docstring says it is checking.

Check: `conda run -n structure python -m pytest -q tests/test_twophoton.py tests/test_twophoton_closure.py` — all pass, including the two new tests and the negative control.

Commit: `git add heff/twophoton.py tests/test_twophoton.py tests/test_twophoton_closure.py && git commit -m "Two-photon dyad: ladder reading by default, named polarizations, K0 identity" -- heff/twophoton.py tests/test_twophoton.py tests/test_twophoton_closure.py`

### Task 2: `heff/transition.py` and its gates

Files: `heff/transition.py` (new), `tests/test_transition.py` (new), one line in `heff/__init__.py` exporting `transition_matrix, TwoPhotonOperator, DipoleOperator, padded_thf, diagonalize, select`.
Interfaces consumed: Task 1's `dyad_weights(..., reading)`, `leg_weights`, `jones`, `CHANNELS`.
Interfaces produced: everything below, used verbatim by Tasks 3 and 4.

Steps:

1. Write `heff/transition.py`:

```python
"""Transition matrices between eigenstates for any channel-matrix operator.

An operator supplies field-free channel matrices once (``channels``) and a
complex weight per channel for a polarization choice (``weights``). The
sandwich ``sum_k w_k V_a^H M_k V_b`` is done here, amplitudes summed over
channels before squaring, so paths interfere. Design:
docs/plans/2026-09-15-transition-matrix-pipeline-design.md.
"""
from dataclasses import dataclass
from functools import partial

import numpy as np

from .assemble import build_term_matrices, hamiltonian
from .conventions import ef_label, parity_operator
from .elements_c2 import REGISTRY_C2, axial_geometry
from .params import thf_v2
from .spec import Spin, StateSpec, block_by_mF, enumerate_kets, thf_spec
from .spectra import dipole_matrix
from .terms import ctx_from
from .track import assign, pin_col
from .twophoton import CHANNELS, dyad_weights, jones, leg_weights, two_photon_matrix

LABEL_KEYS = ("J", "F1", "F", "ef", "mF")


def padded_thf(iso, *, J_max):
    """Return (kets, ctx, tm, pset) in the two-spin dtype for '232'|'229'|'227'.

    232Th is padded with a spin-0 inner spin. The I_Th = 0 recoupling line of
    ``axial_geometry`` is exactly +1: the 6j with a zero argument collapses to
    (-1)^(J+J'+k) delta(F1,J) delta(F1',J') / sqrt((2J+1)(2J'+1)), which cancels
    that line's own phase and sqrt factor. Gate: test_padding_matches_one_spin.
    """
    pset = thf_v2(iso)
    if iso == "232":
        v1 = thf_spec("232", J_max=J_max)
        spec = StateSpec(case="c", electronic=v1.electronic, I=0.5, J_range=v1.J_range,
                         v=0, M="blocks", frame="rotating",
                         spins=(Spin(label="232Th", I=0.0, couple_to="J"),
                                Spin(label="19F", I=0.5, couple_to="F1")))
    else:
        spec = thf_spec(iso, J_max=J_max)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, pset)
    tm = build_term_matrices(kets, ctx, case="c2", registry=REGISTRY_C2)
    return kets, ctx, tm, pset


@dataclass(frozen=True)
class Eigensystem:
    """Full-basis eigenpairs, block-diagonal in signed m_F, pinned gauge."""
    kets: np.ndarray
    evals: np.ndarray       # (d,)
    evecs: np.ndarray       # (d, d), eigenvector k in column k
    labels: tuple           # per column: dict(J, F1, F, mF, parity, ef, purity)
    field: dict


def _label(kets, v, S, ef_rule, P_op):
    dom = int(np.argmax(np.abs(v)))
    J = float(kets["J"][dom])
    parity = int(round(float(np.real(np.vdot(v, P_op @ v)))))
    return {"J": J, "F1": float(kets["F1"][dom]), "F": float(kets["F"][dom]),
            "mF": float(kets["mF"][dom]), "parity": parity,
            "ef": ef_label(J, parity, rule=ef_rule, S=S, ell=0.0),
            "purity": float(np.abs(v[dom]) ** 2)}


def diagonalize(kets, tm, pset, ctx, *, E_z, B_z, blockwise=True):
    """Eigensystem at one field point. Per signed-m_F block by default.

    Each eigenvector is gauged so its dominant Condon-Shortley basis component
    is real and positive (track.pin_col), the same rule in every block, so
    relative phases between different matrix elements are deterministic.
    ``blockwise=False`` diagonalizes the whole basis at once (gate 5 only).
    """
    H = hamiltonian(tm, pset, {"E_z": E_z, "B_z": B_z})
    d = len(kets)
    evals, evecs = np.zeros(d), np.zeros((d, d))
    groups = block_by_mF(kets).index.values() if blockwise else [np.arange(d)]
    col = 0
    for idx in groups:
        w, v = np.linalg.eigh(H[np.ix_(idx, idx)])
        v = v * pin_col(v.T)[None, :]
        n = len(idx)
        evals[col:col + n] = w
        evecs[np.ix_(idx, np.arange(col, col + n))] = v
        col += n
    order = np.argsort(evals, kind="stable")
    evals, evecs = evals[order], evecs[:, order]
    P_op = parity_operator(kets, ctx.S, ell=0.0, s=0.0)
    labels = tuple(_label(kets, evecs[:, k], ctx.S, pset.conventions.ef_rule, P_op)
                   for k in range(d))
    return Eigensystem(kets, evals, evecs, labels, {"E_z": E_z, "B_z": B_z})


def select(eig, **where):
    """Column indices whose labels match every keyword, sorted by LABEL_KEYS.

    ``select(eig, J=(1, 2), ef="e")``; a scalar or a tuple per key; keys are
    J, F1, F, ef, mF, parity.
    """
    def ok(lab):
        for k, val in where.items():
            vals = val if isinstance(val, (tuple, list, set)) else (val,)
            if lab[k] not in vals:
                return False
        return True
    keep = [k for k, lab in enumerate(eig.labels) if ok(lab)]
    key = lambda k: tuple(eig.labels[k][q] for q in LABEL_KEYS)
    return np.array(sorted(keep, key=key), dtype=int)


class TwoPhotonOperator:
    """Rank-2 effective polarizability, channels keyed (K, |dOmega|, P)."""

    def __init__(self, alphas=None, *, reading="ladder"):
        self.alphas = dict(alphas or {"alpha_K2_dOm0": 1.0, "alpha_K2_dOm2": 1.0})
        self.reading = reading
        bad = set(self.alphas) - {knob for _, _, knob in CHANNELS}
        if bad:
            raise ValueError(f"unknown alphas {sorted(bad)}")

    def _keys(self):
        return [(K, dOm, P) for K, dOm, knob in CHANNELS if knob in self.alphas
                for P in range(-K, K + 1)]

    def channels(self, kets_a, kets_b, ctx):
        return {(K, dOm, P): two_photon_matrix(kets_a, kets_b, ctx, K=K, dOmega=dOm, P=P)
                for K, dOm, P in self._keys()}

    def weights(self, eps1, eps2):
        w = dyad_weights(eps1, eps2, reading=self.reading)
        knob = {(K, dOm): name for K, dOm, name in CHANNELS}
        return {(K, dOm, P): complex(self.alphas[knob[(K, dOm)]]) * w[(K, P)]
                for K, dOm, P in self._keys()}


class DipoleOperator:
    """One-photon E1, channels keyed by p; weights c[p] = (-1)^p eps_{-p}."""

    def channels(self, kets_a, kets_b, ctx):
        geo = partial(axial_geometry, k=1, q=0.0)
        return {p: dipole_matrix(kets_a, kets_b, ctx, p, geometry=geo) for p in (-1, 0, 1)}

    def weights(self, eps):
        return leg_weights(jones(eps))


@dataclass(frozen=True)
class TransitionMatrix:
    amp: np.ndarray         # (n_a, n_b) complex, gauge-fixed
    freqs: np.ndarray       # E_b - E_a, MHz
    labels_a: tuple
    labels_b: tuple
    field: dict
    pol: dict

    @property
    def strength(self):
        return np.abs(self.amp) ** 2


def transition_matrix(op, eig, rows, cols, ctx, *, channels=None, **pol):
    """M_ij = <i|T|j> for eigenvector columns ``rows`` (initial) and ``cols`` (final).

    ``channels`` is the cached output of ``op.channels(eig.kets, eig.kets, ctx)``;
    pass it when calling repeatedly at different fields.
    """
    mats = channels if channels is not None else op.channels(eig.kets, eig.kets, ctx)
    w = op.weights(**pol)
    Va, Vb = eig.evecs[:, rows], eig.evecs[:, cols]
    amp = np.zeros((len(rows), len(cols)), dtype=complex)
    for k, M in mats.items():
        c = complex(w.get(k, 0.0))
        if abs(c) > 0.0:
            amp += c * (Va.conj().T @ M @ Vb)
    freqs = eig.evals[cols][None, :] - eig.evals[rows][:, None]
    return TransitionMatrix(amp, freqs, tuple(eig.labels[k] for k in rows),
                            tuple(eig.labels[k] for k in cols), dict(eig.field), dict(pol))


def sweep_transition_matrix(op, kets, tm, pset, ctx, *, E_z, B_values, rows_at, cols_at, **pol):
    """Tracked M(B) along ``B_values`` at fixed E_z.

    ``rows_at``/``cols_at`` are functions eig -> column indices evaluated at the
    FIRST point; later points are matched to it by overlap within each signed
    m_F block (track.assign), so a curve follows one state through crossings.
    Returns (amps (n_B, n_a, n_b), freqs (n_B, n_a, n_b), labels_a, labels_b).
    """
    channels = op.channels(kets, kets, ctx)
    B_values = np.asarray(B_values, dtype=float)
    eig0 = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=B_values[0])
    rows, cols = rows_at(eig0), cols_at(eig0)
    ref = eig0
    amps, freqs = [], []
    for B in B_values:
        eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=B)
        perm = np.arange(len(kets))
        mref = np.array([l["mF"] for l in ref.labels]); m = np.array([l["mF"] for l in eig.labels])
        for mF in np.unique(mref):
            a, b = np.flatnonzero(mref == mF), np.flatnonzero(m == mF)
            ov = ref.evecs[:, a].T @ eig.evecs[:, b]
            perm[a] = b[assign(ov, mode="overlap")]
        eig = Eigensystem(eig.kets, eig.evals[perm], eig.evecs[:, perm],
                          tuple(eig.labels[k] for k in perm), eig.field)
        t = transition_matrix(op, eig, rows, cols, ctx, channels=channels, **pol)
        amps.append(t.amp); freqs.append(t.freqs)
        ref = eig
    return np.array(amps), np.array(freqs), eig0.labels, eig0.labels
```

   Note for the implementer: `assign(overlap, mode="overlap")` is `heff.track.assign`; read its signature and adjust the call if the strategy keyword is required. `pset.conventions.ef_rule` is how the notebook reads the e/f rule; keep that.

2. Write `tests/test_transition.py` with these gates (all on 232, `padded_thf("232", J_max=3)` unless stated; `SIGMA_P`, `SIGMA_M`, `PI_Z`, `X`, `Y` from `heff.twophoton.POLARIZATIONS`):

```python
import numpy as np
import pytest

from heff.spectra import dipole_matrix  # noqa: F401  (import check)
from heff.transition import (DipoleOperator, TwoPhotonOperator, diagonalize, padded_thf,
                             select, sweep_transition_matrix, transition_matrix)
from heff.twophoton import POLARIZATIONS, two_photon_line_strengths
from heff.wigner import w3j, w6j

SP, SM, PI, X, Y = (POLARIZATIONS[k] for k in ("sigma+", "sigma-", "pi", "x", "y"))
PAIRS = {("sigma+", "sigma+"): 2, ("sigma-", "sigma-"): -2, ("pi", "pi"): 0,
         ("sigma+", "sigma-"): 0, ("sigma+", "pi"): 1, ("sigma-", "pi"): -1}


@pytest.fixture(scope="module")
def sys232():
    kets, ctx, tm, pset = padded_thf("232", J_max=3)
    op = TwoPhotonOperator()
    return kets, ctx, tm, pset, op, op.channels(kets, kets, ctx)


def _eig(s, E_z=0.0, B_z=1.0, **kw):
    kets, ctx, tm, pset, op, ch = s
    return diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=B_z, **kw)


def test_gate1_delta_mF_equals_p1_plus_p2(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2, 3))
    for (e1, e2), dm in PAIRS.items():
        t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1=e1, eps2=e2)
        i, j = np.nonzero(t.strength > 1e-12)
        assert len(i) > 0
        got = {t.labels_b[b]["mF"] - t.labels_a[a]["mF"] for a, b in zip(i, j)}
        assert got == {float(dm)}, (e1, e2, got)


def test_gate2_parity_even_at_zero_E(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=1.0)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2, 3))
    t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="pi", eps2="pi")
    i, j = np.nonzero(t.strength > 1e-12)
    assert all(t.labels_a[a]["ef"] == t.labels_b[b]["ef"] for a, b in zip(i, j))
    eigE = _eig(sys232, E_z=20.0, B_z=1.0)
    tE = transition_matrix(op, eigE, select(eigE, J=1), select(eigE, J=(1, 2, 3)), ctx,
                           channels=ch, eps1="pi", eps2="pi")
    i, j = np.nonzero(tE.strength > 1e-6)
    assert any(tE.labels_a[a]["ef"] != tE.labels_b[b]["ef"] for a, b in zip(i, j))


def test_gate3_sum_rule_in_eigenbasis(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=1e-3)
    rows, cols = select(eig, J=1), np.arange(len(kets))
    for dOm in (0, 2):
        tot = np.zeros(len(rows))
        for P in range(-2, 3):
            M = ch[(2, dOm, P)]
            A = eig.evecs[:, rows].conj().T @ M @ eig.evecs[:, cols]
            tot += np.sum(np.abs(A) ** 2, axis=1)
        assert np.allclose(tot, 1.0, rtol=1e-6), (dOm, tot)


def test_gate4_random_phases_leave_strength_unchanged(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=2)
    t = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="x", eps2="y")
    rng = np.random.default_rng(0)
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi, len(kets)))
    from heff.transition import Eigensystem
    eig2 = Eigensystem(eig.kets, eig.evals, eig.evecs * ph[None, :], eig.labels, eig.field)
    t2 = transition_matrix(op, eig2, rows, cols, ctx, channels=ch, eps1="x", eps2="y")
    assert np.allclose(t.strength, t2.strength)
    assert not np.allclose(t.amp, t2.amp)


def test_gate5_blockwise_equals_full(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    a, b = _eig(sys232, B_z=2.0), _eig(sys232, B_z=2.0, blockwise=False)
    ra, ca = select(a, J=1), select(a, J=(1, 2))
    rb, cb = select(b, J=1), select(b, J=(1, 2))
    ta = transition_matrix(op, a, ra, ca, ctx, channels=ch, eps1="sigma+", eps2="pi")
    tb = transition_matrix(op, b, rb, cb, ctx, channels=ch, eps1="sigma+", eps2="pi")
    assert np.allclose(a.evals, b.evals, atol=1e-8)
    assert np.allclose(ta.strength, tb.strength, atol=1e-10)


def test_gate6_x_is_the_phased_sigma_sum(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    rows, cols = select(eig, J=1), select(eig, J=(1, 2))
    tx = transition_matrix(op, eig, rows, cols, ctx, channels=ch, eps1="x", eps2="x")
    parts = {}
    for e1 in ("sigma+", "sigma-"):
        for e2 in ("sigma+", "sigma-"):
            parts[(e1, e2)] = transition_matrix(op, eig, rows, cols, ctx, channels=ch,
                                                eps1=e1, eps2=e2).amp
    # x = (sigma- - sigma+)/sqrt2 in the Condon-Shortley vectors above
    coh = 0.5 * (parts[("sigma-", "sigma-")] - parts[("sigma-", "sigma+")]
                 - parts[("sigma+", "sigma-")] + parts[("sigma+", "sigma+")])
    assert np.allclose(tx.amp, coh, atol=1e-12)
    incoh = 0.25 * sum(np.abs(v) ** 2 for v in parts.values())
    assert not np.allclose(tx.strength, incoh)


def test_gate7_matches_two_photon_line_strengths(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232, E_z=0.0, B_z=0.0)
    ma, mb = select(eig, mF=0.5), select(eig, mF=2.5)
    t = transition_matrix(op, eig, ma, mb, ctx, channels=ch, eps1="sigma+", eps2="sigma+")
    from heff.assemble import hamiltonian
    from heff.spec import block_by_mF
    blk = block_by_mF(kets)
    ia, ib = blk.index[0.5], blk.index[2.5]
    H = hamiltonian(tm, pset, {"E_z": 0.0, "B_z": 0.0})
    wa, va = np.linalg.eigh(H[np.ix_(ia, ia)])
    wb, vb = np.linalg.eigh(H[np.ix_(ib, ib)])
    _, S = two_photon_line_strengths(wa, va, kets[ia], wb, vb, kets[ib], ctx,
                                     eps1=SP, eps2=SP, alphas=op.alphas)
    assert np.allclose(np.sort(t.strength.ravel()), np.sort(S.ravel()), atol=1e-10)


def test_gate8_padding_matches_one_spin(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    I = 0.5
    def geo1(bra, ket, k, q, p):
        J1, F1, m1, O1 = (float(ket[x]) for x in ("J", "F", "mF", "Om"))
        J2, F2, m2, O2 = (float(bra[x]) for x in ("J", "F", "mF", "Om"))
        if abs(m2 - m1 - p) > 1e-9 or abs(O2 - O1 - q) > 1e-9:
            return 0.0
        return ((-1) ** (F2 - m2) * w3j(F2, k, F1, -m2, p, m1)
                * (-1) ** (F1 + J2 + k + I) * np.sqrt((2 * F2 + 1) * (2 * F1 + 1)) * w6j(J1, F1, I, F2, J2, k)
                * (-1) ** (J2 - O2) * np.sqrt((2 * J2 + 1) * (2 * J1 + 1)) * w3j(J2, k, J1, -O2, q, O1))
    for (K, dOm, P), A in ch.items():
        B = np.zeros_like(A)
        for i in range(len(kets)):
            for j in range(len(kets)):
                q = float(kets["Om"][i]) - float(kets["Om"][j])
                if abs(q) == dOm:
                    B[i, j] = geo1(kets[i], kets[j], K, q, P)
        assert np.max(np.abs(A - B)) < 1e-14, (K, dOm, P)
    from heff.assemble import build_term_matrices, hamiltonian
    from heff.params import thf_v1
    from heff.spec import enumerate_kets, thf_spec
    from heff.terms import REGISTRY, ctx_from
    v1 = thf_spec("232", J_max=3); k1 = enumerate_kets(v1)
    tm1 = build_term_matrices(k1, ctx_from(v1, thf_v1()), case="c", registry=REGISTRY)
    for E, B in ((0.0, 0.0), (10.0, 5.0)):
        H1 = hamiltonian(tm1, thf_v1(), {"E_z": E, "B_z": B})
        H2 = hamiltonian(tm, pset, {"E_z": E, "B_z": B})
        assert np.max(np.abs(H1 - H2)) < 1e-12


def test_sweep_is_continuous(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    B = np.linspace(1e-3, 20.0, 41)
    amps, freqs, la, lb = sweep_transition_matrix(
        op, kets, tm, pset, ctx, E_z=0.0, B_values=B,
        rows_at=lambda e: select(e, J=1), cols_at=lambda e: select(e, J=2),
        eps1="sigma+", eps2="sigma+")
    S = np.abs(amps) ** 2
    assert np.max(np.abs(np.diff(S, axis=0))) < 0.05 * max(S.max(), 1e-12)


def test_dipole_operator_selection_rule(sys232):
    kets, ctx, tm, pset, op, ch = sys232
    eig = _eig(sys232)
    d = DipoleOperator()
    t = transition_matrix(d, eig, select(eig, J=1), select(eig, J=2), ctx, eps="sigma+")
    i, j = np.nonzero(t.strength > 1e-12)
    assert {t.labels_b[b]["mF"] - t.labels_a[a]["mF"] for a, b in zip(i, j)} == {1.0}
```

   Implementer notes: gate 7 uses the module-level `hamiltonian` import in a cleaner form (`from heff.assemble import hamiltonian` at the top); the inline `__import__` above is only to keep the plan self-contained, replace it. In gate 6 the sign convention of "x in terms of sigma±" follows from `POLARIZATIONS`: x = (sigma- − sigma+)/√2 with the vectors as defined; if the assertion fails by a sign, recompute the expansion from the vectors, do not change the operator. Gate 3 tolerance: 1e-6 relative admits the second-order J = 3 admixture at 1e-3 G; if it fails at 1e-6 but passes at 1e-4, report the number and stop.

Check: `conda run -n structure python -m pytest -q tests/test_transition.py tests/test_twophoton.py` — all pass. Expected runtime under 2 min (J_max = 3, 60 kets).

Commit: `git add heff/transition.py heff/__init__.py tests/test_transition.py && git commit -m "transition.py: operator-agnostic transition matrices, per-mF blocks, pinned gauge" -- heff/transition.py heff/__init__.py tests/test_transition.py`

### Task 3: plots and the 232 script

Files: `heff/plot_transition.py` (new), `scripts/plot_thf_twophoton.py` (new), `results/thf-twophoton-2026-09-15/` (written by the script, committed).
Interfaces consumed: Task 2 API verbatim.

Numerical choices (approved by Arian 2026-09-15): 232ThF+, basis J_max = 7. The plotted J = 1..5 to J = 1..5 matrix would be exact at J_max = 5 already (|Delta J| <= 2 keeps every final state inside the basis); J_max = 7 is so that the sum-rule gate stays exact for every initial J up to 5 (each row must see all its J' = J + 2 partners) and so that the Stark admixture of J = 6 into the J = 5 rows is converged. Manifolds J = 1..5 both sides; E_z = 0; fixed-B panels at B = 0.001, 1, 3, 5, 10 G; sweep B = 0.001..20 G, 201 points; alphas = 1.0 placeholders; pairs = the six unordered {sigma+, sigma-, pi} pairs.

Steps:

1. `heff/plot_transition.py`:

```python
"""Heatmaps and B-curves for TransitionMatrix objects. matplotlib only."""
import numpy as np
import matplotlib.pyplot as plt

def _ticks(labels, fine_max=48):
    """Major ticks at J boundaries; minor per-state labels if few enough states."""
    J = np.array([l["J"] for l in labels])
    edges = np.flatnonzero(np.diff(J)) + 0.5
    centers = [np.mean(np.flatnonzero(J == j)) for j in np.unique(J)]
    major = (centers, [f"J={int(j)}" for j in np.unique(J)])
    fine = None
    if len(labels) <= fine_max:
        fine = (np.arange(len(labels)),
                [f"F={l['F']:g} {l['ef']} m={l['mF']:+g}" for l in labels])
    return edges, major, fine


def heatmap(ax, t, *, log=True, floor=1e-6):
    S = t.strength
    Z = np.where(S > floor * S.max(), S, np.nan)
    im = ax.imshow(np.log10(Z) if log else Z, origin="upper", aspect="auto",
                   cmap="viridis")
    for axis, labels in (("x", t.labels_b), ("y", t.labels_a)):
        edges, (pos, names), fine = _ticks(labels)
        for e in edges:
            (ax.axvline if axis == "x" else ax.axhline)(e, color="w", lw=0.5)
        if fine is not None:
            (ax.set_xticks if axis == "x" else ax.set_yticks)(fine[0])
            (ax.set_xticklabels if axis == "x" else ax.set_yticklabels)(fine[1], fontsize=5,
                                                                         rotation=90 if axis == "x" else 0)
        else:
            (ax.set_xticks if axis == "x" else ax.set_yticks)(pos)
            (ax.set_xticklabels if axis == "x" else ax.set_yticklabels)(names)
    ax.set_xlabel("final"); ax.set_ylabel("initial")
    return im


def heatmap_grid(mats, *, title):
    """mats: {(pol1, pol2): TransitionMatrix}. One panel per pair."""
    n = len(mats)
    fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(5 * ((n + 1) // 2), 9), squeeze=False)
    for ax, ((e1, e2), t) in zip(axes.ravel(), mats.items()):
        im = heatmap(ax, t)
        dm = int(round(sum({"sigma+": 1, "sigma-": -1, "pi": 0}[e] for e in (e1, e2))))
        ax.set_title(f"({e1}, {e2})  Delta m_F = {dm:+d}")
        fig.colorbar(im, ax=ax, label="log10 |M|^2 (alpha^2)")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def curves_vs_B(B, amps, labels_a, labels_b, rows, cols, picks, *, title):
    """picks: list of (i, j) index pairs into rows/cols. |M_ij|^2 vs B."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, j in picks:
        la, lb = labels_a[rows[i]], labels_b[cols[j]]
        ax.plot(B, np.abs(amps[:, i, j]) ** 2,
                label=f"J={la['J']:g} F={la['F']:g}{la['ef']} m={la['mF']:+g} -> "
                      f"J={lb['J']:g} F={lb['F']:g}{lb['ef']} m={lb['mF']:+g}")
    ax.set_xlabel("B_z (G)"); ax.set_ylabel("|M|^2 (alpha^2, placeholder)")
    ax.set_title(title); ax.legend(fontsize=6)
    fig.tight_layout()
    return fig
```

2. `scripts/plot_thf_twophoton.py`:

```python
"""Two-photon transition matrices for 232ThF+ X3Delta1, J = 1..5, six polarization pairs.

Run: conda run -n structure python scripts/plot_thf_twophoton.py
"""
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from heff.plot_transition import curves_vs_B, heatmap_grid
from heff.transition import TwoPhotonOperator, diagonalize, padded_thf, select, sweep_transition_matrix, transition_matrix

OUT = ROOT / "results/thf-twophoton-2026-09-15"
ISO, J_MAX, J_MANIFOLD = "232", 7, (1, 2, 3, 4, 5)
E_Z = 0.0
B_PANELS = (0.001, 1.0, 3.0, 5.0, 10.0)
B_SWEEP = np.linspace(0.001, 20.0, 201)
PAIRS = [p for p in itertools.combinations_with_replacement(("sigma+", "sigma-", "pi"), 2)]

OUT.mkdir(parents=True, exist_ok=True)
kets, ctx, tm, pset = padded_thf(ISO, J_max=J_MAX)
op = TwoPhotonOperator()
channels = op.channels(kets, kets, ctx)
np.savez(OUT / "channels.npz", **{f"K{K}_dOm{d}_P{P:+d}": M for (K, d, P), M in channels.items()})

for B in B_PANELS:
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_Z, B_z=B)
    rows, cols = select(eig, J=J_MANIFOLD), select(eig, J=J_MANIFOLD)
    mats = {pair: transition_matrix(op, eig, rows, cols, ctx, channels=channels, eps1=pair[0], eps2=pair[1])
            for pair in PAIRS}
    fig = heatmap_grid(mats, title=f"{ISO}ThF+ X3Delta1 two-photon |M|^2, J=1..5, E_z={E_Z} V/cm, B_z={B} G (alphas placeholder, closure form)")
    fig.savefig(OUT / f"heatmaps_B{B:g}G.png", dpi=200); matplotlib.pyplot.close(fig)
    np.savez(OUT / f"matrices_B{B:g}G.npz", **{f"{a}_{b}": t.amp for (a, b), t in mats.items()},
             evals=eig.evals, rows=rows, cols=cols)
    with open(OUT / f"labels_B{B:g}G.csv", "w") as f:
        f.write("index,J,F1,F,ef,mF,parity,purity,E_MHz\n")
        for k in rows:
            l = eig.labels[k]
            f.write(f"{k},{l['J']:g},{l['F1']:g},{l['F']:g},{l['ef']},{l['mF']:+g},{l['parity']:+d},{l['purity']:.4f},{eig.evals[k]:.6f}\n")

rows_at = lambda e: select(e, J=1)
cols_at = lambda e: select(e, J=(1, 2))
for pair in PAIRS:
    amps, freqs, la, lb = sweep_transition_matrix(op, kets, tm, pset, ctx, E_z=E_Z, B_values=B_SWEEP,
                                                  rows_at=rows_at, cols_at=cols_at, eps1=pair[0], eps2=pair[1])
    eig0 = diagonalize(kets, tm, pset, ctx, E_z=E_Z, B_z=B_SWEEP[0])
    rows, cols = rows_at(eig0), cols_at(eig0)
    S0 = np.abs(amps[0]) ** 2
    picks = [tuple(x) for x in np.argwhere(S0 > 0.2 * S0.max())[:8]]
    fig = curves_vs_B(B_SWEEP, amps, la, lb, rows, cols, picks, title=f"{ISO}ThF+ ({pair[0]}, {pair[1]}) J=1 -> J=1,2 strongest elements vs B_z")
    fig.savefig(OUT / f"curves_{pair[0]}_{pair[1]}.png", dpi=200); matplotlib.pyplot.close(fig)
    np.savez(OUT / f"sweep_{pair[0]}_{pair[1]}.npz", B=B_SWEEP, amps=amps, freqs=freqs, rows=rows, cols=cols)

(OUT / "README.md").write_text(f"""# 232ThF+ two-photon transition matrices, {OUT.name}

Model: heff.transition.padded_thf('232', J_max={J_MAX}); params heff.params.thf_v2('232').
Manifolds: J = {J_MANIFOLD} initial and final (dominant zero-field J). E_z = {E_Z} V/cm.
Fixed-B panels: B_z = {B_PANELS} G. Sweep: {B_SWEEP[0]}..{B_SWEEP[-1]} G, {len(B_SWEEP)} points, J=1 -> J=1,2.
Operator: rank-2 effective polarizability, K=2 only, channels (2,0,P) and (2,2,P), alphas = 1.0 PLACEHOLDER
(no ThF+ value exists). Strengths are geometry in units of alpha^2, not rates. Closure form: valid for
detuning >> ~7 GHz, which is not the JILA regime (docs/open-questions.md OPEN-23). K=1 absent (OPEN-21).
Both photons absorbed: Delta m_F = p1 + p2. Gauge: dominant Condon-Shortley component real positive, per state.
Files: heatmaps_B*.png, matrices_B*.npz (complex amps, keyed pol1_pol2), labels_B*.csv, curves_*.png, sweep_*.npz, channels.npz.
Generated by scripts/plot_thf_twophoton.py.
""")
print("wrote", OUT)
```

Check: `conda run -n structure python scripts/plot_thf_twophoton.py` exits 0, `ls results/thf-twophoton-2026-09-15 | wc -l` gives at least 23 files, and the implementer opens `heatmaps_B1G.png` with the Read tool and confirms: six panels, J-block gridlines on both axes, the (sigma+, sigma+) panel populated only at Delta m_F = +2 (visible as an off-diagonal band). Report the manifold size printed (expected 140 states per side).

Commit: `git add heff/plot_transition.py scripts/plot_thf_twophoton.py results/thf-twophoton-2026-09-15 && git commit -m "232ThF+ two-photon transition matrices: heatmap grid and B sweeps" -- heff/plot_transition.py scripts/plot_thf_twophoton.py results/thf-twophoton-2026-09-15`

### Task 4: docs and notebook sweep to the ladder reading

Files: `docs/thf-plus-x3delta1-effective-hamiltonian.md` (S9.5 only), `docs/open-questions.md`, `docs/architecture.md`, `notebooks/build_isotopologues.py`, `notebooks/ThF_plus_Isotopologues.ipynb` (regenerated).
Interfaces consumed: Task 1 (`reading`, `POLARIZATIONS`).

Steps:

1. `docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5.1(3): the polarization table's row (d) becomes the ladder reading (both absorbed, `(sigma+, sigma+) -> Delta m_F = +2`, `(sigma+, sigma-) -> 0`, `(sigma+, pi) -> +1`); row (e) becomes the Raman reading, marked "available as `reading='raman'`". S9.5.2: add one sentence "The K = 0 channel matrix is exactly the identity on any basis (probe and gate `test_K0_is_identity`, 2026-09-15); transitions are carried by K = 2 alone." S9.5 header: add "Audited 2026-09-15: docs/lit/2026-09-15-twophoton-operator-audit.md". Delete any sentence in S9.5 that says the package follows the Raman reading.

2. `docs/open-questions.md`: OPEN-21 add "The transition-matrix pipeline (`heff/transition.py`) uses K = 2 only." OPEN-23 add "Every figure from `scripts/plot_thf_twophoton.py` carries this caveat in its README."

3. `docs/architecture.md`: after the paragraph on `spectra`/`twophoton`, add one paragraph: "`transition.py` is the operator-agnostic sandwich: an operator supplies field-free channel matrices and per-polarization complex weights; the module diagonalizes per signed m_F block with a pinned gauge, sums amplitudes over channels, then squares. `TwoPhotonOperator` and `DipoleOperator` are the two shipped operators; a new operator is a class with `channels(kets_a, kets_b, ctx)` and `weights(**pol)`."

4. `notebooks/build_isotopologues.py` S11-13: replace the Raman paragraphs with the ladder statement; `TWOP_PAIRS` becomes `{'(sig+,sig+) dmF=+2': (SIGMA_P, SIGMA_P, 2), '(sig+,sig-) dmF=0': (SIGMA_P, SIGMA_M, 0), '(sig-,sig-) dmF=-2': (SIGMA_M, SIGMA_M, -2)}`; the Ng p.102 sentence becomes "`m_F = +3/2 -> -1/2` is `Delta m_F = -2`, reached by a same-helicity `sigma- sigma-` pair"; `_v2_shaped_232` is replaced by `from heff.transition import padded_thf` and its docstring's V16 sentence moves to a one-line comment. Then run:

```
conda run -n structure python notebooks/build_isotopologues.py && conda run -n structure python -m jupyter nbconvert --to notebook --execute --inplace notebooks/ThF_plus_Isotopologues.ipynb
```

Check: the nbconvert command above exits 0 with no `--allow-errors`; `rg -n -i 'raman' docs/thf-plus-x3delta1-effective-hamiltonian.md notebooks/build_isotopologues.py heff/twophoton.py` returns only lines that name the option `reading='raman'` or the Raman row (e) of the table.

Commit: `git add docs/thf-plus-x3delta1-effective-hamiltonian.md docs/open-questions.md docs/architecture.md notebooks/build_isotopologues.py notebooks/ThF_plus_Isotopologues.ipynb && git commit -m "Docs and notebook: ladder reading, K0 identity, transition.py" -- docs/thf-plus-x3delta1-effective-hamiltonian.md docs/open-questions.md docs/architecture.md notebooks/build_isotopologues.py notebooks/ThF_plus_Isotopologues.ipynb`

## Execution log (2026-09-15)

- Task 1: committed e5b4d4d; Check rerun by orchestrator, 17 passed.
- Task 2 corrections to the plan's own code, found by the gates and verified by probe before editing:
  1. `transition_matrix` sandwich was transposed. Channel matrices are indexed `[bra, ket]` (the `(2, 0, +2)` channel has `mF[row] = mF[col] + 2`), so the amplitude for initial `i` -> final `j` is `<j|T|i> = (V_final^H M V_initial)^T`. As written, `(sigma+, sigma+)` gave `Delta m_F = -2` and `DipoleOperator` with `sigma+` gave `-1`. Fixed by the transpose; the Delta m_F = p1 + p2 constraint is unchanged.
  2. Gate 2 compared e/f labels across `Delta J = 1`, where e/f flips at fixed parity. The conserved quantity is parity (68/68 connected pairs at E = 0, B = 1 G). The gate now asserts parity equality (and e/f equality only at equal J), and parity mixing at E = 20 V/cm.
  3. Gate 7 passed the initial block as `kets_a` of `two_photon_line_strengths`, whose `kets_a` is the bra (final) side; with the corrected sandwich that gave an all-zero comparison. The final block is now passed first.
  Task 2 committed 87eb5ba; Check 23 passed (10 gates + 12 twophoton + import).
- Task 4a (docs): committed 20aea77 by the haiku worker; rg checks rerun by orchestrator.
- Task 4b (notebook): committed 2f1bc7c. S13 of the notebook still paired `(sigma+, sigma-)` with `Delta m_F = +2`; under the ladder reading that is `Delta m_F = 0`, so the S13 pair is now `(sigma+, sigma+)`. nbconvert exit 0, zero error outputs.
- Task 3: committed b26d0bb. 140 states per side at J_max = 7; `heatmaps_B1G.png` inspected (six panels, J gridlines); `matrices_B1G.npz` checked numerically: each pair populated only at its `Delta m_F = p1 + p2`.
- Full suite after Task 3: 385 passed, 2 skipped.
- Process note: the sonnet/haiku workers died on the account spend limit mid-task; Tasks 2 (gate 7 fix, `__init__` export), 4b and 3 were finished in place by the orchestrator from the plan's verbatim code after Arian raised the cap.

## Orchestration

The orchestrator is an opus-class agent. It never writes package code itself; it dispatches, verifies by rerunning every task's Check command in its own shell, and decides. Every piece of code a worker needs is in this file, so workers need no physics judgment: sonnet for anything that edits Python, haiku for text edits with an exact before/after and a grep check.

### Split of Task 4 for delegation

Task 4 above is two deliverables with different skill needs; dispatch them as:

- **Task 4a (haiku): doc text edits.** Steps 1-3 of Task 4 only (`docs/thf-plus-x3delta1-effective-hamiltonian.md` S9.5, `docs/open-questions.md`, `docs/architecture.md`). The brief carries the exact sentences to delete and the exact sentences to insert, quoted from Task 4. Check: the `rg -n -i 'raman'` command from Task 4 plus `rg -n 'test_K0_is_identity' docs/` returning one hit. Commit those three files only.
- **Task 4b (sonnet): notebook builder and rebuild.** Step 4 of Task 4. Check: the nbconvert command exits 0. Commit `notebooks/build_isotopologues.py` and `notebooks/ThF_plus_Isotopologues.ipynb` only.

### Agent per task

| Task | Model | Why this tier |
|---|---|---|
| 1 | sonnet | Python edits and test flips, code given verbatim |
| 2 | sonnet | `transition.py` and `test_transition.py` are written out in full above; the work is transcription plus running the gates |
| 3 | sonnet | plotting against Task 2's fixed API, code given |
| 4a | haiku | exact-text replacements with a grep check |
| 4b | sonnet | builder edits then a notebook execution that must not error |

### Order

1. Dispatch Task 1. When it reports, rerun its Check yourself. Do not proceed on the worker's word.
2. Dispatch Tasks 2, 4a and 4b in ONE message (independent: different files, no shared output). 4a and 4b depend on Task 1 only through the word "ladder"; they do not import anything.
3. When Task 2 reports and its Check passes under your own rerun, dispatch Task 3.
4. After Task 3: `conda run -n structure python -m pytest -q` (full suite) and open `results/thf-twophoton-2026-09-15/heatmaps_B1G.png` with the Read tool yourself. Confirm six panels and that the (sigma+, sigma+) panel is populated only off the m_F diagonal (Delta m_F = +2).
5. Report to Arian: commits made, gates passed with the numbers printed, anything a worker could not make pass, the heatmap you looked at.

### Brief template (copy per task, fill the brackets)

```
Task [N] of docs/plans/2026-09-15-transition-matrix-pipeline.md in /Users/arianjadbabaie/Code/heff.
Read that task in full first; the design is docs/plans/2026-09-15-transition-matrix-pipeline-design.md.
Env: conda run -n structure. Tests: conda run -n structure python -m pytest -q <file>.
Write ONLY these files: [list from the task's Files line]. Do not edit any other file.
No git push. No process kills. No edits to heff/elements_c2.py, heff/conventions.py,
heff/models/thf_plus.toml, heff/params.py.
Numerical constants: read them from the plan text, never from memory.
Checkpoint: at the first green run of the task's Check, write
docs/plans/checkpoints/2026-09-15-task[N].md with the command you ran and its last 20 lines.
Update it if you change anything afterwards.
Stop condition: the task's Check passes and the task's Commit command has been run
(commit by the exact pathspec given). Then reply with: the Check command, its output tail,
the commit hash, and any step you could not complete verbatim and why.
If a gate fails after two honest attempts, stop, do not weaken the assertion, and report
the failing values.
```

### Escalation to Arian (orchestrator stops and asks)

- Any gate in Task 2 that fails after the worker's two attempts and your own rerun.
- Gate 3 passing only at a looser tolerance than 1e-6 (report the number).
- Gate 6 failing by an overall sign (the x = (sigma- - sigma+)/sqrt2 expansion): recompute from `POLARIZATIONS`, and if the operator would have to change, stop.
- Any worker that touched a file outside its list (check with `git status` and `git diff --name-only` before accepting a commit).
