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
    evals, evecs = np.zeros(d), np.zeros((d, d), dtype=H.dtype)
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

    def __init__(self, alphas=None, *, reading="raman"):
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

    def lines(self, *, floor=1e-12):
        """Enumerate the transitions in M_ij: one dict per (i, j) with |M_ij|^2 > floor.

        Each dict carries ``i``, ``j`` (positions in rows/cols), ``freq_MHz``
        (E_final - E_initial), ``strength``, and the initial/final labels as
        ``a_<key>``/``b_<key>`` for J, F1, F, ef, mF, parity. Sorted by frequency.
        """
        S = self.strength
        out = []
        for i, j in zip(*np.nonzero(S > floor)):
            la, lb = self.labels_a[i], self.labels_b[j]
            rec = {"i": int(i), "j": int(j), "freq_MHz": float(self.freqs[i, j]),
                   "strength": float(S[i, j])}
            rec.update({f"a_{k}": la[k] for k in LABEL_KEYS + ("parity",)})
            rec.update({f"b_{k}": lb[k] for k in LABEL_KEYS + ("parity",)})
            out.append(rec)
        return sorted(out, key=lambda r: r["freq_MHz"])


def transition_matrix(op, eig, rows, cols, ctx, *, channels=None, **pol):
    """M_ij = <j|T|i>: row i is the initial eigenvector (``rows``), column j the final one (``cols``).

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
            amp += c * (Vb.conj().T @ M @ Va).T
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
            ov = ref.evecs[:, a].conj().T @ eig.evecs[:, b]
            perm[a] = b[assign(ov, mode="overlap")]
        eig = Eigensystem(eig.kets, eig.evals[perm], eig.evecs[:, perm],
                          tuple(eig.labels[k] for k in perm), eig.field)
        t = transition_matrix(op, eig, rows, cols, ctx, channels=channels, **pol)
        amps.append(t.amp); freqs.append(t.freqs)
        ref = eig
    return np.array(amps), np.array(freqs), eig0.labels, eig0.labels
