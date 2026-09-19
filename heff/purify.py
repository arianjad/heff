"""State purification by sideband pulse, phonon readout and reset, as a classical belief update.

Model, every knob explicit:
- A pulse is one polarization pair and a center frequency ``f0`` on a bw/2 grid.
  It is calibrated as a sideband pi pulse on the strongest line within +-bw/2 of
  f0 (S_ref): that line gets eta*Omega_ref and pulse length t = 1/(2 eta Omega_ref);
  a line of strength S drives at eta*Omega_ref*sqrt(S/S_ref) (power set per pulse).
- Two transfer models, both stored as a ``Pulse`` with joint tables over the
  addressed subspace ``u``: Mc[i, j] (Mn[i, j]) = probability that initial u[i]
  ends in u[j] with (without) a phonon. Rows sum to 1.
  ``pulse_library``: top-hat window |f - f0| < bw/2, every addressed line on
  resonance, one destination per initial state, P = sin^2(pi/2 sqrt(S/S_ref)).
  ``propagator_library``: exact propagator on (addressed states) x (phonons 0..nph)
  with detunings, all lines of the pair that can matter, ladders and diagonal
  (state-dependent force) drives included. See its docstring.
- Readout is binary: click = phonon seen, with efficiency ``eta`` and false-click
  probability ``dark``. Reset removes the phonon; the molecule keeps its new state.
- The thermal state is diagonal and the belief ``b`` is a probability vector over
  the manifold; molecular coherences created within a pulse are dropped at readout.
A unitary cannot lower entropy, so per cycle the mean entropy drop is at most the
outcome entropy, <= 1 bit for a binary click. H(b0) bits is the cycle floor.
"""
from collections import defaultdict, namedtuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

Pulse = namedtuple("Pulse", "pair f0 u Mc Mn")   # u: positions in ``states``; Mc/Mn: k x k joint tables


def _lines(G, states):
    """{pair: array of (f, i, j, S)} with i, j positions in ``states``."""
    idx = {int(k): n for n, k in enumerate(states)}
    by_pair = defaultdict(list)
    for a, b, d in G.edges(data=True):
        if a in idx and b in idx:
            by_pair[d["pair"]].append((d["f"], idx[a], idx[b], d["S"]))
    return {p: np.array(E, dtype=float) for p, E in by_pair.items()}


def _add(lib, seen, pair, f0, u, Mc, Mn):
    key = (pair, tuple(u), Mc.round(6).tobytes(), Mn.round(6).tobytes())
    if key not in seen:
        seen.add(key)
        lib.append(Pulse(pair, float(f0), u, Mc, Mn))


def pi_pulse(pair, f0, u, v, P):
    """Classical sideband pi pulse: initial u[i] -> v[i] with a phonon with probability P[i], else stays."""
    u, v = np.asarray(u, dtype=int), np.asarray(v, dtype=int)
    sub = np.unique(np.concatenate([u, v]))
    pos = {s: n for n, s in enumerate(sub)}
    iu, iv = [pos[s] for s in u], [pos[s] for s in v]
    Mc, Mn = np.zeros((len(sub), len(sub))), np.eye(len(sub))
    Mc[iu, iv] = P
    Mn[iu, iu] = 1 - np.asarray(P, dtype=float)
    return Pulse(pair, float(f0), sub, Mc, Mn)


def pulse_library(G, states, *, bw, closed=False):
    """Classical pi-pulse candidates over ``states`` (graph node ids): centers on a bw/2 grid per pair.

    ``closed`` keeps only pulses that are closed two-level systems: no addressed
    state is both a source and a destination, and no diagonal (light-shift) drive.
    An open ladder (|m,0> -> |m+1,1> -> |m+2,2>, or |u,0> -> |u,1> -> |u,2>) stays
    resonant up the phonon ladder and is outside this model; ``propagator_library``
    treats it.
    """
    lib, seen = [], set()
    for pair, E in _lines(G, states).items():
        for f0 in np.unique(np.round(E[:, 0] / (bw / 2)) * (bw / 2)):
            sub = E[np.abs(E[:, 0] - f0) < bw / 2]
            if not len(sub):
                continue
            sub = sub[np.argsort(-sub[:, 3])]
            _, first = np.unique(sub[:, 1], return_index=True)   # ponytail: strongest edge per initial state
            sub = sub[first]
            u, v = sub[:, 1].astype(int), sub[:, 2].astype(int)
            if closed and (np.all(u == v) or np.isin(v[v != u], u).any()):
                continue
            P = np.sin(np.pi / 2 * np.sqrt(sub[:, 3] / sub[:, 3].max())) ** 2
            p = pi_pulse(pair, f0, u, v, P)
            _add(lib, seen, pair, f0, p.u, p.Mc, p.Mn)
    return lib


def propagator_library(G, states, *, bw, eta_omega=1.8e-4, tol=1e-3, nph=8):
    """Pulses as propagators on (addressed states) x (phonons 0..nph), after Pipi et al. 2026.

    Same centers and calibration as ``pulse_library``: the strongest line within
    +-bw/2 of f0 (S_ref) gets eta*Omega_ref = ``eta_omega`` (MHz; 0.09 x 2 kHz) and
    a pi pulse, t = 1/(2 eta_omega). Every line l of the pair drives at
    g_l = eta_omega sqrt(S_l/S_ref) with detuning D_l = f_l - f0 and enters the
    Hamiltonian if its peak off-resonant transfer g^2/(g^2 + D^2) exceeds ``tol``.
    Rotating frame: diag(s, n) = E_s - n f0; coupling (g_l/2) sqrt(n+1) on
    |a,n> <-> |b,n+1> (blue sideband; carrier and red sideband are detuned by the
    mode frequency, MHz, and dropped). Square pulse. A diagonal line is a
    state-dependent force: coherent displacement up the phonon ladder, so ``nph``
    must hold the displaced state (|alpha|^2 = (pi/2)^2 S/S_ref; nph = 8 keeps
    the n > nph tail below 1e-3). The subspace splits into connected components of
    the included lines, each propagated by dense eigh; states the pulse leaves
    untouched to within ``tol`` are pruned and rows renormalized.
    """
    Es = np.array([G.nodes[s]["E"] for s in states], dtype=float)
    t, m = 1 / (2 * eta_omega), nph + 1                       # us, and phonon levels
    lib, seen = [], set()
    for pair, E in _lines(G, states).items():
        for f0 in np.unique(np.round(E[:, 0] / (bw / 2)) * (bw / 2)):
            win = E[np.abs(E[:, 0] - f0) < bw / 2]
            if not len(win):
                continue
            g = eta_omega * np.sqrt(E[:, 3] / win[:, 3].max())
            keep = g ** 2 / (g ** 2 + (E[:, 0] - f0) ** 2) > tol
            sub, gs = E[keep], g[keep]
            u = np.unique(sub[:, 1:3].astype(int).ravel())
            pos = {s: n for n, s in enumerate(u)}
            k = len(u)
            ia, ib = np.array([pos[int(a)] for a in sub[:, 1]]), np.array([pos[int(b)] for b in sub[:, 2]])
            ncomp, comp = connected_components(coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(k, k)), directed=False)
            Mc, Mn = np.zeros((k, k)), np.zeros((k, k))
            for c in range(ncomp):
                members = np.flatnonzero(comp == c)
                loc = {s: n for n, s in enumerate(members)}
                kc = len(members)
                H = np.zeros((kc * m, kc * m))
                H[np.diag_indices(kc * m)] = ((Es[u[members]] - Es[u[members]].mean())[:, None] - np.arange(m)[None, :] * f0).ravel()
                for li in np.flatnonzero(comp[ia] == c):
                    ra, rb = loc[ia[li]] * m, loc[ib[li]] * m
                    for q in range(nph):
                        cpl = 0.5 * gs[li] * np.sqrt(q + 1)
                        H[rb + q + 1, ra + q] += cpl
                        H[ra + q, rb + q + 1] += cpl
                w, V = np.linalg.eigh(2 * np.pi * H)
                U = (V * np.exp(-1j * w * t)) @ V.conj().T
                P = (np.abs(U[:, ::m]) ** 2).reshape(kc, m, kc)   # [final s', n, initial s]
                blk = np.ix_(members, members)
                Mc[blk], Mn[blk] = P[:, 1:, :].sum(1).T, P[:, 0, :].T
            touched = (1 - np.diag(Mn) > tol) | (Mc.sum(0) + Mn.sum(0) - np.diag(Mn) > tol)
            u, Mc, Mn = u[touched], Mc[np.ix_(touched, touched)], Mn[np.ix_(touched, touched)]
            norm = (Mc + Mn).sum(1, keepdims=True)
            _add(lib, seen, pair, f0, u, Mc / norm, Mn / norm)
    return lib


def entropy(b):
    b = b[b > 0]
    return float(-(b * np.log2(b)).sum())


def outcome_split(b, p, *, eta=1.0, dark=0.0):
    """(q_click, b_click, b_noclick) for belief ``b`` under pulse ``p``."""
    bu = b[p.u]
    mc, mn = bu @ p.Mc, bu @ p.Mn      # population landing in each addressed state with / without a phonon
    q = float(mc.sum())
    bc = np.zeros_like(b)
    bc[p.u] = mc
    bn = b.copy()
    bn[p.u] = mn
    qc = eta * q + dark * (1.0 - q)
    b_click = (eta * bc + dark * bn) / qc if qc > 0 else b
    b_no = ((1 - eta) * bc + (1 - dark) * bn) / (1 - qc) if qc < 1 else b
    return qc, b_click, b_no


def info_gain(b, p, **det):
    q, bc, bn = outcome_split(b, p, **det)
    return entropy(b) - (q * entropy(bc) + (1 - q) * entropy(bn))


def greedy(b, lib, **det):
    return max(lib, key=lambda p: info_gain(b, p, **det))


def run(b0, lib, rng, *, target=0.99, max_cycles=60, eta=1.0, dark=0.0):
    """One trajectory. Returns (history of (pulse, click), final belief, true state)."""
    s = int(rng.choice(len(b0), p=b0))
    b, hist = b0.copy(), []
    while b.max() < target and len(hist) < max_cycles:
        p = greedy(b, lib, eta=eta, dark=dark)
        qc, bc, bn = outcome_split(b, p, eta=eta, dark=dark)
        i = np.flatnonzero(p.u == s)
        phonon = False
        if i.size:
            row = np.concatenate([p.Mn[i[0]], p.Mc[i[0]]])
            j = int(rng.choice(len(row), p=row / row.sum()))
            phonon, s = j >= len(p.u), int(p.u[j % len(p.u)])
        click = rng.random() < (eta if phonon else dark)
        b = bc if click else bn
        hist.append((p, click))
    return hist, b, s
