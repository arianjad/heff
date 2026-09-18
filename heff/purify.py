"""State purification by sideband pulse, phonon readout and reset, as a classical belief update.

Model, every knob explicit:
- A pulse is one polarization pair, a center frequency ``f0`` and a top-hat
  bandwidth ``bw``; it addresses every graph edge with |f - f0| < bw/2.
- It is a sideband pi pulse calibrated on the strongest addressed line S_ref, so
  an addressed line of strength S transfers with P = sin^2(pi/2 sqrt(S/S_ref))
  and adds one phonon. One destination per initial state (strongest edge kept).
- Readout is binary: click = phonon seen, with efficiency ``eta`` and false-click
  probability ``dark``. Reset removes the phonon; the molecule keeps its new state.
- The thermal state is diagonal and population transfer keeps it diagonal, so the
  conditional state is a probability vector ``b`` over the manifold.
A unitary cannot lower entropy, so per cycle the mean entropy drop is at most the
outcome entropy, <= 1 bit for a binary click. H(b0) bits is the cycle floor.
"""
from collections import defaultdict, namedtuple

import numpy as np

Pulse = namedtuple("Pulse", "pair f0 u v P")   # u, v: positions in ``states``; P: transfer prob per u


def pulse_library(G, states, *, bw, closed=False):
    """Candidate pulses over ``states`` (graph node ids): centers on a bw/2 grid per pair.

    ``closed`` keeps only pulses that are closed two-level systems: no addressed
    state is both a source and a destination, and no diagonal (light-shift) drive.
    An open ladder (|m,0> -> |m+1,1> -> |m+2,2>, or |u,0> -> |u,1> -> |u,2>) stays
    resonant up the phonon ladder and is outside the pi-pulse model above.
    """
    idx = {int(k): n for n, k in enumerate(states)}
    by_pair = defaultdict(list)
    for u, v, d in G.edges(data=True):
        if u in idx and v in idx:
            by_pair[d["pair"]].append((d["f"], idx[u], idx[v], d["S"]))
    lib, seen = [], set()
    for pair, E in by_pair.items():
        E = np.array(E, dtype=float)
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
            key = (pair, tuple(sub[:, 1].astype(int)), tuple(sub[:, 2].astype(int)), tuple(np.round(P, 6)))
            if key in seen:
                continue
            seen.add(key)
            lib.append(Pulse(pair, float(f0), sub[:, 1].astype(int), sub[:, 2].astype(int), P))
    return lib


def entropy(b):
    b = b[b > 0]
    return float(-(b * np.log2(b)).sum())


def outcome_split(b, p, *, eta=1.0, dark=0.0):
    """(q_click, b_click, b_noclick) for belief ``b`` under pulse ``p``."""
    moved = b[p.u] * p.P
    q = float(moved.sum())
    bc = np.zeros_like(b)
    np.add.at(bc, p.v, moved)          # transferred population, now in the final states
    bn = b.copy()
    bn[p.u] -= moved                   # untransferred population
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
        k = np.flatnonzero(p.u == s)
        moved = bool(k.size) and rng.random() < p.P[k[0]]
        if moved:
            s = int(p.v[k[0]])
        click = rng.random() < (eta if moved else dark)
        b = bc if click else bn
        hist.append((p, click))
    return hist, b, s
