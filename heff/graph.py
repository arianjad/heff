"""Eigenstates as nodes, lines as polarization-tagged directed edges (networkx)."""
import networkx as nx
import numpy as np

from .transition import LABEL_KEYS

K_B_MHZ_PER_K = 20836.619   # k_B / h


def thermal_populations(eig, T):
    """Boltzmann weight of every eigenstate at ``T`` kelvin (v = 0, one electronic state)."""
    p = np.exp(-(eig.evals - eig.evals.min()) / (K_B_MHZ_PER_K * T))
    return p / p.sum()


def transition_graph(eig, mats, *, floor=1e-5, keep_self=False):
    """MultiDiGraph over the states touched by ``mats`` = {(eps1, eps2): TransitionMatrix}.

    Node = eigensystem column index, with the label keys, parity, E (MHz).
    Edge u -> v (initial -> final), one per polarization pair, with S = |M|^2,
    f = E_v - E_u (MHz) and pair. ``floor`` is relative to the largest strength
    over all of ``mats``. Self-transitions (the diagonal, a state-dependent
    light shift; on a sideband, Chou et al. Nature 545, 203 (2017)'s readout)
    are dropped unless ``keep_self``. Graph attrs: field, Smax.
    """
    Smax = float(max(t.strength.max() for t in mats.values()))
    G = nx.MultiDiGraph(field=dict(eig.field), Smax=Smax)
    for pair, t in mats.items():
        for k in np.unique(np.concatenate([t.rows, t.cols])):
            if k not in G:
                G.add_node(int(k), E=float(eig.evals[k]),
                           **{q: eig.labels[k][q] for q in LABEL_KEYS + ("parity",)})
        S = t.strength
        for i, j in zip(*np.nonzero(S > floor * Smax)):
            u, v = int(t.rows[i]), int(t.cols[j])
            if keep_self or u != v:
                G.add_edge(u, v, S=float(S[i, j]), f=float(t.freqs[i, j]), pair=pair)
    return G
