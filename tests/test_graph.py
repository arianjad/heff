import networkx as nx
import numpy as np
import pytest

from heff.graph import thermal_populations, transition_graph
from heff.purify import entropy, greedy, info_gain, outcome_split, pi_pulse, propagator_library, pulse_library, run
from heff.transition import TwoPhotonOperator, diagonalize, padded_thf, select, transition_matrix

SIX = [("sigma+", "sigma-"), ("sigma-", "sigma+"), ("sigma+", "sigma+"),
       ("sigma+", "pi"), ("sigma-", "pi"), ("pi", "pi")]
SIGMA = SIX[:3]


@pytest.fixture(scope="module")
def graph232():
    kets, ctx, tm, pset = padded_thf("232", J_max=3)
    op = TwoPhotonOperator()
    ch = op.channels(kets, kets, ctx)
    eig = diagonalize(kets, tm, pset, ctx, E_z=0.0, B_z=1.0)
    sel = select(eig, J=(1, 2, 3))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=ch, eps1=p[0], eps2=p[1]) for p in SIX}
    return eig, mats, transition_graph(eig, mats)


def test_graph_components_are_parity_and_mF_classes_at_zero_E(graph232):
    eig, mats, G = graph232
    assert G.number_of_nodes() == 60
    assert nx.number_weakly_connected_components(G) == 2          # parity
    Gs = transition_graph(eig, {p: mats[p] for p in SIGMA})
    assert nx.number_weakly_connected_components(Gs) == 4         # parity x (m_F mod 2)
    dm = {"sigma+": 1, "sigma-": -1, "pi": 0}
    for u, v, d in G.edges(data=True):
        assert G.nodes[v]["mF"] - G.nodes[u]["mF"] == dm[d["pair"][0]] - dm[d["pair"][1]]
        assert d["f"] == pytest.approx(G.nodes[v]["E"] - G.nodes[u]["E"])
        assert u != v
    Gd = transition_graph(eig, mats, keep_self=True)
    assert Gd.number_of_edges() > G.number_of_edges()
    assert all(d["f"] == 0.0 and d["pair"][0] == d["pair"][1] for u, v, d in Gd.edges(data=True) if u == v)


def test_thermal_populations_boltzmann(graph232):
    eig, _, _ = graph232
    p = thermal_populations(eig, 4.0)
    assert p.sum() == pytest.approx(1.0)
    i, j = np.argmin(eig.evals), np.argmax(eig.evals)
    assert p[i] / p[j] == pytest.approx(np.exp((eig.evals[j] - eig.evals[i]) / (20836.619 * 4.0)))


def test_outcome_split_conserves_probability_and_selective_pulse_purifies():
    b = np.array([0.5, 0.3, 0.2])
    p = pi_pulse(("pi", "pi"), 0.0, [0, 1], [2, 2], [1.0, 0.5])
    q, bc, bn = outcome_split(b, p)
    assert q == pytest.approx(0.65)
    assert bc.sum() == pytest.approx(1.0) and bn.sum() == pytest.approx(1.0)
    assert q * bc + (1 - q) * bn == pytest.approx(np.array([0.0, 0.15, 0.85]))
    one = pi_pulse(("pi", "pi"), 0.0, [1], [2], [1.0])
    q, bc, bn = outcome_split(b, one)
    assert bc[2] == pytest.approx(1.0) and bn[1] == 0.0
    assert info_gain(b, one) == pytest.approx(entropy(b) - 0.7 * entropy(np.array([5 / 7, 0, 2 / 7])))
    assert 0 < info_gain(b, one, eta=0.9, dark=0.05) < info_gain(b, one)


def test_greedy_purifies_thermal_belief_within_cycle_budget(graph232):
    eig, mats, G = graph232
    states = sorted(G.nodes)
    b0 = thermal_populations(eig, 4.0)[states]; b0 /= b0.sum()
    lib = pulse_library(G, states, bw=0.05)
    assert len(lib) > 10 and all(np.allclose((p.Mc + p.Mn).sum(1), 1) for p in lib)
    Gd = transition_graph(eig, mats, keep_self=True)
    full, closed = pulse_library(Gd, states, bw=0.05), pulse_library(Gd, states, bw=0.05, closed=True)
    assert len(closed) < len(full)                                    # some diagonal or ladder pulses exist
    assert all(not ((p.Mc.sum(1) > 0) & (p.Mc.sum(0) > 0)).any() and not np.diag(p.Mc).any() for p in closed)
    rng = np.random.default_rng(0)
    hist, b, s = run(b0, lib, rng, target=0.99, max_cycles=40)
    assert b.max() >= 0.99 and b.argmax() == s
    assert info_gain(b0, greedy(b0, lib)) > 0.3


def _toy(E, lines):
    """Graph with node energies ``E`` (MHz) and ``lines`` = [(a, b, S)] on one pair."""
    G = nx.MultiDiGraph()
    for s, e in enumerate(E):
        G.add_node(s, E=e)
    for a, b, S in lines:
        G.add_edge(a, b, S=S, f=E[b] - E[a], pair=("pi", "pi"))
    return G


def test_propagator_matches_pi_pulse_on_resonance_and_suppresses_detuned_lines():
    g0 = 1.8e-4                                                     # eta*Omega_ref/2pi, MHz
    lib = propagator_library(_toy([0.0, 10.0], [(0, 1, 1.0)]), [0, 1], bw=0.001, eta_omega=g0)
    (p,) = lib
    assert p.Mc[0, 1] == pytest.approx(1.0, abs=1e-9) and np.allclose((p.Mc + p.Mn).sum(1), 1)
    ref = pi_pulse(("pi", "pi"), 10.0, [0], [1], [1.0])
    assert np.allclose(p.Mc, ref.Mc) and np.allclose(p.Mn, ref.Mn)
    # weaker line on resonance: same sin^2 law as the classical model
    lib = propagator_library(_toy([0.0, 10.0, 0.0, 10.0], [(0, 1, 1.0), (2, 3, 0.25)]), range(4), bw=0.001, eta_omega=g0)
    (p,) = lib
    assert p.Mc[p.u.tolist().index(2), p.u.tolist().index(3)] == pytest.approx(np.sin(np.pi / 4) ** 2, abs=1e-9)
    # a line detuned by 3 g (reference is the resonant line) transfers < g^2/(g^2+D^2) = 0.1
    lib = propagator_library(_toy([0.0, 10.0, 0.0, 10.0 + 3 * g0], [(0, 1, 1.0), (2, 3, 1.0)]), range(4), bw=0.001, eta_omega=g0)
    p = next(q for q in lib if q.f0 == 10.0)
    i, j = p.u.tolist().index(2), p.u.tolist().index(3)
    assert 0.01 < p.Mc[i, j] < 0.1


def test_propagator_ladder_closes_with_anharmonicity_and_displacement_clicks():
    g0 = 1.8e-4
    harmonic = propagator_library(_toy([0.0, 10.0, 20.0], [(0, 1, 1.0), (1, 2, 1.0)]), range(3), bw=0.001, eta_omega=g0)
    p = next(q for q in harmonic if q.f0 == 10.0)
    assert p.Mc[0, 2] > 0.3                                         # |0,0> -> |1,1> -> |2,2> climbs the open ladder
    anharm = propagator_library(_toy([0.0, 10.0, 20.0 + 20 * g0], [(0, 1, 1.0), (1, 2, 1.0)]), range(3), bw=0.001, eta_omega=g0)
    p = next(q for q in anharm if q.f0 == 10.0)
    assert p.Mc[0, 2] < 0.01 and p.Mc[0, 1] > 0.95                  # second rung 20 g off: closed two-level
    # diagonal line = state-dependent force: click probability 1 - exp(-|alpha|^2), alpha = pi/2
    for nph, atol in ((8, 5e-3), (14, 1e-4)):
        (p,) = propagator_library(_toy([0.0], [(0, 0, 1.0)]), [0], bw=0.001, eta_omega=g0, nph=nph)
        assert p.Mc[0, 0] == pytest.approx(1 - np.exp(-(np.pi / 2) ** 2), abs=atol)
