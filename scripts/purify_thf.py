"""Greedy sideband-pulse purification of thermal 232ThF+ X3Delta1 from the two-photon graph.

Run: conda run -n structure python scripts/purify_thf.py [--out DIR] [--ntraj N] [--bw MHz] [--pairs six|pi] [--closed]
                                                        [--eta x] [--dark x] [--seed s] [--no-self]
                                                        [--budget N] [--E V/cm ...] [--B G ...]

Deck (see heff.purify for the model), agreed 2026-09-18:
  T = 4 K thermal over the full basis; manifold = J <= J_MANIFOLD (0.999 of population at 4 K).
  Field points (E_z V/cm, B_z G) in FIELDS. PLACEHOLDER alphas.
  Pulse window BW = 0.001 MHz: Pipi et al. 2026 drive sidebands at Omega_ref/2pi = 2 kHz with Lamb-Dicke
  eta = 0.09, so eta*Omega/2pi = 180 Hz, pi pulses of ~3 ms and sub-kHz resolution (per Arian, eta = 0.09).
  --pairs six: all polarization pairs; pi: Pipi's library, E1 = pi with E2 = sigma+/- (Delta m_F = +-1).
  --closed: closed two-level pulses only (drop open Zeeman ladders and diagonal light-shift drives).
  --model prop: per-pulse propagator on addressed states x phonons (detunings, ladders, displacement drives);
    eta*Omega_ref/2pi = --eta-omega (MHz), inclusion/pruning --tol, phonon truncation --nph.
  Readout efficiency ETA, false-click DARK. Policy: greedy expected-information gain; stop at max belief >= TARGET.
"""
import argparse
import sys
import time
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from heff.graph import thermal_populations, transition_graph
from heff.plot_transition import level_panels
from heff.purify import entropy, greedy, info_gain, propagator_library, pulse_library, run
from heff.transition import TwoPhotonOperator, diagonalize, padded_thf, select, transition_matrix

_p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_p.add_argument("--out", default=str(ROOT / "results/thf-purify-2026-09-18"))
_p.add_argument("--ntraj", type=int, default=100)
_p.add_argument("--bw", type=float, default=0.001)
_p.add_argument("--pairs", choices=("six", "pi"), default="six")
_p.add_argument("--closed", action="store_true", help="closed two-level pulses only (pi model)")
_p.add_argument("--model", choices=("pi", "prop"), default="pi", help="pi: top-hat pi-pulse tables; prop: propagator with detunings and ladders")
_p.add_argument("--eta-omega", dest="eta_omega", type=float, default=1.8e-4, help="eta*Omega_ref/2pi in MHz (prop)")
_p.add_argument("--tol", type=float, default=1e-3, help="line inclusion / pruning tolerance (prop)")
_p.add_argument("--nph", type=int, default=8, help="phonon truncation (prop)")
_p.add_argument("--eta", type=float, default=1.0)
_p.add_argument("--dark", type=float, default=0.0)
_p.add_argument("--seed", type=int, default=0)
_p.add_argument("--budget", type=int, default=60, help="max cycles per trajectory")
_p.add_argument("--E", type=float, nargs="+", default=[0.0, 60.0], help="E_z values, V/cm")
_p.add_argument("--B", type=float, nargs="+", default=[1.0], help="B_z values, G")
_p.add_argument("--no-self", dest="self_loops", action="store_false",
                help="drop diagonal (state-dependent light shift) sideband pulses; default keeps them")
A = _p.parse_args()
OUT = Path(A.out); OUT.mkdir(parents=True, exist_ok=True)
ISO, J_MAX, J_MANIFOLD, T, TARGET, MAX_CYCLES = "232", 10, 8, 4.0, 0.99, A.budget
FIELDS = [(E, B) for B in A.B for E in A.E]
SIX = [("sigma+", "sigma-"), ("sigma-", "sigma+"), ("sigma+", "sigma+"),
       ("sigma+", "pi"), ("sigma-", "pi"), ("pi", "pi")]
PAIRS = SIX if A.pairs == "six" else [("sigma+", "pi"), ("sigma-", "pi")]
VARIANT = f"{A.pairs}{'_closed' if A.closed else ''}{'_prop' if A.model == 'prop' else ''}"
assert not (A.closed and A.model == "prop"), "closed applies to the pi model only"

kets, ctx, tm, pset = padded_thf(ISO, J_max=J_MAX)
op = TwoPhotonOperator()
channels = op.channels(kets, kets, ctx)
rng = np.random.default_rng(A.seed)
report = [f"# Greedy purification, {ISO}ThF+ X3Delta1, T={T} K, library `{VARIANT}`\n",
          f"J_max={J_MAX}, manifold J<={J_MANIFOLD}, bw={A.bw} MHz, pairs={A.pairs}, closed={A.closed}, model={A.model}, eta_omega={A.eta_omega} MHz, tol={A.tol}, nph={A.nph}, eta={A.eta}, dark={A.dark}, "
          f"diagonal pulses {A.self_loops}, target max-belief {TARGET}, budget {MAX_CYCLES} cycles, {A.ntraj} trajectories, seed {A.seed}. "
          f"PLACEHOLDER alphas, K=2 only (no scalar K=0 shift).\n"]
for E_z, B_z in FIELDS:
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=B_z)
    sel = select(eig, J=tuple(range(1, J_MANIFOLD + 1)))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=channels, eps1=p[0], eps2=p[1]) for p in PAIRS}
    G = transition_graph(eig, mats, keep_self=A.self_loops)
    states = sorted(G.nodes)
    b0 = thermal_populations(eig, T)[states]; cover = b0.sum(); b0 /= cover
    t0 = time.time()
    lib = (propagator_library(G, states, bw=A.bw, eta_omega=A.eta_omega, tol=A.tol, nph=A.nph) if A.model == "prop"
           else pulse_library(G, states, bw=A.bw, closed=A.closed))
    kmax = max(len(p.u) for p in lib)
    print(f"  library built in {time.time() - t0:.0f} s, largest subspace {kmax} states, tables {sum(p.Mc.nbytes + p.Mn.nbytes for p in lib) / 1e6:.0f} MB", flush=True)
    tag = f"E{E_z:g}_B{B_z:g}"
    print(f"{tag} {VARIANT}: {len(states)} states cover {cover:.4f} of thermal population, H(b0)={entropy(b0):.3f} bits, "
          f"{G.number_of_edges()} edges, {len(lib)} pulses", flush=True)
    first = greedy(b0, lib, eta=A.eta, dark=A.dark)
    print(f"  first pulse: {first.pair} f0={first.f0:.3f} MHz, {len(first.u)} lines, gain {info_gain(b0, first, eta=A.eta, dark=A.dark):.3f} bits", flush=True)
    ncyc, used = [], Counter()
    for _ in range(A.ntraj):
        hist, b, s = run(b0, lib, rng, target=TARGET, max_cycles=MAX_CYCLES, eta=A.eta, dark=A.dark)
        ncyc.append(len(hist))
        used.update((p.pair, round(p.f0, 1)) for p, _ in hist)
    ncyc = np.array(ncyc)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(ncyc, bins=np.arange(ncyc.min() - 0.5, ncyc.max() + 1.5)); ax.axvline(entropy(b0), color="r", ls="--", label="H(b0) bits floor")
    ax.set_xlabel("cycles to max belief >= 0.99"); ax.set_ylabel("trajectories"); ax.legend(); ax.set_title(f"{tag} {VARIANT}")
    fig.tight_layout(); fig.savefig(OUT / f"cycles_{tag}_{VARIANT}.png", dpi=150); plt.close(fig)
    if A.pairs == "six" and not A.closed:
        fig = level_panels(G, J=(1, 2), title=f"{ISO}ThF+ two-photon graph, E_z={E_z:g} V/cm, B_z={B_z:g} G, edges > 1e-2 Smax")
        fig.savefig(OUT / f"graph_J12_{tag}.png", dpi=130); plt.close(fig)
    report.append(f"\n## E_z={E_z:g} V/cm, B_z={B_z:g} G\n\n{len(states)} states, H(b0)={entropy(b0):.3f} bits, {G.number_of_edges()} edges, {len(lib)} candidate pulses. "
                  f"First pulse {first.pair} @ {first.f0:.3f} MHz, {len(first.u)} lines, {info_gain(b0, first, eta=A.eta, dark=A.dark):.3f} bits.\n"
                  f"Cycles: mean {ncyc.mean():.1f} (s.e. {ncyc.std(ddof=1) / np.sqrt(len(ncyc)):.1f}), median {np.median(ncyc):.0f}, "
                  f"min {ncyc.min()}, max {ncyc.max()}, hit budget {int((ncyc >= MAX_CYCLES).sum())}/{A.ntraj}.\n\nMost used pulses (pair, f0 MHz): count\n")
    report += [f"- {k[0][0]},{k[0][1]} @ {k[1]:+.1f}: {n}\n" for k, n in used.most_common(12)]
    print(f"  cycles mean {ncyc.mean():.2f} median {np.median(ncyc):.0f} max {ncyc.max()}", flush=True)
(OUT / f"README_{VARIANT}.md").write_text("".join(report) + "\nGenerated by scripts/purify_thf.py.\n")
print("wrote", OUT, VARIANT)
