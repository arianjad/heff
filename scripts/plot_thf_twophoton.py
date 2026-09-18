"""Two-photon transition matrices for 232ThF+ X3Delta1, J = 1..5, six polarization pairs.

Run: conda run -n structure python scripts/plot_thf_twophoton.py [--case CASE]

--case selects the polarizability anisotropy (all alphas are PLACEHOLDERS):
  default  both rank-2 channels, alpha_K2_dOm0 = alpha_K2_dOm2 = 1, the
           2026-09-15 configuration; writes the 2026-09-15 directory.
  xxyy0    a_xx = a_yy, so alpha_K2_dOm2 ~ (a_xx - a_yy) = 0. The knob is
           OMITTED, not set to zero: TwoPhotonOperator._keys() is keyed on
           presence, so channels.npz then carries only the (2, 0, P) keys.
           A present-but-zero knob gives bit-identical amplitudes (probe
           2026-09-17, max abs difference 0.0) but keeps the dead keys.
  xxyy_pos a_xx - a_yy > 0, alpha_K2_dOm0 = alpha_K2_dOm2 = 1: same numbers
           as default, written to the 2026-09-17 comparison directory.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from heff.plot_transition import curves_vs_B, heatmap_grid
from heff.transition import TwoPhotonOperator, diagonalize, padded_thf, select, sweep_transition_matrix, transition_matrix

#: case -> (alphas, output directory). See the module docstring.
CASES = {
    "default": ({"alpha_K2_dOm0": 1.0, "alpha_K2_dOm2": 1.0}, "results/thf-twophoton-2026-09-15"),
    "xxyy0": ({"alpha_K2_dOm0": 1.0}, "results/thf-twophoton-2026-09-17-anisotropy/xxyy0"),
    "xxyy_pos": ({"alpha_K2_dOm0": 1.0, "alpha_K2_dOm2": 1.0}, "results/thf-twophoton-2026-09-17-anisotropy/xxyy_pos"),
}
_p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_p.add_argument("--case", choices=sorted(CASES), default="default")
CASE = _p.parse_args().case
ALPHAS, _out = CASES[CASE]
OUT = ROOT / _out
ALPHA_TAG = ", ".join(f"{k}={v:g}" for k, v in ALPHAS.items())
ISO, J_MAX, J_MANIFOLD = "232", 7, (1, 2, 3, 4, 5)
E_Z = 0.0
B_PANELS = (0.001, 1.0, 3.0, 5.0, 10.0)
B_SWEEP = np.linspace(0.001, 20.0, 201)
# (eps1 absorbed, eps2 emitted), Raman reading: Delta m_F = p1 - p2 = +2, -2, 0, +1, -1, 0.
# (sigma-, sigma-) equals (sigma+, sigma+) and (pi, sigma+) equals (sigma-, pi) up to sign
# when K = 1 is absent, so six pairs cover the six distinct Delta m_F channels.
PAIRS = [("sigma+", "sigma-"), ("sigma-", "sigma+"), ("sigma+", "sigma+"),
         ("sigma+", "pi"), ("sigma-", "pi"), ("pi", "pi")]

OUT.mkdir(parents=True, exist_ok=True)
kets, ctx, tm, pset = padded_thf(ISO, J_max=J_MAX)
op = TwoPhotonOperator(ALPHAS)
channels = op.channels(kets, kets, ctx)
np.savez_compressed(OUT / "channels.npz", **{f"K{K}_dOm{d}_P{P:+d}": M for (K, d, P), M in channels.items()})

for B in B_PANELS:
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_Z, B_z=B)
    rows, cols = select(eig, J=J_MANIFOLD), select(eig, J=J_MANIFOLD)
    print(f"B_z={B} G: manifold {len(rows)} initial x {len(cols)} final states")
    mats = {pair: transition_matrix(op, eig, rows, cols, ctx, channels=channels, eps1=pair[0], eps2=pair[1])
            for pair in PAIRS}
    fig = heatmap_grid(mats, title=f"{ISO}ThF+ X3Delta1 two-photon |M|^2, J=1..5, E_z={E_Z} V/cm, B_z={B} G\ncase {CASE}: {ALPHA_TAG} (PLACEHOLDER alphas, closure form valid only for detuning >> ~7 GHz)")
    fig.savefig(OUT / f"heatmaps_B{B:g}G.png", dpi=200); matplotlib.pyplot.close(fig)
    np.savez_compressed(OUT / f"matrices_B{B:g}G.npz", **{f"{a}_{b}": t.amp for (a, b), t in mats.items()},
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
    fig = curves_vs_B(B_SWEEP, amps, la, lb, rows, cols, picks, title=f"{ISO}ThF+ ({pair[0]}, {pair[1]}) J=1 -> J=1,2 strongest elements vs B_z\ncase {CASE}: {ALPHA_TAG} (PLACEHOLDER alphas)")
    fig.savefig(OUT / f"curves_{pair[0]}_{pair[1]}.png", dpi=200); matplotlib.pyplot.close(fig)
    np.savez_compressed(OUT / f"sweep_{pair[0]}_{pair[1]}.npz", B=B_SWEEP, amps=amps, freqs=freqs, rows=rows, cols=cols)

(OUT / "README.md").write_text(f"""# 232ThF+ two-photon transition matrices, {OUT.name}

Case --case {CASE}. Model: heff.transition.padded_thf('232', J_max={J_MAX}); params heff.params.thf_v2('232').
Manifolds: J = {J_MANIFOLD} initial and final (dominant zero-field J). E_z = {E_Z} V/cm.
Fixed-B panels: B_z = {B_PANELS} G. Sweep: {B_SWEEP[0]}..{B_SWEEP[-1]} G, {len(B_SWEEP)} points, J=1 -> J=1,2.
Operator: rank-2 effective polarizability, K=2 only, alphas = {ALPHA_TAG}, all PLACEHOLDER
(no ThF+ value exists). Channels (K, |dOmega|) realized: {sorted({(K, d) for K, d, _ in channels})}.
Strengths are geometry in units of alpha^2, not rates. Closure form: valid for
detuning >> ~7 GHz, which is not the JILA regime (docs/open-questions.md OPEN-23). K=1 absent (OPEN-21).
Raman reading: photon 1 (eps1) absorbed, photon 2 (eps2) emitted, Delta m_F = p1 - p2. Gauge: dominant Condon-Shortley component real positive, per state.
Files: heatmaps_B*.png, matrices_B*.npz (complex amps, keyed pol1_pol2), labels_B*.csv, curves_*.png, sweep_*.npz, channels.npz.
Generated by scripts/plot_thf_twophoton.py.
""")
print("wrote", OUT)
