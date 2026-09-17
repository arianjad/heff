"""Per-B-panel comparison of the two placeholder polarizability-anisotropy cases.

Run after both cases have been generated:
  conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy0
  conda run -n structure python scripts/plot_thf_twophoton.py --case xxyy_pos
  conda run -n structure python scripts/compare_thf_twophoton_anisotropy.py

One figure per B panel, six panels (one per polarization pair). Background is
|M|^2 for xxyy_pos, drawn by heff.plot_transition.heatmap. Overlays mark what the
|dOmega| = 2 pathway does to each element: cells it OPENS (below the floor in
xxyy0), and cells it enhances or suppresses by more than 2x in |M|^2 through
interference with the dOmega = 0 pathway.

TOL is an ABSOLUTE amplitude floor. The 1e-12-relative floor used elsewhere sits
inside float64 dust for these matrices (cancelled elements land at ~1e-13) and
invents hundreds of spurious opened/closed cells; probe 2026-09-17.
"""
import csv
import itertools
import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from heff.plot_transition import heatmap

OUT = ROOT / "results/thf-twophoton-2026-09-17-anisotropy"
PAIRS = list(itertools.combinations_with_replacement(("sigma+", "sigma-", "pi"), 2))
B_PANELS = (0.001, 1.0, 3.0, 5.0, 10.0)
TOL = 1e-10
DM = {"sigma+": 1, "sigma-": -1, "pi": 0}


def read_labels(path):
    """Only the keys heff.plot_transition._ticks reads."""
    with open(path) as f:
        return [{"J": float(r["J"]), "F": float(r["F"]), "ef": r["ef"], "mF": float(r["mF"])}
                for r in csv.DictReader(f)]


for B in B_PANELS:
    tag = f"B{B:g}G"
    d0 = np.load(OUT / f"xxyy0/matrices_{tag}.npz")
    dp = np.load(OUT / f"xxyy_pos/matrices_{tag}.npz")
    labels = read_labels(OUT / f"xxyy_pos/labels_{tag}.csv")
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), squeeze=False)
    print(f"-- {tag}\n{'pair':16s} {'shared':>7s} {'opened':>7s} {'closed':>7s} {'enh>2x':>7s} {'sup<0.5x':>9s}")
    for ax, (e1, e2) in zip(axes.ravel(), PAIRS):
        A, P = d0[f"{e1}_{e2}"], dp[f"{e1}_{e2}"]
        onA, onP = np.abs(A) > TOL, np.abs(P) > TOL
        opened, closed = onP & ~onA, onA & ~onP
        both = onA & onP
        ratio = np.full(A.shape, np.nan)
        ratio[both] = (np.abs(P[both]) ** 2) / (np.abs(A[both]) ** 2)
        enh, sup = ratio > 2.0, ratio < 0.5
        print(f"{e1 + '_' + e2:16s} {int(both.sum()):7d} {int(opened.sum()):7d} "
              f"{int(closed.sum()):7d} {int(enh.sum()):7d} {int(sup.sum()):9d}")
        im = heatmap(ax, SimpleNamespace(strength=np.abs(P) ** 2, labels_a=labels, labels_b=labels))
        for mask, kw, name in (
                (enh, dict(marker="o", s=5, facecolors="none", edgecolors="orangered", linewidths=0.4),
                 f"enhanced >2x ({int(enh.sum())})"),
                (sup, dict(marker="x", s=5, c="deepskyblue", linewidths=0.4),
                 f"suppressed <0.5x ({int(sup.sum())})"),
                (opened, dict(marker="s", s=22, facecolors="none", edgecolors="red", linewidths=0.9),
                 f"opened by |dOmega|=2 ({int(opened.sum())})"),
                (closed, dict(marker="+", s=22, c="magenta", linewidths=0.9),
                 f"closed by |dOmega|=2 ({int(closed.sum())})")):
            ij = np.argwhere(mask)
            ax.scatter(ij[:, 1], ij[:, 0], label=name, **kw)
        ax.set_title(f"({e1}, {e2})  Delta m_F = {DM[e1] + DM[e2]:+d}")
        ax.legend(fontsize=5, loc="upper right", framealpha=0.85)
        fig.colorbar(im, ax=ax, label="log10 |M|^2 (alpha^2), xxyy_pos")
    fig.suptitle(
        f"232ThF+ X3Delta1 two-photon |M|^2, J=1..5, E_z=0 V/cm, B_z={B:g} G: effect of a_xx - a_yy\n"
        "background xxyy_pos (alpha_K2_dOm0=alpha_K2_dOm2=1); overlay vs xxyy0 (alpha_K2_dOm2 omitted, a_xx=a_yy)\n"
        f"PLACEHOLDER alphas, no ThF+ value exists; closure form valid only for detuning >> ~7 GHz (OPEN-23). Floor |M|>{TOL:g}")
    fig.tight_layout()
    fig.savefig(OUT / f"compare_{tag}.png", dpi=200)
    plt.close(fig)
print("wrote", OUT / "compare_B*.png")
