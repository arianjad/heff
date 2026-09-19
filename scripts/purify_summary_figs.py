"""Summary figures for the ThF+ purification study (results/thf-purify-2026-09-18*).

Run: conda run -n structure python scripts/purify_summary_figs.py
Reads the per-run README_<variant>.md statistics and histogram PNGs written by scripts/purify_thf.py,
and rebuilds the two-photon graph at E = 0 and 60 V/cm (B = 1 G) for the strength-matrix heatmap.
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from heff.graph import transition_graph
from heff.transition import LABEL_KEYS, TwoPhotonOperator, diagonalize, padded_thf, select, transition_matrix

R = ROOT / "results"
OUT = R / "thf-purify-2026-09-18"
CYC = re.compile(r"Cycles: mean ([\d.]+) \(s\.e\. ([\d.]+)\), median (\d+)")


def stats(readme):
    """[(mean, se, median)] per field-point section, in file order (E=0 first, then E=60 if present)."""
    return [(float(m), float(s), int(md)) for m, s, md in CYC.findall(Path(readme).read_text())]


# --- 1. histogram montage: rows E=0 / E=60, columns = model x library --------------------------------
cols = [("pi model, six open", R / "thf-purify-2026-09-18-budget200" / "cycles_{tag}_six.png"),
        ("pi model, six closed", R / "thf-purify-2026-09-18-budget200" / "cycles_{tag}_six_closed.png"),
        ("propagator, six", R / "thf-purify-2026-09-18-prop" / "cycles_{tag}_six_prop.png"),
        ("propagator, pi (Pipi)", R / "thf-purify-2026-09-18-prop" / "cycles_{tag}_pi_prop.png")]
fig, axes = plt.subplots(2, len(cols), figsize=(4.2 * len(cols), 6.4))
for r, tag in enumerate(("E0_B1", "E60_B1")):
    for c, (name, pat) in enumerate(cols):
        ax = axes[r, c]
        ax.imshow(plt.imread(str(pat).format(tag=tag)))
        ax.set_axis_off()
        rd = pat.parent / f"README_{pat.name.split('cycles_{tag}_')[1][:-4]}.md"
        m, s, md = stats(rd)[r]
        ax.set_title(f"{name}\n{tag.replace('_', ', ')}: mean {m:.1f} +/- {s:.1f}, median {md}", fontsize=9)
fig.suptitle("Cycles to max belief >= 0.99, 100 trajectories, budget 200 (dashed: H(b0) = 7-cycle floor)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "summary_histograms.png", dpi=130)
plt.close(fig)

# --- 2. B sweep: cycles vs B (pi open / pi closed / propagator) + rung anharmonicity ----------------
B = np.array([1, 3.6, 10, 30, 100])
sw = R / "thf-purify-2026-09-18-bsweep"
pr = R / "thf-purify-2026-09-18-prop"
series = {
    "pi model, six open": [stats(sw / f"B{b:g}" / "README_six.md")[0] for b in B],
    "pi model, six closed": [stats(sw / f"B{b:g}" / "README_six_closed.md")[0] for b in B],
}
Bp, prop = [], []
for b in B:
    rd = (pr / "README_six_prop.md") if b == 1 else (pr / f"B{b:g}" / "README_six_prop.md")
    if rd.exists():
        Bp.append(b); prop.append(stats(rd)[0])
# anharmonicity table from the 2026-09-18 zeeman_anharm.py probe (E=0, J<=8; median / population-weighted mean, Hz)
anh_med = np.array([0.178, 2.31, 18.1, 171, 2030])
anh_wtd = np.array([3.58, 46.4, 358, 3220, 35600])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
for name, st in series.items():
    st = np.array(st); a1.errorbar(B, st[:, 0], st[:, 1], marker="o", capsize=3, label=name)
st = np.array(prop); a1.errorbar(Bp, st[:, 0], st[:, 1], marker="s", capsize=3, label="propagator, six")
a1.axhline(6.953, color="r", ls="--", lw=1, label="H(b0) floor")
a1.set_xscale("log"); a1.set_xlabel("B_z (G)"); a1.set_ylabel("mean cycles (+/- s.e.)"); a1.set_title("E_z = 0, six pairs"); a1.legend(fontsize=8); a1.grid(alpha=0.3)
a2.loglog(B, anh_med, "o-", label="median rung")
a2.loglog(B, anh_wtd, "s-", label="population-weighted mean")
a2.loglog(B, anh_med[0] * (B / B[0]) ** 2, "k:", lw=1, label="B^2")
a2.axhline(180, color="g", ls="--", lw=1, label="eta*Omega = 180 Hz"); a2.axhline(1000, color="m", ls="--", lw=1, label="1 kHz window")
a2.set_xlabel("B_z (G)"); a2.set_ylabel("|f(m+1->m+2) - f(m->m+1)| (Hz)"); a2.set_title("Zeeman ladder anharmonicity, J<=8, E_z = 0"); a2.legend(fontsize=8); a2.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(OUT / "summary_bsweep.png", dpi=130)
plt.close(fig)

# --- 3. strength matrix S_ij = |M_ij|^2 (max over the six pairs) for J<=2, and line-count statistics ---
SIX = [("sigma+", "sigma-"), ("sigma-", "sigma+"), ("sigma+", "sigma+"), ("sigma+", "pi"), ("sigma-", "pi"), ("pi", "pi")]
kets, ctx, tm, pset = padded_thf("232", J_max=10)
op = TwoPhotonOperator(); channels = op.channels(kets, kets, ctx)
fig, axes = plt.subplots(2, 2, figsize=(13, 12), gridspec_kw=dict(height_ratios=[3, 1.2]))
for c, E_z in enumerate((0.0, 60.0)):
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=1.0)
    sel = select(eig, J=tuple(range(1, 9)))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=channels, eps1=p[0], eps2=p[1]) for p in SIX}
    G = transition_graph(eig, mats, keep_self=True)
    Smax = G.graph["Smax"]
    idx = select(eig, J=(1, 2))                       # sorted by LABEL_KEYS: J, F1, F, ef, mF
    pos = {int(k): n for n, k in enumerate(idx)}
    S = np.full((len(idx), len(idx)), np.nan)
    for u, v, d in G.edges(data=True):
        if u in pos and v in pos:
            S[pos[u], pos[v]] = np.nanmax([S[pos[u], pos[v]], d["S"]])
    ax = axes[0, c]
    im = ax.imshow(np.log10(S / Smax), vmin=-5, vmax=0, cmap="viridis")
    labels = [f"{eig.labels[k]['J']:g} {eig.labels[k]['F']:g}{eig.labels[k]['ef']} {eig.labels[k]['mF']:+g}" for k in idx]
    ax.set_xticks(range(len(idx))); ax.set_xticklabels(labels, rotation=90, fontsize=5)
    ax.set_yticks(range(len(idx))); ax.set_yticklabels(labels, fontsize=5)
    ax.set_title(f"log10 S_ij / Smax, max over six pairs, J<=2 states (J F ef mF), E_z={E_z:g} V/cm, B_z=1 G", fontsize=9)
    ax.set_xlabel("final j"); ax.set_ylabel("initial i")
    fig.colorbar(im, ax=ax, fraction=0.04)
    # line-strength distribution and out-degree over the whole J<=8 manifold
    ax = axes[1, c]
    Sall = np.array([d["S"] for u, v, d in G.edges(data=True) if u != v])
    deg = np.array([sum(1 for _, _, d in G.out_edges(n, data=True) if d["S"] > 1e-2 * Smax) for n in G.nodes])
    ax.hist(np.log10(Sall / Smax), bins=40, color="C0", alpha=0.7)
    ax.set_xlabel("log10 S / Smax (off-diagonal lines, six pairs, J<=8)"); ax.set_ylabel("lines")
    ax.set_title(f"{G.number_of_edges()} lines over {G.number_of_nodes()} states; out-degree above 1e-2 Smax: "
                 f"median {np.median(deg):.0f}, max {deg.max()}", fontsize=9)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "summary_strength_matrix.png", dpi=130)
plt.close(fig)
# --- 4. every line as a point (initial i, final j, f_ij) over the J<=8 manifold, color = log strength ---
fig = plt.figure(figsize=(16, 7.5))
zt = lambda f: np.sign(f) * np.log10(1 + np.abs(f) / 1e-3)          # symlog, 1 kHz linear scale
ticks = [-5e5, -3e4, -1e3, -1, 0, 1, 1e3, 3e4, 5e5]
for c, E_z in enumerate((0.0, 60.0)):
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=1.0)
    sel = select(eig, J=tuple(range(1, 9)))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=channels, eps1=p[0], eps2=p[1]) for p in SIX}
    G = transition_graph(eig, mats, keep_self=False)
    Smax = G.graph["Smax"]
    pos = {int(k): n for n, k in enumerate(sel)}                      # J-sorted positions
    e = np.array([(pos[u], pos[v], d["f"], d["S"]) for u, v, d in G.edges(data=True)])
    order = np.argsort(e[:, 3])                                       # strong lines drawn last
    ax = fig.add_subplot(1, 2, c + 1, projection="3d")
    sc = ax.scatter(e[order, 0], e[order, 1], zt(e[order, 2]), c=np.log10(e[order, 3] / Smax), cmap="viridis",
                    vmin=-5, vmax=0, s=2 + 10 * (e[order, 3] / Smax) ** 0.5, alpha=0.6, linewidths=0)
    starts = [n for n, k in enumerate(sel) if n == 0 or eig.labels[k]["J"] != eig.labels[sel[n - 1]]["J"]]
    for a in (ax.xaxis, ax.yaxis):
        a.set_ticks(starts); a.set_ticklabels([f"J={eig.labels[sel[n]]['J']:g}" for n in starts], fontsize=7)
    ax.set_zticks([zt(t) for t in ticks]); ax.set_zticklabels([f"{t:g}" for t in ticks], fontsize=7)
    ax.set_xlabel("initial state i (J-sorted)", fontsize=8); ax.set_ylabel("final state j", fontsize=8); ax.set_zlabel("f_ij = E_j - E_i (MHz, symlog)", fontsize=8)
    ax.set_title(f"{len(e)} two-photon lines, six pairs, {len(sel)} states J<=8, E_z={E_z:g} V/cm, B_z=1 G", fontsize=10)
    ax.view_init(elev=22, azim=-55)
fig.colorbar(sc, ax=fig.axes, fraction=0.02, pad=0.02, label="log10 S_ij / Smax")
fig.savefig(OUT / "summary_lines_3d.png", dpi=130)
plt.close(fig)
# --- 5. the Delta J = 0 plane unfolded: f_ij vs initial state i, colored by |Delta m|, E=0 and E=60 ----------
fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
for ax, E_z in zip(axes, (0.0, 60.0)):
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=1.0)
    sel = select(eig, J=tuple(range(1, 9)))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=channels, eps1=p[0], eps2=p[1]) for p in SIX}
    G = transition_graph(eig, mats, keep_self=False)
    pos = {int(k): n for n, k in enumerate(sel)}
    e = np.array([(pos[u], d["f"], abs(G.nodes[u]["mF"] - G.nodes[v]["mF"]), d["S"] / G.graph["Smax"])
                  for u, v, d in G.edges(data=True) if G.nodes[u]["J"] == G.nodes[v]["J"]])
    for dm, col in ((0, "C0"), (1, "C1"), (2, "C3")):
        m = e[:, 2] == dm
        ax.scatter(e[m, 0], zt(e[m, 1]), s=1 + 12 * e[m, 3] ** 0.5, c=col, alpha=0.5, linewidths=0, label=f"|Delta m_F| = {dm}")
    starts = [n for n, k in enumerate(sel) if n == 0 or eig.labels[k]["J"] != eig.labels[sel[n - 1]]["J"]]
    for s in starts[1:]:
        ax.axvline(s - 0.5, color="k", lw=0.5, alpha=0.4)
    ax.set_xticks(starts); ax.set_xticklabels([f"J={eig.labels[sel[n]]['J']:g}" for n in starts])
    zticks = [-1e3, -100, -10, -1, -0.01, 0, 0.01, 1, 10, 100, 1e3]
    ax.set_yticks([zt(t) for t in zticks]); ax.set_yticklabels([f"{t:g}" for t in zticks])
    ax.set_ylabel("f_ij (MHz, symlog)"); ax.grid(alpha=0.3)
    ax.set_title(f"Delta J = 0 lines, {len(e)} of them, E_z={E_z:g} V/cm, B_z=1 G (marker size ~ sqrt S)", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, markerscale=3)
axes[1].set_xlabel("initial state i (J-sorted; within J: F, e/f, m_F)")
fig.tight_layout()
fig.savefig(OUT / "summary_lines_dJ0.png", dpi=130)
plt.close(fig)
# --- 6. the Delta J = +1 and +2 shelves unfolded: f_ij minus the rotational shelf (block-mean E_J' - E_J) ---
# (Delta J = -1, -2 are the same lines with f -> -f.)
fig, axes = plt.subplots(2, 2, figsize=(16, 9), sharex=True)
for c, E_z in enumerate((0.0, 60.0)):
    eig = diagonalize(kets, tm, pset, ctx, E_z=E_z, B_z=1.0)
    sel = select(eig, J=tuple(range(1, 9)))
    mats = {p: transition_matrix(op, eig, sel, sel, ctx, channels=channels, eps1=p[0], eps2=p[1]) for p in SIX}
    G = transition_graph(eig, mats, keep_self=False)
    pos = {int(k): n for n, k in enumerate(sel)}
    EJ = {J: np.mean([eig.evals[k] for k in sel if eig.labels[k]["J"] == J]) for J in range(1, 9)}
    starts = [n for n, k in enumerate(sel) if n == 0 or eig.labels[k]["J"] != eig.labels[sel[n - 1]]["J"]]
    for r, dJ in enumerate((1, 2)):
        ax = axes[r, c]
        e = np.array([(pos[u], d["f"] - (EJ[G.nodes[v]["J"]] - EJ[G.nodes[u]["J"]]), abs(G.nodes[u]["mF"] - G.nodes[v]["mF"]), d["S"] / G.graph["Smax"])
                      for u, v, d in G.edges(data=True) if G.nodes[v]["J"] - G.nodes[u]["J"] == dJ])
        for dm, col in ((0, "C0"), (1, "C1"), (2, "C3")):
            m = e[:, 2] == dm
            ax.scatter(e[m, 0], zt(e[m, 1]), s=1 + 12 * e[m, 3] ** 0.5, c=col, alpha=0.5, linewidths=0, label=f"|Delta m_F| = {dm}")
        for s in starts[1:]:
            ax.axvline(s - 0.5, color="k", lw=0.5, alpha=0.4)
        ax.set_xticks(starts); ax.set_xticklabels([f"J={eig.labels[sel[n]]['J']:g}" for n in starts])
        zticks = [-1e3, -100, -10, -1, -0.01, 0, 0.01, 1, 10, 100, 1e3]
        ax.set_yticks([zt(t) for t in zticks]); ax.set_yticklabels([f"{t:g}" for t in zticks])
        ax.set_ylabel(f"f_ij - (E_J+{dJ} - E_J) (MHz, symlog)"); ax.grid(alpha=0.3)
        ax.set_title(f"Delta J = +{dJ} lines, {len(e)} of them, E_z={E_z:g} V/cm, B_z=1 G (marker size ~ sqrt S)", fontsize=10)
        ax.legend(loc="upper right", fontsize=8, markerscale=3)
    axes[1, c].set_xlabel("initial state i (J-sorted; within J: F, e/f, m_F)")
fig.tight_layout()
fig.savefig(OUT / "summary_lines_dJ12.png", dpi=130)
plt.close(fig)
print("wrote", OUT / "summary_histograms.png", OUT / "summary_bsweep.png", OUT / "summary_strength_matrix.png", OUT / "summary_lines_3d.png", OUT / "summary_lines_dJ0.png", OUT / "summary_lines_dJ12.png")
