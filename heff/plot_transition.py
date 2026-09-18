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


def delta_mF(e1, e2, *, reading="raman"):
    """Delta m_F reached by named polarizations: p1 - p2 (Raman) or p1 + p2 (ladder)."""
    p = {"sigma+": 1, "sigma-": -1, "pi": 0}
    return p[e1] - p[e2] if reading == "raman" else p[e1] + p[e2]


def heatmap_grid(mats, *, title, reading="raman"):
    """mats: {(pol1, pol2): TransitionMatrix}. One panel per pair."""
    n = len(mats)
    fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(5 * ((n + 1) // 2), 9), squeeze=False)
    for ax, ((e1, e2), t) in zip(axes.ravel(), mats.items()):
        im = heatmap(ax, t)
        dm = delta_mF(e1, e2, reading=reading)
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
