"""heff g^u, g^l (J=1, F=3/2, |M_F|=3/2) vs E, overlaid on Petrov & Skripnikov 2025 Fig. 2.

Petrov curves: digitized from the arXiv:2503.02840 Fig. 2 render (pixel precision ~3e-5 in g),
stored as results/thf-zeeman-vs-petrov-2026-09-14/petrov2025-fig2-digitized.csv.
Petrov Table I: Delta g = Delta g_0 + Delta g_1 E at each tabulated E (Eq. 18).
Run from the repo root: conda run -n structure python scripts/plot_thf_zeeman_vs_petrov.py [out_dir]
Physics: docs/lit/lookup-effective-zeeman-tensor.md.
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from heff.assemble import build_term_matrices
from heff.observe import g_factors
from heff.params import thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import ctx_from

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "thf-zeeman-vs-petrov-2026-09-14"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else RES
J_MAX = 4
MF = 1.5

# Petrov 2025 Table I (E V/cm, 1e7*dg0, 1e7*dg1 cm/V)
TABLE_I = [(40, 349.4, 6.3), (50, 280.2, 7.9), (60, 233.8, 8.7), (70, 200.6, 9.3),
           (80, 175.6, 9.6), (90, 156.1, 9.8), (100, 140.6, 10.0), (110, 127.8, 10.1),
           (120, 117.2, 10.2), (130, 108.2, 10.3), (140, 100.5, 10.3), (150, 93.8, 10.4)]
DG_TEXT_E0 = 2.3e-4   # arXiv:2503.02840 Sec. IV, zero field


def heff_curves(Es):
    spec = thf_spec(J_max=J_MAX)
    kets = enumerate_kets(spec)
    idx = block_by_mF(kets).index[MF]
    pset = thf_v1()
    ctx = ctx_from(spec, pset)
    tm = build_term_matrices(kets[idx], ctx)
    sub = kets[idx]
    gu, gl, Eu, El = [], [], [], []
    for E in Es:
        r = g_factors(tm, pset, {"E_z": float(E), "B_z": 0.0}, ctx=ctx)
        W, V, g = r["W"], r["V"], r["g"]
        # the two J=1, F=3/2 levels: dominant ket J=1, F=3/2, ordered by energy
        sel = [k for k in range(len(W))
               if sub[int(np.argmax(np.abs(V[:, k])))]["J"] == 1.0
               and sub[int(np.argmax(np.abs(V[:, k])))]["F"] == 1.5]
        if len(sel) != 2:
            raise RuntimeError(f"E={E}: found {len(sel)} J=1,F=3/2 levels, expected 2")
        lo, hi = sorted(sel, key=lambda k: W[k])
        gl.append(g[lo]); gu.append(g[hi]); El.append(W[lo]); Eu.append(W[hi])
    return np.array(gu), np.array(gl), np.array(Eu), np.array(El)


def main():
    Es = np.concatenate([np.linspace(0.0, 10.0, 41), np.linspace(10.5, 200.0, 120)])
    gu, gl, Eu, El = heff_curves(Es)
    dg = gu - gl
    pet = np.genfromtxt(RES / "petrov2025-fig2-digitized.csv", delimiter=",", names=True)
    tE = np.array([r[0] for r in TABLE_I]); tdg = np.array([(r[1] + r[2] * r[0]) * 1e-7 for r in TABLE_I])

    with open(OUT / "heff_vs_petrov_fig2.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["E_Vcm", "g_upper", "g_lower", "dg", "E_upper_MHz", "E_lower_MHz"])
        for row in zip(Es, gu, gl, dg, Eu, El):
            w.writerow([f"{x:.10g}" for x in row])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(Es, gu * 1e2, "-", c="C0", lw=1.6, label="heff  g$^u$")
    ax[0].plot(Es, gl * 1e2, "-", c="C3", lw=1.6, label="heff  g$^l$")
    ax[0].plot(pet["E_Vcm"], pet["g_upper"] * 1e2, "o", c="C0", mfc="none", ms=6, label="Petrov Fig. 2  g$^u$ (digitized)")
    ax[0].plot(pet["E_Vcm"], pet["g_lower"] * 1e2, "s", c="C3", mfc="none", ms=6, label="Petrov Fig. 2  g$^l$ (digitized)")
    ax[0].set_xlabel("Lab electric field (V/cm)"); ax[0].set_ylabel("g-factor (10$^{-2}$)")
    ax[0].set_title("J=1, F=3/2, |M$_F$|=3/2"); ax[0].legend(fontsize=8)
    ax[1].plot(Es, dg * 1e4, "-", c="k", lw=1.6, label="heff  $\\Delta$g = g$^u$ − g$^l$")
    ax[1].plot(pet["E_Vcm"], (pet["g_upper"] - pet["g_lower"]) * 1e4, "o", mfc="none", c="C2", label="Petrov Fig. 2 (digitized)")
    ax[1].plot(tE, tdg * 1e4, "^", c="C1", label="Petrov Table I  $\\Delta$g$_0$+$\\Delta$g$_1$E")
    ax[1].plot([0], [DG_TEXT_E0 * 1e4], "*", c="C1", ms=11, label="Petrov text, E=0: 2.3e-4")
    ax[1].axhline(0, c="0.7", lw=0.8)
    ax[1].set_xlabel("Lab electric field (V/cm)"); ax[1].set_ylabel("$\\Delta$g (10$^{-4}$)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "heff_vs_petrov_fig2.png", dpi=150)

    i0 = 0; i60 = int(np.argmin(np.abs(Es - 60.0))); imin = int(np.argmin(np.abs(dg)))
    print(f"E=0     : g_u={gu[i0]:+.6e} g_l={gl[i0]:+.6e} dg={dg[i0]:+.4e}  (Petrov text 2.3e-4; digitized {pet['g_upper'][0]-pet['g_lower'][0]:+.2e})")
    print(f"E=60    : g_u={gu[i60]:+.6e} g_l={gl[i60]:+.6e} dg={dg[i60]:+.4e}  (Table I {tdg[2]:+.4e})")
    print(f"|dg| min: {abs(dg[imin]):.2e} at E={Es[imin]:.1f} V/cm  (Petrov digitized ~6.8e-5 near 24 V/cm)")
    for E, d in zip(tE, tdg):
        k = int(np.argmin(np.abs(Es - E))); print(f"  E={E:5.0f}: heff dg={dg[k]:+.3e}  Table I {d:+.3e}  ratio {dg[k]/d:.3f}")


if __name__ == "__main__":
    main()
