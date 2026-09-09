"""M_F-resolved static level diagrams for ThF+ X3Delta1 isotopologues.

One figure per (isotopologue, field configuration); three axes for J = 1, 2, 3.
Static diagonalisations at fixed field, not the sweeps of plot_thf_isotopes.py.

Run from the worktree: python scripts/plot_thf_mf_resolved.py
"""
from pathlib import Path
from functools import lru_cache
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from heff import load_model, parity_operator
from heff.assemble import build_term_matrices, hamiltonian
from heff.engine import sweep
from heff.terms import ctx_from
from scripts._mf_ticks import format_mF, mf_formatter
from scripts._thf_params import parameters

OUT = ROOT / "results/thf-mf-resolved-2026-09-09"

CONFIGS = ((0.0, 0.0), (0.0, 100.0), (100.0, 0.0), (100.0, 100.0))
JS = (1, 2, 3)
ISOS = ("232", "229", "227")
TITLES = {"232": "232Th19F+ · measured / adopted inputs",
          "229": "229Th19F+ · quadrupole-omitted baseline",
          "227": "227Th19F+ · deformed-nucleus theory estimate"}
PARITY_COLOUR = {1: "#0072B2", -1: "#D55E00"}
HALF = .25                                     # bar half-width, dyadic so centres stay exact


RECORD = np.dtype([("mF", float), ("J", int), ("F", float), ("parity", int),
                   ("index", int), ("E0_MHz", float), ("E_MHz", float)])

# The Omega doublet is Stark-mixed over a few V/cm; a linear ramp resolves the
# adiabatic connection to the zero-field parity states. Energies at the endpoint
# are unchanged at 401 and 801 points, for every block and field configuration.
RAMP = 201


def block(iso, m, p, jmax):
    """Term matrices and zero-field Hamiltonian for one signed-M_F block."""
    problem = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=jmax)
    k = problem.kets[problem.kets["mF"] == m]
    tm = build_term_matrices(k, ctx_from(problem.spec, p), case=problem.backend.case,
                             registry=problem.backend.registry, term_names=problem.term_names)
    return k, tm, hamiltonian(tm, p, {"E_z": 0, "B_z": 0})


def labels(k, v, jmax):
    """Parent J from the largest J weight; F from <F^2>, exact at zero field."""
    weights = np.array([np.sum(abs(v[k["J"] == j])**2, axis=0) for j in range(1, jmax + 1)])
    fsq = (k["F"] * (k["F"] + 1)) @ (abs(v)**2)
    return np.argmax(weights, axis=0) + 1, (np.sqrt(1 + 4 * fsq) - 1) / 2


@lru_cache(maxsize=None)
def zero_field(iso, scenario="baseline", jmax=8):
    """Zero-field parent states of every signed-M_F block, plus their term matrices.

    Parents are the states whose dominant J weight is J <= 3; the M_F range is
    taken from an independently enumerated J_max = 3 basis.
    """
    problem, p = parameters(iso, scenario)
    mmax = float(load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=3).kets["mF"].max())
    blocks = []
    for m in np.arange(-mmax, mmax + .1, 1):
        k, tm, h0 = block(iso, m, p, jmax)
        P = parity_operator(k, 1, ell=0, s=0)
        w0, v0 = np.linalg.eigh(h0)
        js, fs = labels(k, v0, jmax)
        ix = np.flatnonzero(js <= 3)
        parity = np.real(np.sum(v0[:, ix].conj() * (P @ v0[:, ix]), axis=0))
        blocks.append(dict(mF=float(m), kets=k, tm=tm, h0=h0, P=P, ix=ix, J=js[ix],
                           F=np.round(2 * fs[ix]) / 2, parity=parity, E0=w0[ix]))
    return p, tuple(blocks)


@lru_cache(maxsize=None)
def level_records(iso, E_z, B_z, *, scenario="baseline", jmax=8):
    """Per-state records for one isotopologue at one (E_z [V/cm], B_z [G]).

    Field energies follow the zero-field states along a ramped field, so each
    record keeps the zero-field (J, F, parity) label of the state it came from.
    """
    p, blocks = zero_field(iso, scenario, jmax)
    s = np.linspace(0.0, 1.0, RAMP)
    rows = []
    for b in blocks:
        e = sweep(b["tm"], p, {"E_z": E_z * s, "B_z": B_z * s},
                  order="adiabatic_step").evals[-1][b["ix"]]
        for n, i in enumerate(b["ix"]):
            rows.append((b["mF"], b["J"][n], b["F"][n], np.round(b["parity"][n]),
                         i, b["E0"][n], e[n]))
    out = np.array(rows, dtype=RECORD)
    out.setflags(write=False)                      # cached; callers must not mutate
    return out


def centroids(iso, *, scenario="baseline", jmax=8):
    """Per-J zero-field centroid, the constant subtracted in every configuration.

    Taken from the zero-field records, so it does not move with the field and
    panels stay comparable across the four configurations.
    """
    zero = level_records(iso, 0.0, 0.0, scenario=scenario, jmax=jmax)
    return {j: float(zero["E0_MHz"][zero["J"] == j].mean()) for j in JS}


def draw_panel(ax, rec, j, ref):
    """Bars at true M_F, coloured by zero-field parity, F annotated.

    No sideways offset per F: the abscissa is M_F and nothing else, so levels of
    different F at one M_F share an abscissa and separate vertically. HALF is a
    dyadic fraction, which keeps each drawn bar centre exactly on its M_F.
    """
    drawn = rec[rec["J"] == j]
    for r in drawn:
        y = r["E_MHz"] - ref
        ax.plot([r["mF"] - HALF, r["mF"] + HALF], [y, y], lw=1.4,
                color=PARITY_COLOUR[int(r["parity"])], solid_capstyle="butt")
    for f in np.unique(drawn["F"]):
        edge = drawn[drawn["F"] == f]
        outer = edge[edge["mF"] == edge["mF"].max()]
        top = outer[np.argmax(outer["E_MHz"])]
        ax.annotate(f"F = {format_mF(f)}", (top["mF"] + HALF, top["E_MHz"] - ref),
                    xytext=(3, 0), textcoords="offset points", fontsize=7,
                    va="center", color="#444444")
    ticks = np.unique(drawn["mF"])
    ax.set(xticks=ticks, xlim=(ticks.min() - 3 * HALF, ticks.max() + 5 * HALF),
           xlabel="$M_F$", title=f"J = {j}")
    ax.xaxis.set_major_formatter(mf_formatter())
    ax.grid(axis="y", alpha=.15)
    return drawn


def draw_case(iso, E_z, B_z, *, scenario="baseline", jmax=8):
    """One figure per (isotopologue, field configuration); one axis per J."""
    rec = level_records(iso, E_z, B_z, scenario=scenario, jmax=jmax)
    ref = centroids(iso, scenario=scenario, jmax=jmax)
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.5), layout="constrained")
    for ax, j in zip(axes, JS):
        draw_panel(ax, rec, j, ref[j])
    axes[0].set_ylabel("$(E - E_{J,0})/h$  (MHz)")
    axes[0].legend(handles=[Line2D([], [], color=PARITY_COLOUR[1], label="p = +1"),
                            Line2D([], [], color=PARITY_COLOUR[-1], label="p = -1")],
                   fontsize=8, title="zero-field parity", title_fontsize=8)
    fig.suptitle(f"{TITLES[iso]}  ·  $E_z$ = {E_z:g} V/cm, $B_z$ = {B_z:g} G\n"
                 "X $^3\\Delta_1$ · levels correlated with J = 1-3 · basis J <= 8", fontsize=15)
    note = ("Energies referenced to the per-J zero-field centroid, the same constant in all four "
            "configurations. Bars sit at true $M_F$; F levels separate vertically, not sideways.")
    if E_z:
        note += ("\nAt 100 V/cm the molecule is essentially polarized: above the Omega doublet J is a "
                 "correlation label, not a quantum number, and parity is the zero-field label carried "
                 "along the field ramp.")
    if iso != "232":
        note += ("\nMolecular constants transferred unscaled from 232ThF+; unknown Th spin rotation "
                 "omitted; no physical uncertainty band.")
    fig.supxlabel(note, fontsize=8)
    return fig


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for iso in ISOS:
        for E_z, B_z in CONFIGS:
            fig = draw_case(iso, E_z, B_z)
            tag = f"thf-{iso}-E{E_z:g}-B{B_z:g}"
            fig.savefig(OUT / f"{tag}.png", dpi=180)
            fig.savefig(OUT / f"{tag}.pdf")
            plt.close(fig)
            print("wrote", tag, flush=True)
    print("Finished", OUT, flush=True)


if __name__ == "__main__":
    main()
