"""M_F-resolved static level diagrams for ThF+ X3Delta1 isotopologues.

Two figure sets per (isotopologue, field configuration): an overview with one
axis per J, and a zoomed grid with one panel per (J, F1) manifold. Static
diagonalisations at fixed field, not the sweeps of plot_thf_isotopes.py.

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

from heff import MU_B, load_model, parity_operator
from heff.assemble import build_term_matrices, hamiltonian
from heff.engine import sweep
from heff.terms import ctx_from
from scripts._mf_ticks import format_mF, mf_formatter
from scripts._thf_params import parameters

OUT = ROOT / "results/thf-mf-resolved-2026-09-09"
ZOOM = ROOT / "results/thf-mf-resolved-zoom-2026-09-10"

CONFIGS = ((0.0, 0.0), (0.0, 100.0), (100.0, 0.0), (100.0, 100.0))
JS = (1, 2, 3)
ISOS = ("232", "229", "227")
TITLES = {"232": "232Th19F+ · measured / adopted inputs",
          "229": "229Th19F+ · quadrupole-omitted baseline",
          "227": "227Th19F+ · deformed-nucleus theory estimate"}
PARITY_COLOUR = {1: "#0072B2", -1: "#D55E00"}
HALF = .25                                     # bar half-width, dyadic so centres stay exact


RECORD = np.dtype([("mF", float), ("J", int), ("F1", float), ("F", float),
                   ("parity", int), ("index", int), ("E0_MHz", float),
                   ("E_MHz", float)])

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


def _quantum_number(k, v, column):
    """Invert <X^2> = X(X+1) for a basis column, exact at zero field."""
    xsq = (k[column] * (k[column] + 1)) @ (abs(v)**2)
    return (np.sqrt(1 + 4 * xsq) - 1) / 2


def labels(k, v, jmax):
    """Parent J from the largest J weight; F and the intermediate F1 = J + I_Th
    from <F^2> and <F1^2>, exact at zero field. A basis with one nuclear spin
    (232Th19F) carries no F1 column, and F1 falls back to F."""
    weights = np.array([np.sum(abs(v[k["J"] == j])**2, axis=0) for j in range(1, jmax + 1)])
    column = "F1" if "F1" in k.dtype.names else "F"
    return (np.argmax(weights, axis=0) + 1, _quantum_number(k, v, "F"),
            _quantum_number(k, v, column))


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
        js, fs, f1s = labels(k, v0, jmax)
        ix = np.flatnonzero(js <= 3)
        parity = np.real(np.sum(v0[:, ix].conj() * (P @ v0[:, ix]), axis=0))
        blocks.append(dict(mF=float(m), kets=k, tm=tm, h0=h0, P=P, ix=ix, J=js[ix],
                           F1=np.round(2 * f1s[ix]) / 2, F=np.round(2 * fs[ix]) / 2,
                           parity=parity, E0=w0[ix]))
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
            rows.append((b["mF"], b["J"][n], b["F1"][n], b["F"][n],
                         np.round(b["parity"][n]), i, b["E0"][n], e[n]))
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


def manifold_groups(rec):
    """Records split into one group per (J, F1), the panels of the zoomed grid.

    F1 is the intermediate coupling J + I_Th, so a group holds the F values that
    one F1 spawns through the remaining I_F = 1/2: two of them, or one for
    232Th19F, whose F1 falls back to F. The isotopologues differ in how many F1
    a J carries -- 229Th19F reaches six at J = 3 -- not in how many F sit in a
    panel. Keys are exact label values, not tolerance matches.
    """
    keys = sorted({(int(r["J"]), float(r["F1"])) for r in rec})
    return {k: rec[(rec["J"] == k[0]) & (rec["F1"] == k[1])] for k in keys}


def _draw_bars(ax, drawn, ref):
    """Bars at true M_F, coloured by zero-field parity, and the M_F axis.

    No sideways offset per F: the abscissa is M_F and nothing else, so levels of
    different F at one M_F share an abscissa and separate vertically. HALF is a
    dyadic fraction, which keeps each drawn bar centre exactly on its M_F.
    """
    for r in drawn:
        y = r["E_MHz"] - ref
        ax.plot([r["mF"] - HALF, r["mF"] + HALF], [y, y], lw=1.4,
                color=PARITY_COLOUR[int(r["parity"])], solid_capstyle="butt")
    ticks = np.unique(drawn["mF"])
    ax.set(xticks=ticks, xlim=(ticks.min() - 3 * HALF, ticks.max() + 5 * HALF),
           xlabel="$M_F$")
    ax.xaxis.set_major_formatter(mf_formatter())
    ax.grid(axis="y", alpha=.15)


def draw_panel(ax, rec, j, ref):
    """One J of the overview figure: every F of that J on one axis."""
    drawn = rec[rec["J"] == j]
    _draw_bars(ax, drawn, ref)
    for f in np.unique(drawn["F"]):
        edge = drawn[drawn["F"] == f]
        outer = edge[edge["mF"] == edge["mF"].max()]
        top = outer[np.argmax(outer["E_MHz"])]
        ax.annotate(f"F = {format_mF(f)}", (top["mF"] + HALF, top["E_MHz"] - ref),
                    xytext=(3, 0), textcoords="offset points", fontsize=7,
                    va="center", color="#444444")
    ax.set_title(f"J = {j}")
    return drawn


def manifold_reference(iso, *, scenario="baseline", jmax=8):
    """Zero-field centroid of each (J, F1) manifold: the constant its own panel
    is referenced to, which is what keeps a zoomed panel zoomed. The per-J
    centroid of centroids() would drag the window back out to the full
    hyperfine range."""
    zero = level_records(iso, 0.0, 0.0, scenario=scenario, jmax=jmax)
    return {k: float(g["E0_MHz"].mean()) for k, g in manifold_groups(zero).items()}


def _g_factor(band, B_z):
    """g read straight off the drawn bars, or None where B is off.

    dE/dM_F is fitted across both parities of the band, so a parity-dependent
    slope shows up as a bad fit rather than as two numbers to reconcile;
    E = -g mu_B B M_F sets the sign.
    """
    if not B_z or len(np.unique(band["mF"])) < 2:
        return None
    return -float(np.polyfit(band["mF"], band["E_MHz"], 1)[0]) / (MU_B * B_z)


def draw_manifold_grid(iso, E_z, B_z, *, scenario="baseline", jmax=8):
    """Zoomed grid: rows are J, columns the F1 manifolds of that J.

    Every panel in a row shares one energy window, set by the widest manifold in
    that row, so splittings are comparable along a row while each panel stays
    centred on its own manifold. Returns the figure and its {(J, F1): axis} map.
    """
    rec = level_records(iso, E_z, B_z, scenario=scenario, jmax=jmax)
    groups = manifold_groups(rec)
    ref = manifold_reference(iso, scenario=scenario, jmax=jmax)
    per_J = {j: [k for k in groups if k[0] == j] for j in JS}
    ncol = max(len(keys) for keys in per_J.values())
    fig, axes = plt.subplots(len(JS), ncol, squeeze=False, layout="constrained",
                             figsize=(3.7 * ncol + 1, 3.9 * len(JS)))
    panels = {}
    for row, j in enumerate(JS):
        keys = per_J[j]
        half = 1.18 * max(float(np.max(np.abs(groups[k]["E_MHz"] - ref[k]))) for k in keys)
        for col, ax in enumerate(axes[row]):
            if col >= len(keys):
                ax.set_visible(False)                  # ragged: J=1 has fewer F1 than J=3
                continue
            k, panel = keys[col], groups[keys[col]]
            _draw_bars(ax, panel, ref[keys[col]])
            ax.set_ylim(-half, half)
            readout = []
            for f in np.unique(panel["F"]):
                band = panel[panel["F"] == f]
                outer = band[band["mF"] == band["mF"].max()]
                top = outer[np.argmax(outer["E_MHz"])]
                ax.annotate(f"F = {format_mF(f)}",
                            (top["mF"] + HALF, top["E_MHz"] - ref[k]),
                            xytext=(3, 0), textcoords="offset points", fontsize=7,
                            va="center", color="#444444")
                g = _g_factor(band, B_z)
                if g is not None:
                    readout.append(f"F = {format_mF(f)}:  $g_F$ = {g:+.4f}")
            if readout:
                # Boxed inside the axes rather than beside each band: the band
                # labels are short enough to sit at the edge, these are not.
                ax.text(.02, .98, "\n".join(readout), transform=ax.transAxes,
                        va="top", ha="left", fontsize=6.5, color="#444444",
                        bbox=dict(fc="white", alpha=.75, ec="none", pad=1.5))
            # The subtracted centroid, stated because per-panel referencing is
            # exactly what removes it: differences of these recover the spacing
            # between any two manifolds in the figure.
            ax.set_title(f"J = {j},  $F_1$ = {format_mF(k[1])}\n"
                         f"offset {ref[k] / 1e3:.6f} GHz", fontsize=10)
            panels[k] = ax
        axes[row][0].set_ylabel("$(E - E_{J,F_1,0})/h$  (MHz)")
    fig.legend(handles=[Line2D([], [], color=PARITY_COLOUR[1], label="p = +1"),
                        Line2D([], [], color=PARITY_COLOUR[-1], label="p = -1")],
               fontsize=8, title="zero-field parity", title_fontsize=8,
               loc="outside upper right", ncols=2)
    fig.suptitle(f"{TITLES[iso]}  ·  $E_z$ = {E_z:g} V/cm, $B_z$ = {B_z:g} G\n"
                 "X $^3\\Delta_1$ · one panel per $(J, F_1)$ manifold · basis J <= 8",
                 fontsize=15)
    note = ("Each panel is referenced to its own $(J, F_1)$ zero-field centroid, printed as its offset;\n"
            "differences of those offsets give the spacing between any two manifolds.\n"
            "Panels in a row share one energy window, set by the widest manifold in that row.")
    if B_z:
        note += "\n$g_F$ is the fitted d$E$/d$M_F$ of that F band divided by $-\\mu_B B_z$."
    fig.supxlabel(note, fontsize=8)
    return fig, panels


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


def main_zoom():
    """The zoomed companion set: one panel per (J, F1) manifold."""
    ZOOM.mkdir(parents=True, exist_ok=True)
    for iso in ISOS:
        for E_z, B_z in CONFIGS:
            fig, _ = draw_manifold_grid(iso, E_z, B_z)
            tag = f"thf-{iso}-zoom-E{E_z:g}-B{B_z:g}"
            fig.savefig(ZOOM / f"{tag}.png", dpi=180)
            fig.savefig(ZOOM / f"{tag}.pdf")
            plt.close(fig)
            print("wrote", tag, flush=True)
    print("Finished", ZOOM, flush=True)


if __name__ == "__main__":
    main()
    main_zoom()
