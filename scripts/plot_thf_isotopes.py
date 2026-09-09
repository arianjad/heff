"""Reproduce the 2026-09-08 ThF+ field plots; all Hamiltonians come from heff.

Run from the worktree: python scripts/plot_thf_isotopes.py
Uses a separate, fully exported exploratory parameter set; package defaults stay intact.
"""
from pathlib import Path
import sys
import json
import argparse
import csv
from dataclasses import asdict

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

from heff import load_model, parity_operator
from heff.assemble import build_term_matrices, hamiltonian
from heff.terms import ctx_from
from scripts._thf_params import parameters

OUT = ROOT / "results/thf-fields-2026-09-08"
EGRID = np.unique(np.r_[0, np.geomspace(.001, 100, 151), np.linspace(100, 10000, 397),
                        10, 1000, 5000])
BGRID = np.linspace(0, 100, 401)
ISOS = ("232", "229", "227")
TITLES = {"232": "²³²Th¹⁹F⁺ · measured / adopted inputs",
          "229": "²²⁹Th¹⁹F⁺ · quadrupole-omitted baseline",
          "227": "²²⁷Th¹⁹F⁺ · deformed-nucleus theory estimate"}


def matrices(iso, m, p, jmax):
    problem = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=jmax)
    k = problem.kets[problem.kets["mF"] == m]
    # Build from the current source so student edits cannot reuse stale disk matrices.
    tm = build_term_matrices(k, ctx_from(problem.spec, p), case=problem.backend.case,
                             registry=problem.backend.registry, term_names=problem.term_names)
    h0 = hamiltonian(tm, p, {"E_z": 0, "B_z": 0})
    he = hamiltonian(tm, p, {"E_z": 1, "B_z": 0}) - h0
    hb = hamiltonian(tm, p, {"E_z": 0, "B_z": 1}) - h0
    return k, h0, he, hb, tm


def labels(k, v):
    weights = np.array([np.sum(abs(v[k["J"] == j])**2, axis=0) for j in range(1, 9)])
    js = np.argmax(weights, axis=0) + 1
    fsq = (k["F"] * (k["F"] + 1)) @ (abs(v)**2)
    fs = (np.sqrt(1 + 4 * fsq) - 1) / 2
    return js, weights.max(axis=0), fs


def track(h0, op, grid, selected, initial=None):
    """Exact eigenvalues, stepwise character continuation (not dynamical adiabaticity)."""
    w0, v0 = np.linalg.eigh(h0) if initial is None else initial
    energies = [w0]
    prev = v0
    worst = np.ones(len(w0))
    probes = set(grid[::32]) | {0., 10., 100., 1000., 5000., 10000.}
    snapshots = {0.: (w0[selected], v0[:, selected])}
    refinements = 0

    def advance(x0, reference, x, depth=0):
        nonlocal refinements, worst
        w, v = np.linalg.eigh(h0 + x * op)
        overlap = abs(reference.conj().T @ v)**2
        _, idx = linear_sum_assignment(-overlap)
        fidelity = overlap[np.arange(len(idx)), idx]
        if np.min(fidelity[selected]) < .90 and depth < 16:
            refinements += 1
            mid = (x0+x)/2
            _, vm = advance(x0, reference, mid, depth+1)
            return advance(mid, vm, x, depth+1)
        if np.min(fidelity[selected]) < .85:
            raise RuntimeError(f"Unresolved state tracking at E={x:g} V/cm")
        worst = np.minimum(worst, fidelity)
        return w[idx], v[:, idx]

    for x0, x in zip(grid[:-1], grid[1:]):
        w, prev = advance(x0, prev, x)
        energies.append(w)
        if x in probes:
            snapshots[float(x)] = (w[selected], prev[:, selected])
    return np.array(energies), prev, worst, snapshots, refinements


def convergence(k, h0, he, hb, target_vectors, snapshots):
    """Overlap-match embedded J<=7 states to the J<=8 plotted subspace."""
    mask = k["J"] <= 7
    low = np.flatnonzero(mask)
    max_error = 0.
    min_overlap = 1.
    for e, b in [(e, 0) for e in sorted(snapshots)] + [(0, 1), (0, 10), (0, 100)]:
        h = h0 + e * he + b * hb
        n = target_vectors.shape[1]
        if b == 0:
            wt, vt = snapshots[e]
        else:
            w, v = np.linalg.eigh(h)
            _, high_idx = linear_sum_assignment(-abs(target_vectors.conj().T @ v)**2)
            wt, vt = w[high_idx], v[:, high_idx]
        ws, vs = np.linalg.eigh(h[np.ix_(low, low)])
        ov = abs(vt[low].conj().T @ vs)**2
        _, idx = linear_sum_assignment(-ov)
        max_error = max(max_error, float(np.max(abs(wt - ws[idx]))))
        min_overlap = min(min_overlap, float(np.min(ov[np.arange(n), idx])))
    return max_error, min_overlap


def run_case(iso, scenario="baseline", jmax=8):
    problem, p = parameters(iso, scenario)
    dataset = dict(E_V_cm=EGRID, B_G=BGRID)
    rec = {"isotope": iso, "scenario": scenario, "Jmax": jmax,
           "parameters": {n: asdict(q) for n, q in p.params.items()},
           "conventions": p.conventions.stamp(), "blocks": []}
    rows = []
    mmax = {"232": 3.5, "229": 6, "227": 4}[iso]
    mmin = .5 if iso == "232" else 0
    # Only positive m for Stark (negative m is exactly degenerate at B=0).
    # Zeeman uses the time-reversal partner at -B to cover negative m, checked below.
    for m in np.arange(mmin, mmax + .1, 1):
        print(f"{iso} {scenario} m={m:g}: constructing J<={jmax}", flush=True)
        k, h0, he, hb, tm = matrices(iso, m, p, jmax)
        P = parity_operator(k, 1, ell=0, s=0)
        assert np.max(abs(h0 @ P - P @ h0)) < 1e-7
        assert np.max(abs(hb @ P - P @ hb)) < 1e-7
        assert np.max(abs(he @ P + P @ he)) < 1e-7
        w0, v0 = np.linalg.eigh(h0)
        js, weight, fs = labels(k, v0)
        ix = np.flatnonzero(js <= 3)
        assert len(ix) == np.sum(k["J"] <= 3), "Parent subspace is ambiguous"
        es, ve, overlap, snapshots, refinements = track(h0, he, EGRID, ix, (w0, v0))
        err, conv_overlap = convergence(k, h0, he, hb, v0[:, ix], snapshots)
        assert err < .001, f"J cutoff error {err} MHz exceeds the 1 kHz plot target"
        info = dict(mF=float(m), dimension=len(k), parents=len(ix),
                    min_zero_J_weight=float(weight[ix].min()),
                    min_stark_step_overlap=float(overlap[ix].min()),
                    max_J7_to_J8_error_MHz=err, min_cutoff_overlap=conv_overlap,
                    adaptive_tracking_bisections=refinements,
                    convergence_E_V_cm=sorted(snapshots), convergence_B_G=[0, 1, 10, 100])
        rec["blocks"].append(info)
        key = f"m{m:g}"
        dataset[key+"_stark_MHz"] = es[:, ix]
        dataset[key+"_J0"] = js[ix]
        dataset[key+"_F0"] = fs[ix]
        dataset[key+"_E0_MHz"] = w0[ix]
        dataset[key+"_parity0"] = np.real(np.sum(v0[:, ix].conj() * (P @ v0[:, ix]), axis=0))
        # J content at high E is diagnostic, not a conserved quantum number.
        dataset[key+"_J_expect_Emax"] = k["J"] @ abs(ve[:, ix])**2
        for sign in (-1, 1):
            if sign == -1 and m == 0:
                continue
            # In the signed-Omega basis this direct spectral identity is independently checked.
            for parity in (-1, 1):
                pe, pv = np.linalg.eigh(P)
                u = pv[:, np.isclose(pe, parity)]
                hp = u.T @ h0 @ u
                bp = u.T @ hb @ u
                wz, vz = np.linalg.eigh(hp)
                jz, jw, fz = labels(k, u @ vz)
                iz = np.flatnonzero(jz <= 3)
                bz = np.array([np.linalg.eigvalsh(hp + sign*b*bp)[iz] for b in BGRID])
                zkey = f"m{sign*m:g}_p{parity}"
                dataset[zkey+"_zeeman_MHz"] = bz
                dataset[zkey+"_J0"] = jz[iz]
                dataset[zkey+"_E0_MHz"] = wz[iz]
        if m == mmin:
            for i in ix:
                rows.append(dict(J=int(js[i]), F=float(fs[i]), parity=int(round(dataset[key+"_parity0"][list(ix).index(i)])),
                                 energy_MHz=float(w0[i]), J_weight=float(weight[i])))
        print(json.dumps(info), flush=True)
    # Direct negative-m construction verifies the time reversal used above at both requested endpoints.
    test_m = .5 if iso == "232" else 1.
    kp, h, he, hb, _ = matrices(iso, test_m, p, jmax)
    kn, hn, hen, hbn, _ = matrices(iso, -test_m, p, jmax)
    tr_errors = [float(np.max(abs(np.linalg.eigvalsh(hn + e*hen+b*hbn) -
                                 np.linalg.eigvalsh(h + e*he-b*hb))))
                 for e,b in ((0, 100), (10000, 0), (10000, 100))]
    assert max(tr_errors) < 1e-6
    rec["time_reversal_error_MHz"] = max(tr_errors)
    rec["levels"] = rows
    tag = f"{iso}-{scenario}"
    np.savez_compressed(OUT / f"{tag}.npz", **dataset)
    (OUT / f"{tag}.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return dataset, rec


def draw_case(data, rec, combined=None):
    iso = rec["isotope"]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    fig, ax = plt.subplots(3, 3, figsize=(15, 11), layout="constrained")
    mmax = max(b["mF"] for b in rec["blocks"])
    cmap = plt.colormaps["viridis"]
    p = rec["parameters"]
    for row, j in enumerate((1,2,3)):
        rotational = p["B0"]["value"]*j*(j+1) - p["D0"]["value"]*.001*(j*(j+1))**2
        for level in rec["levels"]:
            if level["J"] != j: continue
            x = level["F"] + .10 * level["parity"]
            ax[row,0].plot([x-.07,x+.07], [level["energy_MHz"]-rotational]*2,
                           color="#0072B2" if level["parity"] > 0 else "#D55E00", lw=2)
        for block in rec["blocks"]:
            m = block["mF"]; key=f"m{m:g}"
            pick = data[key+"_J0"] == j
            y = (data[key+"_stark_MHz"] - data[key+"_E0_MHz"])[:,pick]/1000
            ax[row,1].plot(data["E_V_cm"]/1000, y, color=cmap(m/mmax), lw=.8, alpha=.85)
        for key in data:
            if not key.endswith("_zeeman_MHz"): continue
            stem=key.removesuffix("_zeeman_MHz")
            m=float(stem.split("_p")[0][1:]); pick=data[stem+"_J0"] == j
            ax[row,2].plot(data["B_G"], (data[key]-data[stem+"_E0_MHz"])[:,pick],
                           color=cmap(abs(m)/mmax), lw=.8, alpha=.85, ls="--" if m<0 else "-")
        ax[row,0].set(ylabel=f"J₀ = {j}   •   (E − E_rot) / h (MHz)", xlabel="Zero-field F")
        f_ticks=sorted({round(2*level["F"])/2 for level in rec["levels"] if level["J"]==j})
        ax[row,0].set_xticks(f_ticks)
        ax[row,0].set_xlim(min(f_ticks)-.3,max(f_ticks)+.3)
        ax[row,1].set(ylabel="Stark shift ΔE / h (GHz)", xlabel="Electric field (kV/cm)")
        ax[row,2].set(ylabel="Zeeman shift ΔE / h (MHz)", xlabel="Magnetic field (G)")
        for a in ax[row]: a.grid(alpha=.15)
    ax[0,0].set_title("Zero-field hyperfine / parity structure")
    ax[0,1].set_title("Stark: B = 0 · ±mF degenerate")
    ax[0,2].set_title("Zeeman: E = 0 · all signed mF")
    ax[0,0].legend(handles=[Line2D([],[],color="#0072B2",label="p = +1"),
                            Line2D([],[],color="#D55E00",label="p = −1")], fontsize=8)
    cbar=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,mmax),cmap=cmap), ax=ax[:,1:], shrink=.65)
    cbar.set_label("|mF|  (negative-mF Zeeman curves dashed)")
    title = TITLES[iso] if rec["scenario"] == "baseline" else "²²⁹Th¹⁹F⁺ · UNVALIDATED quadrupole sensitivity case"
    fig.suptitle(title+"\nX ³Δ₁ · levels correlated with J = 1–3 · basis J ≤ 8",fontsize=17)
    note="Exploratory model; uncertainty not shown. Zero-field J labels; J mixes in E. Stark: stepwise character tracking; Zeeman: energy order within (mF,p)."
    if iso != "232":
        note += "\nMolecular constants transferred unscaled from ²³²ThF⁺; unknown Th spin rotation omitted; no physical uncertainty band."
    fig.supxlabel(note, fontsize=9)
    tag=f"thf-{iso}-{rec['scenario']}"
    fig.savefig(OUT/f"{tag}.png",dpi=180)
    fig.savefig(OUT/f"{tag}.pdf")
    if combined is not None:
        combined.savefig(fig)
    plt.close(fig)


def draw_overview(records, combined=None):
    fig, axs=plt.subplots(1,3,figsize=(12,7),layout="constrained",sharey=True)
    colors={1:"#0072B2",2:"#D55E00",3:"#009E73"}
    for ax,rec in zip(axs,records):
        lowest=min(x["energy_MHz"] for x in rec["levels"])
        for level in rec["levels"]:
            j=level["J"]; ax.plot([j-.3,j+.3],[(level["energy_MHz"]-lowest)/1000]*2,color=colors[j],alpha=.7)
        ax.set(xticks=[1,2,3],xlabel="Zero-field parent J",title=TITLES[rec["isotope"]].replace(" · ","\n"))
        ax.grid(axis="y",alpha=.2)
    axs[0].set_ylabel("Energy above lowest level / h (GHz)")
    fig.suptitle("ThF⁺ X ³Δ₁ · zero-field rotational structure",fontsize=17)
    fig.supxlabel("Each isotope has its own zero. Odd-isotope molecular constants are unscaled transfers from ²³²ThF⁺.\n"
                  "Unknown Th spin rotation omitted; ²²⁹Th quadrupole omitted. Exploratory models without physical uncertainty bands.",fontsize=10)
    for ext in ("png","pdf"): fig.savefig(OUT/f"thf-level-overview.{ext}",dpi=180)
    if combined is not None:
        combined.savefig(fig)
    plt.close(fig)


def export_tables(cases):
    summaries=[]
    with (OUT/"parameters.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["value","unit","status","uncertainty","source","note"]
        writer=csv.writer(f);writer.writerow(["isotope","scenario","parameter"]+fields)
        for data,rec in cases:
            if rec["scenario"] != "baseline": continue
            for name,q in rec["parameters"].items():
                writer.writerow([rec["isotope"],rec["scenario"],name]+[q.get(k) for k in fields])
            with (OUT/f"{rec['isotope']}-zero-field-levels.csv").open("w",newline="") as levels:
                table=csv.DictWriter(levels,fieldnames=["J","F","parity","energy_MHz","J_weight"])
                table.writeheader();table.writerows(rec["levels"])
            shifts={}
            for field in ("stark","zeeman"):
                values=[]
                suffix=f"_{field}_MHz"
                for key in data:
                    if key.endswith(suffix):
                        values.extend(data[key][-1]-data[key.removesuffix(suffix)+"_E0_MHz"])
                shifts[field]=[min(values),max(values)]
            summaries.append(dict(isotope=rec["isotope"],
                states=sum(b["parents"]*(1 if b["mF"]==0 else 2) for b in rec["blocks"]),
                max_cutoff_error_Hz=1e6*max(b["max_J7_to_J8_error_MHz"] for b in rec["blocks"]),
                min_zero_J_weight=min(b["min_zero_J_weight"] for b in rec["blocks"]),
                min_tracking_overlap=min(b["min_stark_step_overlap"] for b in rec["blocks"]),
                stark_10kV_cm_range_GHz=[x/1000 for x in shifts["stark"]],
                zeeman_100G_range_MHz=shifts["zeeman"],
                time_reversal_error_MHz=rec["time_reversal_error_MHz"]))
    (OUT/"validation.json").write_text(json.dumps(summaries,indent=2),encoding="utf-8")


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--render-only",action="store_true")
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    cases=[]
    for iso,scenario in [(x,"baseline") for x in ISOS]+[("229","legacy-quadrupole")]:
        if args.render_only:
            with np.load(OUT/f"{iso}-{scenario}.npz") as z: data=dict(z)
            rec=json.loads((OUT/f"{iso}-{scenario}.json").read_text())
        else: data,rec=run_case(iso,scenario)
        cases.append((data,rec))
    export_tables(cases)
    with PdfPages(OUT/"thf-isotopes-levels-stark-zeeman.pdf") as combined:
        draw_overview([rec for data,rec in cases if rec["scenario"]=="baseline"],combined)
        for data,rec in cases: draw_case(data,rec,combined)
    print("Finished", OUT, flush=True)


if __name__ == "__main__":
    main()
