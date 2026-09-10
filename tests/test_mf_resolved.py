"""Physical invariants of the M_F-resolved ThF+ level records and figures.

The records come from `scripts/plot_thf_mf_resolved.py`; the checks compare them
against an independently enumerated basis, against exact symmetries of the
Hamiltonian, and against `heff.observe.g_factors`.
"""
import re
from collections import Counter

import numpy as np

from heff import g_factors, load_model
from heff.assemble import hamiltonian
from heff.terms import ctx_from
from scripts._thf_params import parameters
from scripts.plot_thf_mf_resolved import (CONFIGS, JS, centroids, draw_panel,
                                          draw_manifold_grid, level_records,
                                          manifold_groups, manifold_reference,
                                          plt, zero_field)

ISOS = ("232", "229", "227")


def test_parent_J_selection_reproduces_the_full_Jmax3_state_counts():
    """Selecting parent J<=3 out of the J_max=8 basis must reproduce, per (J, F),
    the state count of an independently enumerated J_max=3 basis: one state per
    (J, F, M_F, Omega) with Omega = +/-1. A dropped, duplicated, or mislabelled
    signed-M_F block changes one of these counts."""
    for iso in ISOS:
        kets = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=3).kets
        want = Counter(zip(kets["J"].tolist(), kets["F"].tolist()))
        rec = level_records(iso, 0.0, 0.0)
        got = Counter(zip(rec["J"].tolist(), rec["F"].tolist()))
        assert got == want, iso


def test_zero_field_energies_are_M_F_degenerate_within_a_J_F_parity_group():
    """At E = B = 0 nothing breaks rotational invariance, so a (J, F, parity)
    group's zero-field energies cannot depend on M_F. Compared as a sorted
    spectrum per M_F because (J, F, parity) is itself doubly degenerate for the
    two-nuclear-spin isotopologues (227Th J=1 carries two F=1 states).

    Tolerance 1e-6 MHz = 1 Hz: five orders below the smallest physical structure
    in these manifolds (the ~5 MHz Omega doublet), and about three orders above
    the eigenvalue round-off bound n*eps*||H|| ~ 4e-9 MHz at ||H|| ~ 1e5 MHz.
    """
    for iso in ISOS:
        groups = {}
        for r in level_records(iso, 0.0, 0.0):
            key = (r["J"], r["F"], r["parity"])
            groups.setdefault(key, {}).setdefault(r["mF"], []).append(r["E0_MHz"])
        for key, per_mF in groups.items():
            spectra = [np.sort(v) for v in per_mF.values()]
            assert len({len(s) for s in spectra}) == 1, (iso, key)
            spread = float(np.max(np.abs(np.array(spectra) - spectra[0])))
            assert spread < 1e-6, (iso, key, spread)


def test_stark_energies_are_degenerate_between_plus_and_minus_M_F():
    """At B = 0 the Hamiltonian is invariant under reflection in a plane
    containing z, so E(+M_F) = E(-M_F) at any E_z. Compared per
    (J, F, parity, |M_F|) as sorted spectra.

    Tolerance 1e-7 MHz: an exact symmetry, so the only budget is round-off. The
    eigenvalue bound n*eps*||H|| is ~4e-9 MHz for the largest block (n = 176,
    ||H|| ~ 1e5 MHz); 1e-7 MHz allows a factor ~25 on that and still sits nine
    orders below the ~75 MHz Stark structure being plotted.
    """
    for iso in ISOS:
        groups = {}
        for r in level_records(iso, 100.0, 0.0):
            key = (r["J"], r["F"], r["parity"], abs(r["mF"]))
            groups.setdefault(key, {}).setdefault(np.sign(r["mF"]), []).append(r["E_MHz"])
        for key, per_sign in groups.items():
            if len(per_sign) < 2:
                continue                            # M_F = 0 has no partner
            spectra = [np.sort(v) for v in per_sign.values()]
            assert len({len(s) for s in spectra}) == 1, (iso, key)
            gap = float(np.max(np.abs(np.array(spectra) - spectra[0])))
            assert gap < 1e-7, (iso, key, gap)


def test_time_reversal_maps_plus_M_F_at_plus_B_onto_minus_M_F_at_minus_B():
    """With E_z held fixed, time reversal sends (M_F, B_z) -> (-M_F, -B_z) and
    leaves the energy unchanged. The combined E and B configuration is where a
    sign error in the signed-M_F block or in the Zeeman knob would show, since
    neither the Stark nor the Zeeman symmetry alone constrains it.

    Tolerance 1e-7 MHz: an exact symmetry, budget set by round-off as in the
    Stark check, and nine orders below the plotted structure.
    """
    for iso in ISOS:
        up, down = {}, {}
        for r in level_records(iso, 100.0, 100.0):
            up.setdefault((r["J"], r["F"], r["parity"], r["mF"]), []).append(r["E_MHz"])
        for r in level_records(iso, 100.0, -100.0):
            down.setdefault((r["J"], r["F"], r["parity"], -r["mF"]), []).append(r["E_MHz"])
        assert set(up) == set(down), iso
        for key in up:
            err = float(np.max(np.abs(np.sort(up[key]) - np.sort(down[key]))))
            assert err < 1e-7, (iso, key, err)


def test_zeeman_shift_matches_the_g_factors_slope_computed_independently():
    """The plotted Zeeman shift must reproduce E = -g mu_B B M_F with g taken
    from heff.observe.g_factors, which gets the slope by exact Hellmann-Feynman
    and sum-over-states derivatives at B = 0 rather than by diagonalising at
    100 G. This compares the figure's numbers against a separate package API
    instead of against themselves.

    The residual is genuinely second order in B, so the tolerance is the exact
    second derivative g_factors also returns: |shift - slope*B| <= |d2| B^2
    (twice the second-order term 0.5*|d2|*B^2, leaving room for third order),
    plus 1e-9 MHz of round-off. An injected 1% error in the Zeeman ramp violates
    this by a factor 90, so the bound is not merely absorbing the comparison.
    """
    B = 100.0
    for iso in ISOS:
        problem, p = parameters(iso)
        ctx = ctx_from(problem.spec, p)
        _, blocks = zero_field(iso)
        shift = {(r["mF"], r["index"]): r["E_MHz"] - r["E0_MHz"]
                 for r in level_records(iso, 0.0, B)}
        checked = resolves = 0
        for b in blocks:
            if b["mF"] == 0.0:
                continue                            # g is undefined at M_F = 0
            gf = g_factors(b["tm"], p, {"E_z": 0, "B_z": 0}, ctx=ctx)
            d2 = gf["d2"][("B_z", "B_z")]
            for i in b["ix"]:
                predicted = -gf["g"][i] * ctx.mu_B * B * b["mF"]
                bound = abs(d2[i]) * B ** 2 + 1e-9
                got = shift[(b["mF"], i)]
                assert abs(got - predicted) <= bound, (iso, b["mF"], i, got, predicted)
                # Resolving power: states where the bound is tight enough that a
                # 25% error in g would violate the assertion above.
                resolves += .25 * abs(predicted) > bound
                checked += 1
        assert checked > 0, iso
        # Not vacuous: the bound leaves room for a wrong g only where the exact
        # curvature is itself comparable to the linear term. Measured 80% (232),
        # 100% (227), 92% (229); a global 1% error in g is caught outright.
        assert resolves > .75 * checked, (iso, resolves, checked)


def test_parity_is_exact_at_zero_field_and_commutes_with_B_but_not_E():
    """The colour channel of the figure is zero-field parity, so it has to be a
    genuine eigenvalue: <P> = +/-1 exactly, which holds only if the zero-field
    Hamiltonian really is field free. The operator identities behind that are
    [H0, P] = [H_B, P] = 0 and {H_E, P} = 0 -- the electric dipole is parity odd,
    the magnetic moment parity even -- and they are the same assertions
    plot_thf_isotopes.py makes per block.

    Tolerances are relative to each operator's own largest element because the
    three operators differ by orders of magnitude in scale. These are exact
    algebraic identities, so the budget is round-off: measured 1.1e-15, and
    1e-12 allows a factor 1000 on that.
    """
    for iso in ISOS:
        problem, p = parameters(iso)
        _, blocks = zero_field(iso)
        assert set(np.unique(level_records(iso, 0.0, 0.0)["parity"])) == {-1, 1}, iso
        for b in blocks:
            assert np.max(np.abs(np.abs(b["parity"]) - 1)) < 1e-12, (iso, b["mF"])
            he = hamiltonian(b["tm"], p, {"E_z": 1, "B_z": 0}) - b["h0"]
            hb = hamiltonian(b["tm"], p, {"E_z": 0, "B_z": 1}) - b["h0"]
            for op, sign, name in ((b["h0"], -1, "H0"), (hb, -1, "H_B"), (he, +1, "H_E")):
                resid = float(np.max(np.abs(op @ b["P"] + sign * b["P"] @ op)))
                scale = float(np.max(np.abs(op)))
                assert resid < 1e-12 * scale, (iso, b["mF"], name, resid / scale)


def test_energy_reference_is_the_zero_field_centroid_and_the_same_in_every_config():
    """Panels are comparable across configurations only if the subtracted
    constant does not move with the field. The reference must be the mean
    zero-field energy of that J manifold, and the drawn ordinates must equal
    E_MHz minus exactly that constant in all four configurations."""
    for iso in ISOS:
        ref = centroids(iso)
        zero = level_records(iso, 0.0, 0.0)
        assert set(ref) == set(JS), iso
        for j in JS:
            want = float(zero["E0_MHz"][zero["J"] == j].mean())
            assert ref[j] == want, (iso, j, ref[j], want)
        for E_z, B_z in CONFIGS:
            rec = level_records(iso, E_z, B_z)
            fig, ax = plt.subplots()
            drawn = draw_panel(ax, rec, 2, ref[2])
            ordinates = sorted(float(line.get_ydata()[0]) for line in ax.lines)
            assert ordinates == sorted(float(r["E_MHz"]) - ref[2] for r in drawn), (iso, E_z, B_z)
            plt.close(fig)


def test_bars_are_drawn_at_true_M_F_with_no_F_dependent_offset():
    """The rejected layout shifted each F sideways, which makes the x axis
    misreport M_F. Every bar centre must be its record's M_F exactly, so bars
    belonging to different F at the same M_F share one abscissa. A 1e-12
    tolerance is twelve orders below the ~0.1 offset that layout would apply.

    Also checks the Scope A tick formatter is wired to the axis, since nothing
    else covers that join.
    """
    rec = level_records("229", 100.0, 100.0)
    fig, ax = plt.subplots()
    drawn = draw_panel(ax, rec, 1, centroids("229")[1])
    assert len(ax.lines) == len(drawn) > 0
    centres = np.sort([np.mean(line.get_xdata()) for line in ax.lines])
    assert np.max(np.abs(centres - np.sort(drawn["mF"]))) < 1e-12
    # Several F share each M_F here; an offset layout would give them distinct centres.
    assert len(set(np.round(centres, 9))) == len(set(drawn["mF"].tolist()))
    assert ax.xaxis.get_major_formatter()(1.5, 0) == "3/2"
    plt.close(fig)


def test_the_plotted_energies_are_converged_against_the_J_cutoff():
    """The J = 1-3 parents are drawn out of a J_max = 8 basis, so the figure is
    only meaningful if the cutoff has stopped mattering. Compared at the
    strongest configuration, where the Stark mixing into higher J is largest.

    Tolerance 1e-3 MHz is the 1 kHz plot target plot_thf_isotopes.py already
    holds its own J <= 7 versus J <= 8 comparison to; the observed J_max 6
    versus 8 difference here is ~1e-9 MHz, six orders inside it.
    """
    for iso in ISOS:
        coarse, fine = {}, {}
        for rec, into in ((level_records(iso, 100.0, 100.0, jmax=6), coarse),
                          (level_records(iso, 100.0, 100.0, jmax=8), fine)):
            for r in rec:
                into.setdefault((r["J"], r["F"], r["parity"], r["mF"]), []).append(r["E_MHz"])
        assert set(coarse) == set(fine), iso
        for key in fine:
            err = float(np.max(np.abs(np.sort(coarse[key]) - np.sort(fine[key]))))
            assert err < 1e-3, (iso, key, err)


def test_recovered_F1_lies_on_the_basis_lattice_and_couples_to_F():
    """The zoomed grid groups panels by the intermediate F1 = J + I_Th, so F1
    has to be a genuine label and not a relabelled F. Two independent
    constraints: every recovered F1 is a value the basis itself enumerates, and
    it couples to that record's F through the remaining spin,
    |F1 - I| <= F <= F1 + I. An off-by-one in the F1(F1+1) inversion leaves the
    lattice; reading the F column instead breaks the triangle rule wherever
    F1 != F.

    232Th has no second spin, so its basis carries no F1 column and F1 falls
    back to F; the lattice check still applies.
    """
    for iso in ISOS:
        problem = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=8)
        column = "F1" if "F1" in problem.kets.dtype.names else "F"
        lattice = set(np.unique(problem.kets[column]).tolist())
        rec = level_records(iso, 0.0, 0.0)
        assert set(np.unique(rec["F1"]).tolist()) <= lattice, iso
        I = problem.spec.I
        assert np.all(np.abs(rec["F1"] - I) <= rec["F"] + 1e-9), iso
        assert np.all(rec["F"] <= rec["F1"] + I + 1e-9), iso


def test_manifold_groups_partition_every_record_into_exactly_one_J_F1_panel():
    """The zoomed grid draws one panel per (J, F1), so the grouping has to be a
    partition: no record dropped, none drawn twice, and every panel homogeneous
    in both labels. Grouping on a float F1 is where this would break -- a
    tolerance-based match would merge adjacent F1, and a stale key would drop a
    manifold silently, since the figure gives no count to check against.
    """
    for iso in ISOS:
        rec = level_records(iso, 0.0, 0.0)
        groups = manifold_groups(rec)
        assert set(groups) == {(int(r["J"]), float(r["F1"])) for r in rec}, iso
        assert sum(len(g) for g in groups.values()) == len(rec), iso
        for (j, f1), g in groups.items():
            assert len(g) > 0 and np.all(g["J"] == j) and np.all(g["F1"] == f1), (iso, j, f1)


def test_zoomed_grid_gives_every_panel_in_a_J_row_one_shared_energy_window():
    """Two requirements meet in the y axis. Referencing each panel to its own
    (J, F1) zero-field centroid is what keeps it zoomed; sharing one window
    across a row is what keeps splittings comparable between the F1 manifolds of
    one J. Together they force identical ylim within a row.

    Rows must differ from each other or the check is vacuous -- a single global
    window would satisfy 'identical within a row' while zooming nothing. All
    three rows are required distinct, which no global window can satisfy, plus a
    factor 1.3 between the first and last: 227Th19F spans 32.5, 35.3 and 49.4 MHz
    for J = 1, 2, 3, a measured ratio of 1.52.
    """
    iso = "227"
    groups = manifold_groups(level_records(iso, 0.0, 100.0))
    fig, panels = draw_manifold_grid(iso, 0.0, 100.0)
    assert set(panels) == set(groups)
    rows = {}
    for key, ax in panels.items():
        rows.setdefault(key[0], []).append(ax.get_ylim())
        assert len(ax.lines) == len(groups[key]), key
    for j, lims in rows.items():
        assert len({tuple(np.round(lim, 9)) for lim in lims}) == 1, (j, lims)
    spans = [rows[j][0][1] - rows[j][0][0] for j in JS]
    assert len({round(s, 6) for s in spans}) == len(JS), spans
    assert spans[-1] > 1.3 * spans[0], spans
    plt.close(fig)


def test_each_zoomed_panel_reports_the_offset_its_own_window_hides():
    """Referencing every panel to its own centroid is what removes the
    information needed to compare manifolds, so each panel has to state the
    constant it subtracted. The test reads those numbers back off the titles and
    requires their differences to reproduce the true spacing between manifold
    centroids -- the quantity the labels exist to let a reader compute, rather
    than merely that some number is printed.

    229Th19F, whose F1 manifolds are split by the ~1.5 GHz Th hyperfine, so
    within-J and between-J spacings are both exercised. Tolerance 2 kHz: the
    labels carry 1 kHz (six decimals in GHz) and a difference of two of them
    inherits at most twice that rounding.
    """
    iso = "229"
    ref = manifold_reference(iso)
    fig, panels = draw_manifold_grid(iso, 0.0, 100.0)
    shown = {k: float(re.search(r"([\d.]+)\s*GHz", ax.get_title()).group(1)) * 1e3
             for k, ax in panels.items()}
    assert set(shown) == set(ref)
    keys = sorted(shown)
    for a, b in zip(keys, keys[1:]):
        assert abs((shown[a] - shown[b]) - (ref[a] - ref[b])) < 2e-3, (a, b)
    plt.close(fig)


def test_the_exploratory_parameter_set_is_the_one_actually_diagonalised():
    """The figures are made with the exploratory inputs, not package defaults:
    229 takes g_N_Th from Zitzer 2025 and omits the unknown Th quadrupole, 227
    takes the Minkov 2024 deformed-nucleus estimate. Read off the ParamSet that
    zero_field() actually built its Hamiltonians from, so a silent fallback to
    problem.params would fail here rather than quietly changing every figure.
    """
    for iso, expected in (("229", {"g_N_Th": .1460, "A_par_Th": -1519.495,
                                   "eQq0_Th": 0.0, "eQq2_Th": 0.0}),
                          ("227", {"g_N_Th": -.1720, "A_par_Th": 1790.176})):
        used = zero_field(iso)[0].params
        default = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=8).params.params
        for name, value in expected.items():
            assert used[name].value == value, (iso, name, used[name].value)
        # Guards the guard: these must actually differ from what the package ships.
        assert any(used[n].value != default[n].value for n in expected), iso
