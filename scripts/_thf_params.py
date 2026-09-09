"""Exploratory ThF+ isotopologue parameter sets, shared by the plotting scripts.

Deliberately separate from the package defaults: the 229 baseline omits the
unknown Th quadrupole and 227 uses a deformed-nucleus theory estimate. One home
so both `plot_thf_isotopes.py` and `plot_thf_mf_resolved.py` plot the same inputs.
"""
from dataclasses import replace

from heff import load_model


def parameters(iso, scenario="baseline"):
    problem = load_model("thf_plus", isotope=f"{iso}Th19F").problem(J_max=8)
    p = problem.params
    if iso == "229":
        p = p.with_(g_N_Th=replace(p.params["g_N_Th"], value=.1460, uncertainty=.0012,
                   source="Zitzer et al. PRA111 L050802 (2025), Table I: mu=.365(3) mu_N / I=2.5",
                   note="Plot-specific update from the 2021 input."),
                   A_par_Th=replace(p.params["A_par_Th"], value=-1519.495,
                   source="Skripnikov & Titov PRA91 042504 (2015), Table II: -4163*.365 MHz",
                   note="Single-source negative branch; about 7% theory scale, no calibrated combined interval."))
        if scenario == "baseline":
            p = p.with_(**{name: replace(p.params[name], value=0, status="held-fixed",
                   uncertainty=None, source="Explicit omission for this plotting baseline",
                   note="Unknown physical quadrupole, NOT an estimate of zero; compare sensitivity figure.")
                   for name in ("eQq0_Th", "eQq2_Th")})
    if iso == "227":
        for name, value in (("g_N_Th", -.1720), ("A_par_Th", 1790.176)):
            p = p.with_(**{name: replace(p.params[name], value=value, uncertainty=None,
                status="estimate", source="Minkov et al. PRC110 034327 (2024), Table IV",
                note="mu=-.0860 mu_N, I=.5; g=mu/I, A=(-10408 MHz)*g. Nuclear-model uncertainty unquantified; no Coriolis/collective mixing correction.")})
    return problem, p
