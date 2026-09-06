"""Typed parameter records with units, provenance and a conventions block.

Units live in the data model (spec S3.4). Canonical internal unit is MHz;
a dipole is canonically MHz/(V/cm) and E_eff is canonically MHz per (e.cm), so
that every coefficient in the term catalogue is a plain product of canonical
values and every assembled Hamiltonian is in MHz.

`status` exists because Arian already needs it -- c_I is an estimate, G_par is
derived from a measurement, A_par is measured. The report layer surfaces
statuses; NOTHING gates on them.
"""
import warnings
from dataclasses import dataclass, replace
from typing import Mapping

# Physical constants and conversions
# docs/thf-plus-x3delta1-effective-hamiltonian.md S3; Molecule-Structure
# Source Code/molecule_parameters.py:24-31 carries the same cm-1 and Debye factors.
CM1_TO_MHZ = 29979.2458
DEBYE_TO_MHZ_PER_V_CM = 0.5034118
MU_B = 1.3996245          # MHz/G
MU_N = 7.6225932e-4       # MHz/G  (= MU_B / 1836.15267)
# 1 GV/cm x 1 e.cm = 1e9 eV; e/h = 2.417989242e14 Hz/V exactly (SI 2019),
# so 1 eV = 2.417989242e8 MHz and 1 GV/cm -> 2.417989242e17 MHz per (e.cm).
GV_PER_CM_TO_MHZ_PER_E_CM = 2.417989242e17

_TO_MHZ = {
    "MHz": 1.0,
    "GHz": 1e3,
    "kHz": 1e-3,
    "Hz": 1e-6,
    "cm-1": CM1_TO_MHZ,
    "D": DEBYE_TO_MHZ_PER_V_CM,
    "MHz/(V/cm)": 1.0,
    "GV/cm": GV_PER_CM_TO_MHZ_PER_E_CM,
    "MHz/(e cm)": 1.0,
    "(MHz/(V/cm))^2/MHz": 1.0,
    "": 1.0,
}

STATUSES = frozenset({"measured", "ab-initio", "derived", "estimate",
                      "held-fixed", "stale", "unspecified", "placeholder"})


@dataclass(frozen=True)
class Param:
    """One number with its unit, uncertainty, source, status and convention."""
    value: float
    unit: str = "MHz"
    uncertainty: float = None
    source: str = None
    status: str = "unspecified"
    isotopologue: str = None
    convention: str = None
    note: str = None

    def __post_init__(self):
        if self.status not in STATUSES:
            raise ValueError(f"status {self.status!r} not in {sorted(STATUSES)}")

    @classmethod
    def of(cls, value, **kw):
        """One-argument constructor: a scratch value is not a form-filling exercise."""
        return cls(float(value), **kw)

    @property
    def canonical(self):
        """The value in canonical units (MHz, MHz/(V/cm), MHz/(e cm), or 1)."""
        if self.unit not in _TO_MHZ:
            raise ValueError(f"unknown unit {self.unit!r}; known: {sorted(_TO_MHZ)}")
        return float(self.value) * _TO_MHZ[self.unit]


@dataclass(frozen=True)
class ParamSet:
    """symbol -> Param, plus the conventions block that the results are stamped with."""
    params: Mapping
    conventions: object

    def __post_init__(self):
        # A convention mismatch between two Params in one set is a structural
        # error, so it raises (spec S3.4 "what could go wrong" (iii)).
        seen = {p.convention for p in self.params.values() if p.convention}
        if len(seen) > 1:
            raise ValueError(f"conflicting convention tags in one ParamSet: {sorted(seen)}")
        # d_mf is ORIGIN-DEPENDENT because ThF+ is an ion: 3.37 D about the
        # centre of nuclear mass is 2.74 D about the Th nucleus, so a mistagged
        # dipole is a silent 20 % error on every Stark element ([HAM] S2.7).
        # Checked here rather than in the Stark term because Ctx carries no
        # ParamSet, and this catches it before any assembly can happen.
        d = self.params.get("d_mf")
        if d is not None and not d.convention:
            from .conventions import _ALLOWED
            raise ValueError(
                f"d_mf has no origin tag (Param.convention); tag it with one of "
                f"{_ALLOWED['dipole_origin']} matching the conventions block's "
                "dipole_origin, e.g. Param(3.37, 'D', convention='center_of_mass') "
                "-- an untagged d_mf silently skips the origin check ([HAM] S2.7)")
        if d is not None and d.convention != self.conventions.dipole_origin:
            raise ValueError(
                f"d_mf is tagged convention={d.convention!r} but the conventions "
                f"block says dipole_origin={self.conventions.dipole_origin!r}; the "
                "two origins differ by e.r(Th->c.m.) = 0.72 D ([HAM] S2.7)")

    def value(self, symbol, default=None):
        if symbol not in self.params:
            if default is None:
                raise KeyError(symbol)
            return float(default)
        return self.params[symbol].canonical

    def with_(self, **overrides):
        """Return a new ParamSet. A bare float is accepted with a warning."""
        new = dict(self.params)
        for symbol, val in overrides.items():
            if isinstance(val, Param):
                new[symbol] = val
            else:
                old = self.params.get(symbol)
                unit = old.unit if old is not None else "MHz"
                warnings.warn(
                    f"bare float for {symbol!r}: wrapping as Param({val}, {unit!r}, "
                    "status='unspecified'); pass a Param to keep provenance",
                    UserWarning, stacklevel=2)
                new[symbol] = replace(old, value=float(val), status="unspecified",
                                      uncertainty=None, source=None) if old is not None \
                    else Param.of(val)
        return ParamSet(new, self.conventions)

    def table(self):
        """Human-readable parameter table for a notebook cell."""
        w = max(len(s) for s in self.params)
        lines = [f"{'symbol'.ljust(w)}  {'value':>14}  {'unit':<12} {'status':<11} source"]
        for sym in sorted(self.params):
            p = self.params[sym]
            lines.append(f"{sym.ljust(w)}  {p.value:>14.7g}  {p.unit:<12} "
                         f"{p.status:<11} {p.source or ''}")
        return "\n".join(lines)


def thf_v1():
    """The v1 232Th19F+ X 3Delta1 parameter set.

    Every value, unit, status and source is copied from
    docs/thf-plus-x3delta1-effective-hamiltonian.md S3.
    """
    from .conventions import Conventions

    P = Param
    params = {
        "B0": P(7274.3325, "MHz", uncertainty=0.0010, status="measured",
                isotopologue="232Th19F+",
                source="Ng 2022 SII B p.2 (4B = 29.09733(4) GHz); Ng thesis S4.1.1 p.78"),
        "D0": P(3.897, "kHz", uncertainty=0.120, status="measured",
                source="Gresh 2016 Table 1, X3Delta1 D'' column (1.30(4)e-7 cm-1)"),
        "omega_ef": P(5.29, "MHz", uncertainty=0.05, status="measured",
                      source="Ng 2022 Table I p.5; Ng thesis Table 4.1 p.87"),
        "A_par": P(-20.1, "MHz", uncertainty=0.1, status="measured",
                   source="Ng 2022 Table I p.5",
                   note="may absorb 2 c_I .. 6 c_I ~ 40-120 kHz of unmodelled "
                        "nuclear spin-rotation -- [HAM] S2.6, OPEN-6"),
        "c_I": P(20.0, "kHz", status="estimate",
                 source="[HAM] S2.6: B&C Table 8.12 c2(19F, CsF) = 15.1 kHz at "
                        "B_v = 0.183782 cm-1, scaled linearly in B",
                 note="factor-of-3 uncertainty; ThF+ is open-shell, CsF is not"),
        "d_mf": P(3.37, "D", uncertainty=0.09, status="measured",
                  convention="center_of_mass",
                  source="Ng 2022 Table I p.5",
                  note="ORIGIN-DEPENDENT for an ion: 2.74 D w.r.t. the Th nucleus "
                       "(Skripnikov & Titov 2015 Table II) = 3.46 D at c.m."),
        "G_par": P(0.04756, "", uncertainty=0.002, status="derived",
                   source="Ng thesis p.87 (0.048(2) for g_F < 0); reproduces "
                          "|g_F=3/2| = 0.0149 exactly -- [HAM] S3",
                   note="sign of g_F is NOT measured; theory forces g_F < 0. OPEN-4"),
        "g_N": P(5.25773, "", status="held-fixed",
                 source="19F nuclear g-factor, Petrov et al. arXiv:1704.06631 Eq. 3"),
        "E_eff": P(35.0, "GV/cm", uncertainty=2.45, status="ab-initio",
                   source="Ng thesis Eq. C.8 discussion p.322 (JILA adopts ~35); "
                          "Skripnikov & Titov 2015 give 37.3(7 %), Denis 2015 35.2",
                   note="OPEN-14: 35.2 vs 37.3 GV/cm. Uncertainty is the quoted "
                        "+-7 % of the adopted value (2.45 GV/cm), not a spread "
                        "over the three calculations"),
        "W_TP": P(50.0, "kHz", uncertainty=3.5, status="ab-initio",
                  source="Skripnikov & Titov 2015 Table II p.8 (+-7 %)",
                  note="uncertainty is that quoted +-7 % of 50 kHz"),
        "d_e": P(0.0, "", status="held-fixed",
                 note="electron EDM in e.cm; 0 turns the PT-odd block off"),
        "k_TP": P(0.0, "", status="held-fixed",
                  note="scalar-pseudoscalar coupling, dimensionless; 0 turns it off"),
    }
    return ParamSet(params, Conventions())


_ISOTOPOLOGUES_V2 = ("232", "229", "227")


def thf_v2(isotopologue, *, a_par_th_sign="negative"):
    """The 232/229/227 ThF+ X 3Delta1 parameter set (spec-v2 S5).

    Every value, unit, uncertainty, status and source/note is copied
    verbatim from docs/superpowers/specs/2026-09-05-heff-v2-isotopologues-
    two-photon.md S5, except the 227Th A_par_Th/g_N_Th placeholders, which
    are Arian's Schmidt-moment ruling ([HAM] S9.6, OPEN-20) superseding the
    spec's mu(227Th) = mu(229Th) row.

    `isotopologue` is one of '232' | '229' | '227'. Every isotopologue
    carries thf_v1()'s shared 19F/rotational/Stark/EDM knobs unchanged (S5.1)
    plus the two-photon alphas (S5.4); '232' is spin-0 so its Th knobs are
    held at zero (a gate requires thf_v2('232') to agree with thf_v1() on
    every v1 symbol); '229' and '227' add the isotope-specific Th hyperfine/
    quadrupole knobs. eQq0_Th/eQq2_Th are STRUCTURALLY ABSENT for 227ThF+
    (I_Th = 1/2 has no rank-2 nuclear matrix element), not zero-valued.

    `a_par_th_sign` is one of 'negative' | 'positive'. The 229Th A_par_Th
    Param is a SIGNED physical value, and this keyword selects which ab
    initio calculation it is signed from: 'negative' (the default) trusts
    Skripnikov & Titov 2015, recommended by docs/lit/lookup-apar-th-sign-
    convention.md; 'positive' trusts Denis 2015. This is a parameter choice
    -- which calculation you trust -- not a convention, so it lives here in
    params, not in conventions ([HAM] S9.4.3). For '227' and '232' the
    keyword is accepted but ignored: 227ThF+'s A_par_Th placeholder already
    carries a physical sign ([HAM] S9.6), and 232Th has no nuclear spin.
    """
    from .conventions import Conventions

    if isotopologue not in _ISOTOPOLOGUES_V2:
        raise ValueError(
            f"isotopologue must be one of {_ISOTOPOLOGUES_V2}, got {isotopologue!r}")
    if a_par_th_sign not in ("negative", "positive"):
        raise ValueError(
            f"a_par_th_sign must be one of ('negative', 'positive'), "
            f"got {a_par_th_sign!r}")

    P = Param
    two_photon = {
        "alpha_K0_dOm0": P(1.0, "(MHz/(V/cm))^2/MHz", status="placeholder",
            source="[SPEC-v2] S5.4",
            note="no ThF+ two-photon polarisability exists in any source read "
                 "([2gamma] S3.0, S5 gap 1: no published Raman/two-photon rate "
                 "or line strength for a transition within X 3Delta1, for "
                 "either molecule, with a positive control confirming the "
                 "query shape); the value 1.0 exists to make the geometry "
                 "plottable"),
        "alpha_K2_dOm0": P(1.0, "(MHz/(V/cm))^2/MHz", status="placeholder",
            source="[SPEC-v2] S5.4",
            note="no ThF+ two-photon polarisability exists in any source read "
                 "([2gamma] S3.0, S5 gap 1); the value 1.0 exists to make the "
                 "geometry plottable"),
        "alpha_K2_dOm2": P(1.0, "(MHz/(V/cm))^2/MHz", status="placeholder",
            source="[SPEC-v2] S5.4",
            note="the Delta-Omega = +-2 channel, which requires an Omega = 0 "
                 "intermediate ([2gamma] S3.3); no ThF+ two-photon "
                 "polarisability exists in any source read"),
    }

    if isotopologue == "232":
        th = {
            "A_par_Th": P(0.0, "MHz", status="held-fixed",
                source="[SPEC-v2] S5.1: 232Th is spin-0, no Th hyperfine term"),
            "g_N_Th": P(0.0, "", status="held-fixed",
                source="[SPEC-v2] S5.1: 232Th is spin-0"),
            "eQq0_Th": P(0.0, "MHz", status="held-fixed",
                source="[SPEC-v2] S5.1: 232Th is spin-0, no Th quadrupole term"),
            "eQq2_Th": P(0.0, "MHz", status="held-fixed",
                source="[SPEC-v2] S5.1: 232Th is spin-0"),
            "c_I_Th": P(0.0, "kHz", status="held-fixed",
                source="[SPEC-v2] S5.1: 232Th is spin-0"),
        }
    elif isotopologue == "229":
        th = {
            "A_par_Th": P(-1510 if a_par_th_sign == "negative" else 1510,
                "MHz", uncertainty=60, status="ab-initio",
                source="Skripnikov & Titov 2015 Table II FINAL(ThF+) -4163 "
                       "(mu/mu_N) MHz and Denis 2015 +1833 MHz, both rescaled "
                       "to mu = 0.366(6) mu_N => -1524 / +1491 MHz; mean of "
                       "the two rescalings with a spread-based uncertainty "
                       "([TH] S2.2)",
                note="This value is SIGNED. sign UNVERIFIED, gap G4; add the "
                     "authors' 7 % in quadrature for a hard bar. The sign "
                     "selects which ab initio calculation is trusted: "
                     "'negative' (default) trusts Skripnikov & Titov 2015, "
                     "recommended by docs/lit/lookup-apar-th-sign-"
                     "convention.md; 'positive' trusts Denis 2015"),
            "g_N_Th": P(0.1464, "", uncertainty=0.0024, status="derived",
                source="mu(229Th)/I = 0.366(6)/(5/2) ([TH] S1.2, Porsev 2021 "
                       "arXiv:2107.14723)",
                note="the 1974 value 0.46(4) still in ENSDF is superseded and "
                     "must never be used to rescale a published A_par"),
            "eQq0_Th": P(-2600, "MHz", uncertainty=1000, status="estimate",
                source="HfF+ anchor: eQq0(177HfF+) = -2100 MHz (Petrov 2018, "
                       "CCSD(T)) x Q(229Th)/Q(177Hf) = 3.11/3.365 x R_el in "
                       "[1, 1.65] => -2 to -3.3 GHz ([TH] S4.3)",
                note="No ThF+ or ThO eQq0 is published, for any isotope or "
                     "state -- gap G2"),
            "eQq2_Th": P(300, "MHz", uncertainty=100, status="estimate",
                source="Petrov 2018 Eqs. (24)-(25) route with w(ThF+) = "
                       "G_par + 0.002319 = 0.0499 against w(HfF+) = 0.014 "
                       "=> ~200-400 MHz ([TH] S4.4)",
                note="Inherits the UNVERIFIED normalisation bridge of [TH] "
                     "S4/S9.4 (OPEN-17): the value inherits an unresolved "
                     "factor sqrt(2) AND a sign between the B&C (9.52) "
                     "q = +-2 normalisation and Petrov 2018 Eq. (23)"),
            "c_I_Th": P(0.0, "kHz", status="held-fixed",
                source="No value anywhere -- gap G3",
                note="[TH] S4.6 declines to pick one and brackets it at "
                     "~1 kHz to ~1 MHz: g_N(Th)/g_N(F) = 0.0279 pushes down, "
                     "the 2700x larger electronic hyperfine factor on Th "
                     "pushes up. OPEN-19"),
            "Q_Th": P(3.11, "e*b", uncertainty=0.02, status="measured",
                source="Porsev 2021, weighted average over four Th3+ states "
                       "([TH] S1.3)",
                note="Carried for provenance; the Hamiltonian consumes "
                     "eQq0_Th/eQq2_Th, not Q"),
        }
    else:  # '227'
        schmidt_note = (
            "PLACEHOLDER, Arian's ruling 2026-09-05 ([HAM] S9.6, OPEN-20): "
            "the Schmidt single-particle moment for the tentative ENSDF "
            "(1/2+) odd-neutron ground state, an s1/2 orbital (j = l + 1/2, "
            "l = 0), gives g_s(n) = -3.826 and mu_Schmidt = 1/2 g_s(n) = "
            "-1.913 mu_N, so g_N = mu/I = -3.826 and A_par = G_el x g_N = "
            "(-10408 MHz) x (-3.826) = +39821 MHz. Schmidt values for "
            "deformed actinides are typically wrong by a factor ~2 "
            "(quenched g_s^eff ~ 0.6 g_s^free would move it to ~+24 GHz), "
            "and the (1/2+) spin assignment is itself tentative (the 5/2+ "
            "level sits only 9.3 keV away). Recorded for comparison: the "
            "alternative mu(227Th) = mu(229Th) = 0.366 mu_N assumption gives "
            "g_N = 0.732 and A_par = -7619 MHz (~-7.62 GHz) -- opposite "
            "sign, 5x smaller. See docs/lit/lookup-227th-nuclear-moment.md: "
            "no measured or estimated mu(227Th) exists in Stone's "
            "compilations or the IAEA NDS moments database.")
        th = {
            "A_par_Th": P(39821, "MHz", status="placeholder",
                source="[HAM] S9.6; docs/lit/lookup-227th-nuclear-moment.md",
                note=schmidt_note),
            "g_N_Th": P(-3.826, "", status="placeholder",
                source="[HAM] S9.6; docs/lit/lookup-227th-nuclear-moment.md",
                note=schmidt_note),
            "c_I_Th": P(0.0, "kHz", status="held-fixed",
                source="gap G3, as for 229Th ([TH] S4.6)",
                note="[TH] S4.6 brackets it at ~1 kHz to ~1 MHz, as for "
                     "229Th. OPEN-19"),
            # eQq0_Th/eQq2_Th deliberately absent: a rank-2 nuclear operator
            # has no matrix element for I_Th = 1/2 -- structurally absent,
            # not zero-valued ([SPEC-v2] S5.3, [TH] S1.4).
        }

    params = {**thf_v1().params, **th, **two_photon}
    return ParamSet(params, Conventions(version="thf-v2"))
