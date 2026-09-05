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
    "": 1.0,
}

STATUSES = frozenset({"measured", "ab-initio", "derived", "estimate",
                      "held-fixed", "stale", "unspecified"})


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
