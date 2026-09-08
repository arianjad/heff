"""The term registry: one decorated function plus one parameter entry per term.

``rules`` serves as the assembler's sparsity mask, identifies terms that can be
nonzero in a basis, and makes declared selection rules available for validation.

A term's `param` is a TUPLE of knob symbols whose PRODUCT is the coefficient --
('d_mf', 'E_z') for the Stark term, ('G_par', 'B_z') for the Zeeman one. Fields
and Hamiltonian parameters form one flat catalogue, allowing
``heff.assemble.vertex`` to evaluate derivatives with respect to a knob.
"""
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Rules:
    """Declared selection rules, as DATA, separate from the formula.

    ``dF1`` applies only to two-spin kets. It is ignored when it
    is None (the default) or when the ket dtype has no F1 field, so a v1
    Rules() and a v1 KET_C basis are both untouched.
    """
    dJ: tuple = (0,)
    dOm: tuple = (0.0,)
    dF: tuple = (0,)
    dmF: tuple = (0,)
    dF1: tuple | None = None

    def allows(self, bra, ket):
        ok = (float(bra["J"] - ket["J"]) in tuple(float(x) for x in self.dJ)
              and float(bra["Om"] - ket["Om"]) in tuple(float(x) for x in self.dOm)
              and float(bra["F"] - ket["F"]) in tuple(float(x) for x in self.dF)
              and float(bra["mF"] - ket["mF"]) in tuple(float(x) for x in self.dmF))
        if ok and self.dF1 is not None and "F1" in bra.dtype.names:
            ok = float(bra["F1"] - ket["F1"]) in tuple(float(x) for x in self.dF1)
        return ok


@dataclass(frozen=True)
class Ctx:
    """Everything an element needs that is not a named knob.

    ``S``, ``Lambda``, ``I``, ``mu_B``, and ``mu_N`` have no defaults so every
    element receives its molecular and unit conventions explicitly.

    `frame` is the StateSpec label carried through so the term-matrix manifest
    can record it; no element reads it (see StateSpec.frame).

    ``spins`` is the coupled-nuclear-spin chain carried through from
    ``StateSpec.spins``. Its empty default represents a one-spin context.
    """
    S: float
    Lam: float
    I: float
    mu_B: float
    mu_N: float
    conventions: object
    frame: str
    spins: tuple = ()


def ctx_from(spec, pset):
    """Build the element context from a StateSpec and a ParamSet."""
    from .params import MU_B, MU_N

    es = spec.electronic[0]
    return Ctx(S=es.S, Lam=es.Lam, I=spec.I, mu_B=MU_B, mu_N=MU_N,
               conventions=pset.conventions, frame=spec.frame, spins=spec.spins)


@dataclass(frozen=True)
class Term:
    name: str
    param: tuple
    cases: tuple
    rules: Rules
    hermitian: bool
    real: bool
    cite: str
    fn: Callable


REGISTRY = {}


def term(*, name, param, cases, rules, hermitian, real, cite, registry=REGISTRY):
    """Register an element function. Returns the plain function unchanged."""
    def deco(fn):
        if name in registry:
            raise ValueError(f"term {name!r} is already registered")
        registry[name] = Term(name=name, param=tuple(param), cases=tuple(cases),
                              rules=rules, hermitian=hermitian, real=real,
                              cite=cite, fn=fn)
        return fn
    return deco


def terms_for_case(case, *, names=None, registry=REGISTRY):
    """Compatible registered terms, all sorted or an explicit ordered subset."""
    if names is None:
        return tuple(registry[n] for n in sorted(registry) if case in registry[n].cases)
    requested = tuple(names)
    if not requested:
        raise ValueError("empty selected term names; omit names for the complete case catalogue")
    unknown = [name for name in requested if name not in registry]
    if unknown:
        raise ValueError(f"unknown selected term(s) {unknown}; known: {sorted(registry)}")
    incompatible = [name for name in requested if case not in registry[name].cases]
    if incompatible:
        raise ValueError(f"term(s) {incompatible} do not support case {case!r}")
    if len(set(requested)) != len(requested):
        raise ValueError(f"duplicate selected term names: {requested}")
    return tuple(registry[name] for name in requested)


def check_selection_rules(t, kets, ctx, *, tol=1e-12):
    """Check whether a formula is nonzero inside its declared rules and zero
    outside them.

    A rule set NARROWER than the formula fails loudly (non-zero outside). A rule
    set WIDER than the formula reports the nonzero values it finds inside the
    declared region.
    """
    n_allowed = n_inside = n_outside = 0
    max_outside = 0.0
    for i in range(len(kets)):
        for j in range(len(kets)):
            val = complex(t.fn(kets[i], kets[j], ctx))
            allowed = t.rules.allows(kets[i], kets[j])
            n_allowed += int(allowed)
            if abs(val) > tol:
                if allowed:
                    n_inside += 1
                else:
                    n_outside += 1
                    max_outside = max(max_outside, abs(val))
    return {"n_allowed": n_allowed, "n_nonzero_inside": n_inside,
            "n_nonzero_outside": n_outside, "max_outside": max_outside}
