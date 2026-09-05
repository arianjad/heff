"""The term registry: one decorated function plus one parameter entry per term.

No dispatch dicts. This is the single largest simplification versus
Molecule-Structure, where a state outside the pre-existing envelope needs new
keys in nine collect_* dicts (spec S3.2).

`rules` earns its place three ways: the assembler uses it as a sparsity mask so
only allowed (i, j) are evaluated (the dominant build cost); it answers "which
terms are non-zero in this basis?" without building anything; and it is one half
of gate A5, the formula-vs-declared-rules consistency check that catches dead
operators.

A term's `param` is a TUPLE of knob symbols whose PRODUCT is the coefficient --
('d_mf', 'E_z') for the Stark term, ('G_par', 'B_z') for the Zeeman one. Fields
and Hamiltonian parameters are therefore one flat catalogue (spec S2.1(2)), which
is what makes dH/d(knob) available for free (heff.assemble.vertex).
"""
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Rules:
    """Declared selection rules, as DATA, separate from the formula.

    dF1 is v2-only (spec-v2 S2.2): it is ignored -- not enforced -- when it
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

    No defaults for S, Lambda, I, mu_B, mu_N (spec S3.2): Molecule-Structure's
    `S = 1/2` keyword default is the entire reason S is never threaded there and
    every element silently evaluates at S = 1/2. A signature with no defaults
    turns that whole bug class into a TypeError at first call.

    `frame` is the StateSpec label carried through so the term-matrix manifest
    can record it; no element reads it (see StateSpec.frame).

    `spins` is the v2 coupled-nuclear-spin chain, carried through from
    StateSpec.spins (spec-v2 S2.2). Empty is the v1 default, so a v1 Ctx() is
    unchanged.
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


def terms_for_case(case, *, registry=REGISTRY):
    """Every registered term applicable to a coupling case, sorted by name."""
    return tuple(registry[n] for n in sorted(registry) if case in registry[n].cases)


def check_selection_rules(t, kets, ctx, *, tol=1e-12):
    """Gate A5: is the formula non-zero somewhere inside its declared rules,
    and exactly zero outside them?

    A rule set NARROWER than the formula fails loudly (non-zero outside). A rule
    set WIDER than the formula degrades the gate to "non-zero somewhere", which
    is still the check that catches the dead-operator case -- graceful in the
    direction that matters (spec S3.2).
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
