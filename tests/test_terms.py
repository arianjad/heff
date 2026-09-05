"""Gate A5 -- every term must be non-zero somewhere inside its declared
selection rules and exactly zero outside them.

Uniquely catches DEAD OPERATORS. The real one this gate exists for:
Molecule-Structure Source Code/matrix_elements.py:644
LambdaDoubling_q_even_aBJ gates on Delta-Lambda = +2q while its 3j
wigner_3j(J0,2,J1,-P0,2*q,P1) forces Delta-P = -2q, so with Delta-Sigma = 0 it
is identically zero everywhere -- and BaF A0, CaOH A000 and YbOH A000 all carry
a non-zero q doing nothing. QuantumStates.jl's version of the same operator uses
-2q and works. Both PASS and FAIL are reachable on real code; here they are
reachable on fabricated terms in a private registry, below.
"""
import numpy as np
import pytest

from heff.conventions import Conventions
from heff.params import MU_B, MU_N, thf_v1
from heff.spec import enumerate_kets, thf_spec
from heff.terms import (REGISTRY, Ctx, Rules, Term, check_selection_rules,
                        ctx_from, term, terms_for_case)


@pytest.fixture
def basis():
    return enumerate_kets(thf_spec())


@pytest.fixture
def ctx():
    return ctx_from(thf_spec(), thf_v1())


def test_ctx_has_no_defaults_for_physical_quantities():
    """Ctx() with no args raises TypeError; ctx_from() populates all fields."""
    with pytest.raises(TypeError):
        Ctx()  # S, Lam, I, mu_B, mu_N, conventions are all required
    c = ctx_from(thf_spec(), thf_v1())
    assert (c.S, c.Lam, c.I) == (1.0, 2.0, 0.5)
    assert (c.mu_B, c.mu_N) == (MU_B, MU_N)
    assert isinstance(c.conventions, Conventions)


def test_rules_allows_is_a_pure_predicate_on_two_kets(basis):
    """Rules.allows is a pure function of two kets' quantum numbers."""
    r = Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,))
    assert r.allows(basis[0], basis[0])
    other = basis[np.flatnonzero(basis["J"] != basis["J"][0])[0]]
    assert not r.allows(basis[0], other)


def test_registration_records_metadata_and_rejects_a_duplicate_name():
    """@term records a Term with the given metadata and rejects a duplicate name."""
    reg = {}

    @term(name="toy", param=("p",), cases=("c",), rules=Rules(), hermitian=True,
          real=True, cite="none", registry=reg)
    def toy(bra, ket, ctx):
        return 1.0

    assert set(reg) == {"toy"}
    t = reg["toy"]
    assert isinstance(t, Term) and t.param == ("p",) and t.cite == "none"
    assert t.fn is toy  # the decorator returns the plain function

    with pytest.raises(ValueError, match="already registered"):
        term(name="toy", param=("p",), cases=("c",), rules=Rules(), hermitian=True,
             real=True, cite="none", registry=reg)(toy)


def test_A5_passes_for_an_honest_term(basis, ctx):
    """PASS side of gate A5: a formula that matches its declared rules is
    non-zero somewhere inside them and zero everywhere outside.

    Companion FAIL demos below (test_A5_catches_a_dead_operator,
    test_A5_catches_a_term_that_leaks_outside_its_rules) show the two failure
    modes this gate catches: a dead operator (rules and formula disagree, so
    the formula is zero everywhere) and a leaky operator (formula wider than
    its declared rules).
    """
    reg = {}

    @term(name="honest_rotation", param=("B0",), cases=("c",),
          rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
          hermitian=True, real=True, cite="toy", registry=reg)
    def honest(bra, ket, ctx):
        if (bra["J"], bra["Om"], bra["F"], bra["mF"]) != (
                ket["J"], ket["Om"], ket["F"], ket["mF"]):
            return 0.0
        return ket["J"] * (ket["J"] + 1.0)

    report = check_selection_rules(reg["honest_rotation"], basis, ctx)
    assert report["n_nonzero_inside"] > 0
    assert report["n_nonzero_outside"] == 0


def test_A5_catches_a_dead_operator(basis, ctx):
    """FAIL demo 1: the Molecule-Structure LambdaDoubling_q bug, transplanted.

    The rules declare Delta-Omega = +-2, but the formula's selector and its
    3j-analogue disagree in sign, so the element is identically zero.
    """
    reg = {}

    @term(name="dead_doubling", param=("omega_ef",), cases=("c",),
          rules=Rules(dJ=(0,), dOm=(-2.0, 2.0), dF=(0,), dmF=(0,)),
          hermitian=True, real=True, cite="toy", registry=reg)
    def dead(bra, ket, ctx):
        if bra["Om"] != ket["Om"] + 2.0:      # selector says +2
            return 0.0
        if bra["Om"] != ket["Om"] - 2.0:      # the 3j needs -2: never both
            return 0.0
        return 1.0

    report = check_selection_rules(reg["dead_doubling"], basis, ctx)
    assert report["n_allowed"] > 0
    assert report["n_nonzero_inside"] == 0     # <-- the gate's failing condition


def test_A5_catches_a_term_that_leaks_outside_its_rules(basis, ctx):
    """FAIL demo 2: a formula wider than its declared rules."""
    reg = {}

    @term(name="leaky", param=("B0",), cases=("c",),
          rules=Rules(dJ=(0,), dOm=(0.0,), dF=(0,), dmF=(0,)),
          hermitian=True, real=True, cite="toy", registry=reg)
    def leaky(bra, ket, ctx):
        if abs(bra["J"] - ket["J"]) > 1.0:
            return 0.0
        return 1.0

    report = check_selection_rules(reg["leaky"], basis, ctx)
    assert report["n_nonzero_outside"] > 0
    assert report["max_outside"] == pytest.approx(1.0)


def test_terms_for_case_filters_and_is_sorted():
    """terms_for_case returns only matching-case terms, sorted by name."""
    reg = {}
    for name, cases in (("b_only", ("b",)), ("c_only", ("c",)), ("both", ("a", "c"))):
        term(name=name, param=("p",), cases=cases, rules=Rules(), hermitian=True,
             real=True, cite="toy", registry=reg)(lambda bra, ket, ctx: 0.0)
    got = [t.name for t in terms_for_case("c", registry=reg)]
    assert got == ["both", "c_only"]
