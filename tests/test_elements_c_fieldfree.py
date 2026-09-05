"""Field-free case (c) elements against the closed forms of
docs/thf-plus-x3delta1-effective-hamiltonian.md S6 (V3, V4, V15) and S2.5.

Every check here is either an exact closed form or an exact structural
identity. None compares against a stored spectrum.
"""
import numpy as np
import pytest

from _helpers import _elem, _find
from heff import elements_c  # importing registers the six terms
from heff.params import thf_v1
from heff.spec import enumerate_kets, thf_spec
from heff.terms import check_selection_rules, ctx_from, terms_for_case


@pytest.fixture
def basis():
    return enumerate_kets(thf_spec())


@pytest.fixture
def ctx():
    return ctx_from(thf_spec(), thf_v1())


def test_A5_holds_for_every_registered_field_free_term(basis, ctx):
    """Gate A5 over the real registry, not toys.

    Catches, per term, the two failures a hand-written element table cannot
    self-report: a DEAD operator (formula zero everywhere inside its own
    declared rules -- a mistyped guard, a wrong field name) and a LEAKING one
    (non-zero where the rules say zero, which silently breaks the assembler's
    sparsity mask). Both outcomes are reachable: gate A5's PASS/FAIL pair is
    demonstrated on fabricated terms in tests/test_terms.py.
    """
    for t in terms_for_case("c"):
        report = check_selection_rules(t, basis, ctx)
        assert report["n_nonzero_inside"] > 0, f"{t.name} is a DEAD OPERATOR"
        assert report["n_nonzero_outside"] == 0, (
            f"{t.name} is non-zero outside its declared rules "
            f"(max {report['max_outside']:.3e})")


def test_ph_raises_on_a_half_integer_exponent():
    """(-1.0)**2.5 is complex in Python; _ph must raise instead of poisoning H."""
    assert elements_c._ph(2.0) == 1.0
    assert elements_c._ph(-3) == -1.0
    with pytest.raises(ValueError):
        elements_c._ph(2.5)


def test_every_term_carries_a_citation():
    """Every registered case (c) term names a primary source."""
    for t in terms_for_case("c"):
        assert len(t.cite) > 10 and any(k in t.cite for k in ("Ng", "B&C", "Leanhardt"))


def test_rotation_and_centrifugal_are_the_case_c_forms(basis, ctx):
    """No -Omega^2 term: Ng Ch. 2.4 Eq. 2.2a, [HAM] S2.1."""
    for J in (1, 2, 3, 4):
        i = _find(basis, J, 1.0, J + 0.5, 0.5)
        assert _elem("rotation", basis, ctx, i, i) == pytest.approx(J * (J + 1))
        assert _elem("centrifugal", basis, ctx, i, i) == pytest.approx(-(J * (J + 1.0)) ** 2)
    # Omega-independence is the case (c) signature: flip Omega, same energy
    a = _find(basis, 2, 1.0, 2.5, 0.5)
    b = _find(basis, 2, -1.0, 2.5, 0.5)
    assert _elem("rotation", basis, ctx, a, a) == _elem("rotation", basis, ctx, b, b)


def test_V3_hyperfine_J_scaling_and_the_nine_fifths_ratio(basis, ctx):
    """[HAM] V3: splitting = A_par (2J+1) / (2 J(J+1)); J=1 : J=2 is exactly 9/5.

    Uniquely catches dividing by 2J^2 instead of 2J(J+1) (Ng writes the operator
    as (F^2 - I^2 - J^2)/(2J^2)), and mis-coupling I to J. FAIL is reachable: a
    2J^2 denominator gives splitting (2J+1)/(2J^2), i.e. a J=1:J=2 ratio of
    12/5 = 2.4, not 9/5.
    """
    A = thf_v1().value("A_par")
    split = {}
    for J in (1, 2, 3, 4):
        up = _find(basis, J, 1.0, J + 0.5, 0.5)
        lo = _find(basis, J, 1.0, J - 0.5, 0.5)
        split[J] = A * (_elem("hyperfine_A_par", basis, ctx, up, up)
                        - _elem("hyperfine_A_par", basis, ctx, lo, lo))
        assert split[J] == pytest.approx(A * (2 * J + 1) / (2 * J * (J + 1)))
    assert split[1] / split[2] == pytest.approx(9.0 / 5.0)
    assert split[1] == pytest.approx(0.75 * A)
    assert split[2] == pytest.approx(A * 5.0 / 12.0)


def test_V15_cI_and_A_par_have_opposite_J_dependence(basis, ctx):
    """[HAM] V15: c_I splitting = c_I (2J+1)/2 GROWS with J; A_par's FALLS.

    Uniquely catches conflating the two hyperfine operators, which are
    indistinguishable at a single J. Both share the
    [F(F+1) - I(I+1) - J(J+1)] factor and differ only in the denominator, so a
    c_I element that kept A_par's 1/[2J(J+1)] passes at no J beyond a scale and
    FAILS the monotonicity and the 3.0 ratio here.
    """
    cI = thf_v1().value("c_I")
    sp = []
    for J in (1, 2, 3, 4):
        up = _find(basis, J, 1.0, J + 0.5, 0.5)
        lo = _find(basis, J, 1.0, J - 0.5, 0.5)
        s = cI * (_elem("spin_rotation_cI", basis, ctx, up, up)
                  - _elem("spin_rotation_cI", basis, ctx, lo, lo))
        assert s == pytest.approx(cI * (2 * J + 1) / 2.0)
        sp.append(abs(s))
    assert sp == sorted(sp) and sp[0] < sp[-1]          # grows with J
    assert sp[3] / sp[0] == pytest.approx(3.0)          # 90 kHz / 30 kHz


def test_V4_omega_doubling_splitting_law(basis, ctx):
    """[HAM] V4 (magnitude half): splitting = omega_ef J(J+1)/2, ratio J=2:J=1 = 3.

    Uniquely catches a J(J+1) that is really J^2 (ratio 4, not 3), and a factor
    2 or 4 slip between Ng Eq. C.3's [J(J+1)/2] and its 1/2 prefactor (which
    moves the J=1 splitting off 5.29 MHz).
    """
    w = thf_v1().value("omega_ef")
    split = {}
    for J in (1, 2, 3, 4):
        p = _find(basis, J, 1.0, J + 0.5, 0.5)
        m = _find(basis, J, -1.0, J + 0.5, 0.5)
        off = w * _elem("omega_doubling", basis, ctx, p, m)
        split[J] = 2.0 * abs(off)
        assert split[J] == pytest.approx(w * J * (J + 1) / 2.0)
    assert split[1] == pytest.approx(5.29)
    assert split[2] / split[1] == pytest.approx(3.0)
    assert split[4] == pytest.approx(52.90, abs=1e-9)


def test_omega_doubling_is_hermitian_and_only_connects_the_doublet(basis, ctx):
    """Symmetric in (bra, ket), zero on the diagonal, zero across F."""
    p = _find(basis, 2, 1.0, 2.5, 0.5)
    m = _find(basis, 2, -1.0, 2.5, 0.5)
    assert _elem("omega_doubling", basis, ctx, p, m) == pytest.approx(
        _elem("omega_doubling", basis, ctx, m, p))
    assert _elem("omega_doubling", basis, ctx, p, p) == 0.0
    other = _find(basis, 2, -1.0, 1.5, 0.5)
    assert _elem("omega_doubling", basis, ctx, p, other) == 0.0


def test_V4_upper_doublet_parity_alternates_as_minus_one_to_the_J(basis, ctx):
    """[HAM] V4 (parity half): the UPPER Omega-doublet component has parity (-1)^J.

    This is the convention-free statement, read off Ng 2022 Fig. 2 as an image
    ([HAM] S2.3): negative parity above at J = 1, positive parity above at J = 2;
    equivalently e above f uniformly, which is Gresh's k'' < 0.

    Uniquely catches a J-DEPENDENT phase on the Omega-doubling element in a
    package whose parity operator is E*|J,Om> = (-1)^(J-S+s)|J,-Om>: Ng Eq. C.3's
    (-1)^J prefactor, transcribed literally, passes at J = 1 and 3 and FAILS at
    J = 2 and 4 (the two alternations cancel, leaving upper parity -1 at every
    J). A globally flipped sign fails at odd J instead. Both outcomes reachable;
    see docs/open-questions.md OQ-A.
    """
    from heff.conventions import superposition_parity

    w = thf_v1().value("omega_ef")
    for J in (1, 2, 3, 4):
        p = _find(basis, J, 1.0, J + 0.5, 0.5)
        m = _find(basis, J, -1.0, J + 0.5, 0.5)
        c = w * _elem("omega_doubling", basis, ctx, p, m)
        sym_is_upper = c > 0                       # eigenvalues of c*sigma_x are +-c
        upper = superposition_parity(J, +1 if sym_is_upper else -1, S=ctx.S,
                                     ell=0.0, s=0.0)
        assert upper == int((-1) ** J), f"J={J}: upper parity {upper}"


def test_hyperfine_dJ1_does_not_connect_F_equal_J_minus_half_to_the_next_J(basis, ctx):
    """J = 1, F = 1/2 has no Delta-J = +-1 partner: no J = 2 ket shares F = 1/2."""
    i = _find(basis, 1, 1.0, 0.5, 0.5)      # J=1, F=1/2 = J - 1/2
    for j in np.flatnonzero(basis["J"] == 2.0):
        assert _elem("hyperfine_A_par_dJ1", basis, ctx, i, int(j)) == 0.0
    # and the partner that does exist is non-zero, so the test is not vacuous
    k = _find(basis, 2, 1.0, 1.5, 0.5)
    j1 = _find(basis, 1, 1.0, 1.5, 0.5)
    assert abs(_elem("hyperfine_A_par_dJ1", basis, ctx, j1, k)) > 0.0


def test_hyperfine_dJ1_matrix_element_matches_the_derived_value(basis, ctx):
    """Arithmetic cross-check against [HAM] S2.5's derived table: the
    J = 1 <-> 2, F = 3/2 element is +8.704 MHz at A_par = -20.1 MHz, Omega = +1.

    Not a spectrum snapshot: one matrix element against a value computed by hand
    from B&C Eq. (9.51) in the physics document. Catches a transcription slip in
    the radical (which of the four factors carries the +1) and an index swap
    between the upper and the lower J there: at J = 1 <-> 2, F = 3/2 the correct
    J = 2 gives factors (3, 5, 1, 1), radical 15 and +8.704 MHz, while J = 1
    gives (2, 4, 0, 2) -- the third factor vanishes and the element is 0.
    """
    A = thf_v1().value("A_par")
    i = _find(basis, 1, 1.0, 1.5, 0.5)
    j = _find(basis, 2, 1.0, 1.5, 0.5)
    assert A * _elem("hyperfine_A_par_dJ1", basis, ctx, i, j) == pytest.approx(8.704, abs=1e-3)
    assert _elem("hyperfine_A_par_dJ1", basis, ctx, i, j) == pytest.approx(
        _elem("hyperfine_A_par_dJ1", basis, ctx, j, i))
    # the other two rows of the same [HAM] S2.5 table, which pin the J
    # dependence of the (J^2 - Omega^2)^(1/2) / 2J prefactor
    for Jlo, want in ((2, 9.475), (3, 9.731)):
        F = Jlo + 0.5
        a = _find(basis, Jlo, 1.0, F, 0.5)
        b = _find(basis, Jlo + 1, 1.0, F, 0.5)
        assert A * _elem("hyperfine_A_par_dJ1", basis, ctx, a, b) == pytest.approx(
            want, abs=1e-3)


def test_hyperfine_dJ1_is_odd_in_omega(basis, ctx):
    """B&C 9.51 carries {a Lambda + (b_F + 2c/3) Sigma} = A_par / Omega, so the
    element flips sign with Omega -- the property that makes this term
    parity-EVEN once the (-1)^(J+J') from Delta-J = +-1 is included."""
    ip = _find(basis, 1, 1.0, 1.5, 0.5)
    jp = _find(basis, 2, 1.0, 1.5, 0.5)
    im = _find(basis, 1, -1.0, 1.5, 0.5)
    jm = _find(basis, 2, -1.0, 1.5, 0.5)
    assert _elem("hyperfine_A_par_dJ1", basis, ctx, ip, jp) == pytest.approx(
        -_elem("hyperfine_A_par_dJ1", basis, ctx, im, jm))
