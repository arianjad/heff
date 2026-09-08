import pytest

from heff.assemble import build_term_matrices
from heff.params import thf_v1
from heff.spec import block_by_mF, enumerate_kets, thf_spec
from heff.terms import REGISTRY, Rules, Term, ctx_from, terms_for_case


def test_omitting_names_preserves_the_complete_sorted_case_catalogue():
    expected = tuple(
        REGISTRY[name] for name in sorted(REGISTRY) if "c" in REGISTRY[name].cases
    )

    assert terms_for_case("c") == expected


def test_explicit_term_names_preserve_user_order():
    terms = terms_for_case("c", names=("rotation", "omega_doubling"))

    assert tuple(term.name for term in terms) == ("rotation", "omega_doubling")


def test_unknown_selected_term_fails():
    with pytest.raises(ValueError, match="unknown selected term"):
        terms_for_case("c", names=("rotation", "rotatoin"))


def test_incompatible_selected_term_fails():
    registry = {
        "c2_only": Term("c2_only", ("p",), ("c2",), Rules(), True, True,
                        "test", lambda bra, ket, ctx: 0.0)
    }

    with pytest.raises(ValueError, match="do not support case"):
        terms_for_case("c", names=("c2_only",), registry=registry)


def test_duplicate_selected_term_fails():
    with pytest.raises(ValueError, match="duplicate"):
        terms_for_case("c", names=("rotation", "rotation"))


def test_build_term_matrices_selects_requested_terms_in_user_order():
    spec = thf_spec()
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    block = kets[block_by_mF(kets).index[1.5]]

    matrices = build_term_matrices(
        block, ctx, term_names=("rotation", "omega_doubling")
    )

    assert matrices.names == ("rotation", "omega_doubling")
