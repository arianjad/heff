"""TOML must reproduce native ThF isotope calculations and their caveats."""
import numpy as np
import pytest
from heff import load_model
from heff.spec import thf_spec, enumerate_kets
from heff.params import thf_v2
from heff.elements_c2 import REGISTRY_C2
from heff.assemble import build_term_matrices, hamiltonian
from heff.terms import ctx_from


@pytest.mark.parametrize('isotope', ['229', '227'])
def test_bundled_odd_thorium_matches_native_two_spin_calculation(isotope):
    problem = load_model('thf_plus', isotope=f'{isotope}Th19F').problem(J_max=2)
    spec = thf_spec(isotope, J_max=2)
    pset = thf_v2(isotope)
    names = tuple(sorted(name for name, term in REGISTRY_C2.items()
                         if all(symbol in pset.params or symbol in ('E_z', 'B_z')
                                for symbol in term.param)))
    kets = enumerate_kets(spec)
    native = build_term_matrices(kets, ctx_from(spec, pset), case='c2',
                                 registry=REGISTRY_C2, term_names=names)
    assert np.array_equal(problem.kets, kets)
    assert problem.term_names == names
    for got, expected in zip(problem.term_matrices.mats, native.mats):
        np.testing.assert_array_equal(got, expected)
    for fields in ({'E_z': 0., 'B_z': 0.}, {'E_z': 20., 'B_z': .03}):
        np.testing.assert_allclose(problem.hamiltonian(**fields),
            hamiltonian(native, pset, fields), atol=1e-11, rtol=0)
    for symbol in {s for name in names for s in REGISTRY_C2[name].param} - {'E_z', 'B_z'}:
        assert problem.params.params[symbol] == pset.params[symbol]
    assert problem.params.conventions == pset.conventions


def test_227_has_no_nuclear_quadrupole_and_preserves_placeholder_status():
    problem = load_model('thf_plus', isotope='227Th19F').problem(J_max=1)
    assert not any('quadrupole' in name for name in problem.term_names)
    assert 'eQq0_Th' not in problem.params.params
    assert problem.params.params['A_par_Th'].status == 'placeholder'
    assert problem.params.params['g_N_Th'].status == 'placeholder'


def test_229_selection_does_not_change_default_232():
    before = load_model('thf_plus').problem(J_max=1).hamiltonian()
    load_model('thf_plus', isotope='229Th19F').validate()
    after = load_model('thf_plus').problem(J_max=1).hamiltonian()
    np.testing.assert_array_equal(before, after)


def test_232_toml_preserves_native_parameter_provenance():
    from heff.params import thf_v1
    actual = load_model('thf_plus').problem().params
    assert actual.params == thf_v1().params
