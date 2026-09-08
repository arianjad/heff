"""Amide TOML execution: equivalent protons and a separate metal nuclear spin."""
import numpy as np
import pytest

from heff import load_model
from heff.params import MU_N


def amide_model(tmp_path, *, metal=0.5, terms=('hyperfine_M_contact',),
                parameters='a_M = 1.0', nmax=0, kvals=(0,)):
    path = tmp_path / 'amide.toml'
    path.write_text(f'''schema_version = 1
model_id = "synthetic_amide"
default_manifold = "X"
default_isotopologue = "test"

[[isotopologues.test.spins]]
label = "nitrogen"
I = 1.0
couple_to = "J"
[[isotopologues.test.spins]]
label = "protons"
I = 0.5
equivalent_count = 2
couple_to = "F_N"
[[isotopologues.test.spins]]
label = "metal"
I = {metal}
couple_to = "F_core"

[manifolds.X]
backend = "amide_c2v"
terms = {list(terms)!r}
[manifolds.X.electronic]
S = 0.5
vibronic_sign = 1
[manifolds.X.basis]
N_min = 0
N_max = {nmax}
K = {list(kvals)!r}
M = "all"
frame = "lab"
[manifolds.X.parameters]
{parameters}
''', encoding='utf-8')
    return path


@pytest.mark.parametrize('metal', [0.5, 3.5, 4.5])
def test_metal_contact_has_the_electron_nucleus_coupled_spectrum(tmp_path, metal):
    # N=K=0 => proton singlet. Nitrogen I=1 is a spectator of S.I_M.
    problem = load_model(amide_model(tmp_path, metal=metal)).problem()
    expected = np.repeat([-(metal+1)/2, metal/2],
                         [int(3*(2*metal)), int(3*(2*metal+2))])
    np.testing.assert_allclose(np.linalg.eigvalsh(problem.hamiltonian()), expected,
                               atol=1e-12, rtol=0)


def test_metal_zeeman_contains_every_metal_projection(tmp_path):
    problem = load_model(amide_model(tmp_path, metal=1.5,
        terms=('zeeman_metal',), parameters='g_M = 1.0')).problem()
    expected = np.sort(np.repeat(-MU_N * np.array([-1.5, -0.5, 0.5, 1.5]), 6))
    np.testing.assert_allclose(np.linalg.eigvalsh(problem.hamiltonian(B_z=1)),
                               expected, atol=1e-14, rtol=0)


def test_amide_basis_keeps_exchange_allowed_proton_sectors(tmp_path):
    problem = load_model(amide_model(tmp_path, metal=0.0,
        terms=('rotation_A',), parameters='A = 1.0', nmax=1, kvals=(0, 1))).problem()
    kets = problem.kets
    assert np.all(kets['I_T'][kets['K'] == 0] == 0)
    assert np.all(kets['I_T'][abs(kets['K']) == 1] == 1)
    assert problem.hamiltonian().shape == (len(kets), len(kets))


def test_amide_missing_selected_metal_parameter_is_not_silently_zero(tmp_path):
    with pytest.raises(ValueError, match='parameters.a_M'):
        load_model(amide_model(tmp_path, parameters='')).validate()


def test_amide_requires_explicit_equivalent_pair_declaration(tmp_path):
    path = amide_model(tmp_path)
    path.write_text(path.read_text().replace('equivalent_count = 2', 'equivalent_count = 1'))
    with pytest.raises(ValueError, match='equivalent'):
        load_model(path).validate()


def test_metal_coupled_matrices_match_uncoupled_spin_operators(tmp_path):
    # Independent product basis |m_S,m_N,m_M>, with N=K=I_T=0.
    # This checks signed matrix elements, not only invariant eigenvalues.
    from itertools import product
    from sympy import Rational
    from sympy.physics.wigner import clebsch_gordan
    from heff.params import MU_B

    def cg(a, b, c, ma, mb, mc):
        return float(clebsch_gordan(*(Rational(round(2*x), 2)
                                     for x in (a, b, c, ma, mb, mc))))

    states = list(product((-.5, .5), (-1., 0., 1.), (-.5, .5)))
    problem = load_model(amide_model(tmp_path,
        terms=('hyperfine_M_contact', 'zeeman_metal', 'zeeman_electron'),
        parameters='a_M = 1.0\ng_M = 1.0\ng_s = 1.0')).problem()
    U = np.array([[cg(.5, 1, q['F_N'], ms, mn, ms+mn)
                   * cg(q['F_core'], .5, q['F'], ms+mn, mm, q['mF'])
                   for q in problem.kets] for ms, mn, mm in states])
    np.testing.assert_allclose(U.T @ U, np.eye(len(states)), atol=1e-14)
    contact = np.zeros((len(states), len(states)))
    for col, (ms, mn, mm) in enumerate(states):
        contact[col, col] = ms * mm
        for ds in (-1, 1):
            dest = (ms+ds, mn, mm-ds)
            if dest in states:
                contact[states.index(dest), col] = .5 * np.sqrt(
                    (.75-ms*(ms+ds)) * (.75-mm*(mm-ds)))
    expected = {
        'hyperfine_M_contact': contact,
        'zeeman_metal': np.diag([-MU_N*s[2] for s in states]),
        'zeeman_electron': np.diag([MU_B*s[0] for s in states]),
    }
    matrices = problem.term_matrices
    for name, actual in zip(matrices.names, matrices.mats):
        np.testing.assert_allclose(actual, U.T @ expected[name] @ U,
                                   atol=1e-13, rtol=0)


def test_proton_triplet_core_tensors_match_explicit_spectator_cg_sum():
    from math import sqrt
    from heff.amide_basis import AmideBasisSpec, enumerate_amide_kets
    from heff.backends.amide import AmideContext, AmideConventions, REGISTRY
    from heff.wigner import w3j
    from dataclasses import replace

    spec = AmideBasisSpec(S=.5, I_N=1, i_H=.5, I_M=.5,
        N_range=(1, 1), K_values=(1,), vibronic_sign=1)
    kets = enumerate_amide_kets(spec)
    kets = kets[kets['mF'] == 0]
    ctx = AmideContext(spec, AmideConventions())
    core_ctx = AmideContext(replace(spec, I_M=0), ctx.conventions)

    def cg(g, f, mg, mi):
        return (-1)**int(g-.5) * sqrt(2*f+1) * w3j(g, .5, f, mg, mi, 0)

    for term in REGISTRY.values():
        if term.rules.rank != 1 or term.name == 'zeeman_metal':
            continue
        for bra in kets:
            for ket in kets:
                expected = 0.0
                for mi in (-.5, .5):
                    b, k = bra.copy(), ket.copy()
                    b['F'], k['F'] = b['F_core'], k['F_core']
                    b['mF'] = k['mF'] = -mi
                    expected += (cg(b['F'], bra['F'], -mi, mi)
                                 * cg(k['F'], ket['F'], -mi, mi)
                                 * term.fn(b, k, core_ctx))
                assert term.fn(bra, ket, ctx) == pytest.approx(expected, abs=1e-13)
