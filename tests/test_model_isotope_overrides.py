"""Isotope selection must choose execution inputs without cross-isotope leakage."""
import numpy as np
import pytest
from heff import load_model


def model_file(tmp_path, override):
    p = tmp_path / 'isotopes.toml'
    p.write_text('''schema_version = 1
model_id = "isotopes"
default_manifold = "X"
default_isotopologue = "base"
[[isotopologues.base.spins]]
label = "F"
I = 0.5
couple_to = "J"
[[isotopologues.other.spins]]
label = "F"
I = 0.5
couple_to = "J"
''' + override + '''
[manifolds.X]
backend = "case_c"
terms = ["rotation"]
[manifolds.X.electronic]
label = "X"
Omega = 1.0
S = 1.0
Lambda = 2.0
[manifolds.X.basis]
J_min = 1
J_max = 1
M = "all"
frame = "lab"
[manifolds.X.parameters]
B0 = 2.0
D0 = 0.1
''', encoding='utf-8')
    return p


def test_isotope_parameters_and_terms_override_without_mutating_base(tmp_path):
    p = model_file(tmp_path, '''[isotopologues.other.manifolds.X]
terms = ["rotation", "centrifugal"]
[isotopologues.other.manifolds.X.parameters]
B0 = { value = 3000.0, unit = "kHz", status = "estimate", source = "test source" }
''')
    base = load_model(p)
    other = load_model(p, isotope='other')
    assert other.problem().params.params['B0'].canonical == 3.0
    assert other.problem().params.params['B0'].status == 'estimate'
    assert other.problem().params.params['B0'].source == 'test source'
    assert other.problem().term_names == ('rotation', 'centrifugal')
    assert base.problem().params.params['B0'].canonical == 2.0
    assert base.problem().term_names == ('rotation',)
    assert other.definition.manifolds['X'].parameters['B0'].raw == 2.0
    assert not np.array_equal(base.problem().hamiltonian(), other.problem().hamiltonian())


def test_unknown_override_field_is_rejected_instead_of_ignored(tmp_path):
    p = model_file(tmp_path, '[isotopologues.other.manifolds.X]\nparamters = { B0 = 4.0 }\n')
    with pytest.raises(ValueError, match='paramters'):
        load_model(p, isotope='other')


def test_override_unknown_manifold_is_rejected(tmp_path):
    p = model_file(tmp_path, '[isotopologues.other.manifolds.typo]\nbackend = "case_c"\n')
    with pytest.raises(ValueError, match='typo'):
        load_model(p, isotope='other')

