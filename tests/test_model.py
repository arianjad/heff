from dataclasses import FrozenInstanceError
import subprocess
import sys
import warnings

import numpy as np
import pytest

from heff.assemble import build_term_matrices, hamiltonian
from heff.engine import sweep
import heff.model as model_module
from heff.model import load_model
from heff.params import thf_v1
from heff.spec import enumerate_kets, thf_spec
from heff.terms import ctx_from


def _write_thf_model(tmp_path, *, model_id="temporary_thf", extra=""):
    path = tmp_path / f"{model_id}.toml"
    path.write_text(
        f'''schema_version = 1
model_id = "{model_id}"
default_manifold = "X3Delta1"
default_isotopologue = "232Th19F+"
{extra}

[isotopologues."232Th19F+"]

[[isotopologues."232Th19F+".spins]]
label = "19F"
I = 0.5
couple_to = "J"

[manifolds.X3Delta1]
backend = "case_c"

[manifolds.X3Delta1.electronic]
label = "X3Delta1"
Omega = 1.0
S = 1.0
Lambda = 2.0
T0 = 0.0

[manifolds.X3Delta1.basis]
J_min = 1
J_max = 4
M = "blocks"
frame = "rotating"

[manifolds.X3Delta1.parameters]
B0 = {{ value = 7274.3325, unit = "MHz" }}
D0 = {{ value = 3.897, unit = "kHz" }}
omega_ef = {{ value = 5.29, unit = "MHz" }}
A_par = {{ value = -20.1, unit = "MHz" }}
c_I = {{ value = 20.0, unit = "kHz" }}
d_mf = {{ value = 3.37, unit = "D", convention = "center_of_mass" }}
G_par = {{ value = 0.04756, unit = "" }}
g_N = {{ value = 5.25773, unit = "" }}
E_eff = {{ value = 35.0, unit = "GV/cm" }}
W_TP = {{ value = 50.0, unit = "kHz" }}
d_e = {{ value = 0.0, unit = "" }}
k_TP = {{ value = 0.0, unit = "" }}
''', encoding="utf-8")
    return path


def _replace(path, old, new):
    text = path.read_text(encoding="utf-8")
    assert old in text
    path.write_text(text.replace(old, new), encoding="utf-8")


def test_loaded_thf_problem_matches_the_low_level_hamiltonian(tmp_path):
    """Catches composition that changes basis order, terms, or coefficients."""
    path = _write_thf_model(tmp_path)
    model = load_model(path, manifold="X3Delta1", isotope="232Th19F+")
    problem = model.problem(J_max=2)

    spec = thf_spec(J_max=2)
    kets = enumerate_kets(spec)
    ctx = ctx_from(spec, thf_v1())
    legacy_tm = build_term_matrices(kets, ctx)
    legacy_H = hamiltonian(
        legacy_tm, thf_v1(), {"E_z": 20.0, "B_z": 0.01})

    assert np.array_equal(problem.kets, kets)
    assert problem.term_matrices.names == legacy_tm.names
    for got, expected in zip(problem.term_matrices.mats, legacy_tm.mats):
        assert np.array_equal(got, expected)
    assert np.array_equal(problem.hamiltonian(E_z=20.0, B_z=0.01), legacy_H)


def test_basis_override_changes_only_the_problem_and_composition_is_frozen(tmp_path):
    """Catches an override mutating the source definition or model/problem."""
    model = load_model(_write_thf_model(tmp_path))
    source_basis = dict(model.manifold.basis)

    problem = model.problem(J_max=2)

    assert problem.spec.J_range == (1, 2)
    assert dict(model.manifold.basis) == source_basis
    assert model.problem().spec.J_range == (1, 4)
    with pytest.raises(FrozenInstanceError):
        model.backend = None
    with pytest.raises(FrozenInstanceError):
        problem.spec = None


def test_describe_includes_execution_inputs_without_optional_metadata_warnings(tmp_path):
    """Catches inspectability depending on optional provenance metadata."""
    model = load_model(_write_thf_model(tmp_path))

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        description = model.describe()

    assert caught == []
    assert "backend: case_c" in description
    assert "J_max=4" in description
    assert "rotation" in description
    assert "B0" in description and "7274.333" in description


def test_validate_checks_structure_without_building_term_matrices(tmp_path, monkeypatch):
    """Catches validation accidentally paying the matrix-construction cost."""
    model = load_model(_write_thf_model(tmp_path))

    def forbidden_build(*args, **kwargs):
        raise AssertionError("validate built term matrices")

    monkeypatch.setattr(model_module, "build_term_matrices", forbidden_build)
    assert model.validate() is None


def test_validate_accepts_runtime_knobs_and_absent_optional_zero_parameters(tmp_path):
    """Catches E_z or optional d_e being misclassified as missing constants."""
    path = _write_thf_model(tmp_path)
    _replace(
        path,
        '[manifolds.X3Delta1]\nbackend = "case_c"',
        '[manifolds.X3Delta1]\nbackend = "case_c"\n'
        'terms = ["stark_z", "pt_odd_edm"]')
    _replace(path, 'd_e = { value = 0.0, unit = "" }\n', "")

    model = load_model(path)

    assert model.validate() is None
    assert model.problem().term_names == ("stark_z", "pt_odd_edm")


def test_validate_rejects_a_missing_active_term_parameter_with_full_path(tmp_path):
    """Catches an active term silently receiving zero for a required constant."""
    path = _write_thf_model(tmp_path)
    _replace(
        path,
        '[manifolds.X3Delta1]\nbackend = "case_c"',
        '[manifolds.X3Delta1]\nbackend = "case_c"\nterms = ["rotation"]')
    _replace(path, 'B0 = { value = 7274.3325, unit = "MHz" }\n', "")

    with pytest.raises(
            ValueError,
            match=r"manifolds\.X3Delta1\.parameters\.B0.*rotation"):
        load_model(path).validate()


def test_unknown_backend_error_names_its_full_path_and_source_file(tmp_path):
    """Catches backend lookup errors detached from their configuration key."""
    path = _write_thf_model(tmp_path)
    _replace(path, 'backend = "case_c"', 'backend = "case_z"')

    with pytest.raises(ValueError) as caught:
        load_model(path)

    message = str(caught.value)
    assert str(path.resolve()) in message
    assert "manifolds.X3Delta1.backend" in message
    assert "case_z" in message and "case_c" in message


def test_runtime_unknown_basis_key_names_its_full_path_and_source_file(tmp_path):
    """Catches a misspelled runtime truncation key losing path context."""
    path = _write_thf_model(tmp_path)
    model = load_model(path)

    with pytest.raises(ValueError) as caught:
        model.problem(J_stop=2)

    message = str(caught.value)
    assert str(path.resolve()) in message
    assert "manifolds.X3Delta1.basis.J_stop" in message


def test_invalid_runtime_basis_range_names_the_failing_key(tmp_path):
    """Catches backend range failures being reported without a useful key."""
    model = load_model(_write_thf_model(tmp_path))

    with pytest.raises(
            ValueError, match=r"manifolds\.X3Delta1\.basis\.J_min"):
        model.problem(J_min=5)


def test_unknown_selected_term_names_its_full_path(tmp_path):
    """Catches selected-term lookup errors detached from manifolds.<id>.terms."""
    path = _write_thf_model(tmp_path)
    _replace(
        path,
        '[manifolds.X3Delta1]\nbackend = "case_c"',
        '[manifolds.X3Delta1]\nbackend = "case_c"\nterms = ["rotatoin"]')

    with pytest.raises(
            ValueError, match=r"manifolds\.X3Delta1\.terms.*rotatoin"):
        load_model(path).validate()


def test_manifold_and_isotopologue_lookup_is_case_insensitive(tmp_path):
    """Catches configuration identifiers being matched case-sensitively."""
    model = load_model(
        _write_thf_model(tmp_path),
        manifold="x3delta1", isotope="232th19f+")

    assert model.manifold.id == "X3Delta1"
    assert model.isotopologue.id == "232Th19F+"


def test_problem_sweep_matches_the_existing_engine(tmp_path):
    """Catches the high-level sweep changing arrays or engine defaults."""
    problem = load_model(_write_thf_model(tmp_path)).problem(J_max=1)
    knobs = {"E_z": np.array([0.0, 20.0]), "B_z": np.array([0.0, 0.01])}

    got = problem.sweep(**knobs)
    expected = sweep(problem.term_matrices, problem.params, knobs)

    assert np.array_equal(got.evals, expected.evals)
    assert np.array_equal(got.evecs, expected.evecs)
    assert got.active_terms == expected.active_terms


def test_importing_model_does_not_import_or_initialize_backend_adapters():
    """Catches the high-level module defeating lazy backend loading."""
    code = (
        "import sys; import heff.model; "
        "from heff.backend_registry import list_backends; "
        "assert 'heff.backends' not in sys.modules; "
        "assert list_backends() == ()"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
