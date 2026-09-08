from heff.backends import load_bundled_backends
from heff.backend_registry import get_backend
from heff.model_io import read_model_toml, resolve_param_set
import pytest


def _minimal_model(parameters, *, top_metadata=""):
    return f'''schema_version = 1
model_id = "minimal"
default_manifold = "X"
default_isotopologue = "test"
{top_metadata}

[isotopologues.test]

[[isotopologues.test.spins]]
label = "19F"
I = 0.5
couple_to = "J"

[manifolds.X]
backend = "case_c"

[manifolds.X.electronic]
label = "X3Delta1"
Omega = 1.0
S = 1.0
Lambda = 2.0

[manifolds.X.basis]
J_min = 1
J_max = 1
M = "blocks"
frame = "rotating"

[manifolds.X.parameters]
{parameters}
'''


def test_bare_numbers_need_no_source_status_or_unit(tmp_path):
    path = tmp_path / "minimal.toml"
    path.write_text(_minimal_model("B0 = 7274.3325"), encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()
    params = resolve_param_set(definition.manifolds["X"], get_backend("case_c"))
    assert params.params["B0"].value == 7274.3325
    assert params.params["B0"].unit == "MHz"
    assert params.params["B0"].source is None


def test_rich_dipole_without_source_uses_explicit_unit_and_convention(tmp_path):
    path = tmp_path / "rich.toml"
    path.write_text(
        _minimal_model('d_mf = { value = 3.37, unit = "D", convention = "center_of_mass" }'),
        encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()
    params = resolve_param_set(definition.manifolds["X"], get_backend("case_c"))

    assert params.value("d_mf") == pytest.approx(1.6964978, abs=1e-7)
    assert params.params["d_mf"].source is None


def test_explicit_unsupported_unit_fails_when_parameters_are_resolved(tmp_path):
    path = tmp_path / "unsupported-unit.toml"
    path.write_text(_minimal_model('B0 = { value = 1.0, unit = "bananas" }'),
                    encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()

    with pytest.raises(ValueError, match="bananas"):
        resolve_param_set(definition.manifolds["X"], get_backend("case_c"))


def test_unknown_optional_status_warns_and_normalizes_to_unspecified(tmp_path):
    path = tmp_path / "unknown-status.toml"
    path.write_text(_minimal_model('B0 = { value = 1.0, status = "unreviewed" }'),
                    encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()

    with pytest.warns(UserWarning, match="unknown optional status"):
        params = resolve_param_set(definition.manifolds["X"], get_backend("case_c"))
    assert params.params["B0"].status == "unspecified"


def test_unknown_parameter_metadata_is_retained_without_becoming_a_load_gate(tmp_path):
    path = tmp_path / "descriptive-metadata.toml"
    path.write_text(_minimal_model('B0 = { value = 1.0, calibration = "run-17" }'),
                    encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()

    assert definition.manifolds["X"].parameters["B0"].metadata["calibration"] == "run-17"
    with pytest.warns(UserWarning, match="preserving unknown metadata"):
        params = resolve_param_set(definition.manifolds["X"], get_backend("case_c"))
    assert params.value("B0") == 1.0


def test_recognized_convention_override_is_applied_before_dipole_origin_check(tmp_path):
    path = tmp_path / "heavy-nucleus.toml"
    path.write_text(
        _minimal_model(
            'd_mf = { value = 2.74, unit = "D", convention = "heavy_nucleus" }')
        + '''
[conventions]
dipole_origin = "heavy_nucleus"
''', encoding="utf-8")
    definition = read_model_toml(path)
    load_bundled_backends()
    params = resolve_param_set(definition.manifolds["X"], get_backend("case_c"))

    assert params.conventions.dipole_origin == "heavy_nucleus"


def test_unknown_top_level_descriptive_metadata_is_retained_and_does_not_block_load(tmp_path):
    path = tmp_path / "display-name.toml"
    path.write_text(
        _minimal_model("B0 = 1.0", top_metadata='display_name = "Minimal model"'),
        encoding="utf-8")

    definition = read_model_toml(path)

    assert definition.metadata["display_name"] == "Minimal model"


def test_missing_case_c_basis_key_names_its_toml_path(tmp_path):
    path = tmp_path / "missing-j-max.toml"
    path.write_text(
        '''schema_version = 1
model_id = "minimal"
default_manifold = "X"
default_isotopologue = "test"

[isotopologues.test]
spins = []

[manifolds.X]
backend = "case_c"

[manifolds.X.electronic]
label = "X3Delta1"
Omega = 1.0
S = 1.0
Lambda = 2.0

[manifolds.X.basis]
J_min = 1
M = "blocks"
frame = "rotating"

[manifolds.X.parameters]
B0 = 7274.3325
''', encoding="utf-8")

    with pytest.raises(ValueError, match=r"manifolds\.X\.basis\.J_max"):
        read_model_toml(path)
