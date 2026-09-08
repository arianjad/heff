from importlib.resources import files


def test_thf_toml_is_packaged():
    """Catches a source or wheel installation missing the bundled model."""
    resource = files("heff").joinpath("models/thf_plus.toml")
    assert resource.is_file()
    assert 'model_id = "thf_plus"' in resource.read_text(encoding="utf-8")
