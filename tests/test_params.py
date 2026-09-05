"""Parameter record: units in the data model, provenance on every value.

Uniquely catches the C2V-Molecules d_0 incident (digest wart #3): a
reduced-vs-physical dipole default living in one script while the physical
override lived in another, silently scaling every answer by 6.4 %. A unit
field on the value is what prevents it. Nothing here gates on `status`
(per Arian's standing rule: estimator parameters are ordinary parameters).
"""
import pytest

from heff.conventions import Conventions
from heff.params import (DEBYE_TO_MHZ_PER_V_CM, GV_PER_CM_TO_MHZ_PER_E_CM,
                         MU_B, MU_N, Param, ParamSet, thf_v1)


def test_param_of_is_a_one_argument_constructor():
    p = Param.of(7274.3325)
    assert (p.value, p.unit, p.status) == (7274.3325, "MHz", "unspecified")


def test_unit_conversions_are_the_documented_ones():
    assert Param(0.2426456, "cm-1").canonical == pytest.approx(7274.33, abs=0.01)
    assert Param(3.37, "D").canonical == pytest.approx(1.6964978, abs=1e-7)
    assert Param(3.897, "kHz").canonical == pytest.approx(3.897e-3)


def test_edm_shift_reproduces_the_documented_micro_hertz():
    """[HAM] S2.12: d_e = 1e-31 e.cm at E_eff = 35 GV/cm is 0.846 uHz."""
    shift_mhz = 1e-31 * Param(35.0, "GV/cm").canonical
    assert shift_mhz * 1e12 == pytest.approx(0.846, abs=0.002)  # MHz -> uHz


def test_canonical_is_idempotent_in_the_sense_that_matters():
    """Conversion happens once, at read time; the record itself never mutates."""
    p = Param(3.37, "D")
    first, second = p.canonical, p.canonical
    assert first == second
    assert (p.value, p.unit) == (3.37, "D")


def test_unknown_unit_and_unknown_status_raise():
    with pytest.raises(ValueError, match="unit"):
        Param(1.0, "furlongs").canonical
    with pytest.raises(ValueError, match="status"):
        Param(1.0, "MHz", status="probably-fine")


def test_conflicting_conventions_between_params_raise():
    ps = {"d_mf": Param(3.37, "D", convention="center_of_mass"),
          "d_other": Param(2.74, "D", convention="heavy_nucleus")}
    with pytest.raises(ValueError, match="convention"):
        ParamSet(ps, Conventions())


def test_thf_v1_carries_the_documented_values_and_statuses():
    ps = thf_v1()
    assert ps.value("B0") == pytest.approx(7274.3325)
    assert ps.value("D0") == pytest.approx(3.897e-3)
    assert ps.value("omega_ef") == pytest.approx(5.29)
    assert ps.value("A_par") == pytest.approx(-20.1)
    assert ps.value("d_mf") == pytest.approx(1.6964978, abs=1e-7)
    assert ps.value("G_par") == pytest.approx(0.04756)
    assert ps.value("g_N") == pytest.approx(5.25773)
    assert ps.value("c_I") == pytest.approx(0.020)
    assert ps.params["c_I"].status == "estimate"
    assert ps.params["A_par"].status == "measured"
    assert ps.params["G_par"].status == "derived"
    assert ps.value("d_e") == 0.0 and ps.value("k_TP") == 0.0


def test_constants_match_the_document():
    assert MU_B == pytest.approx(1.3996245)
    assert MU_N == pytest.approx(7.6225932e-4)
    assert MU_N / MU_B == pytest.approx(1 / 1836.15267, rel=1e-5)
    assert DEBYE_TO_MHZ_PER_V_CM == pytest.approx(0.5034118)
    assert GV_PER_CM_TO_MHZ_PER_E_CM == pytest.approx(2.417989242e17)


def test_with_overrides_does_not_mutate_the_original():
    ps = thf_v1()
    with pytest.warns(UserWarning):
        ps2 = ps.with_(A_par=-21.5)
    assert ps.value("A_par") == pytest.approx(-20.1)
    assert ps2.value("A_par") == pytest.approx(-21.5)
    assert ps2.params["A_par"].unit == "MHz"


def test_bare_float_override_warns_once():
    with pytest.warns(UserWarning, match="bare float"):
        thf_v1().with_(A_par=-21.5)


def test_missing_symbol_returns_the_default_not_a_keyerror():
    assert thf_v1().value("no_such_knob", default=0.0) == 0.0
    with pytest.raises(KeyError):
        thf_v1().value("no_such_knob")


def test_table_names_every_symbol_with_its_status():
    text = thf_v1().table()
    for sym in ("B0", "D0", "omega_ef", "A_par", "d_mf", "G_par", "g_N", "c_I"):
        assert sym in text
    assert "estimate" in text and "measured" in text
