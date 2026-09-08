"""Parameter units, provenance, and convention contracts.

Each value carries its unit, so reduced and physical dipoles cannot be mixed.
Parameter status remains metadata; it does not gate calculation.
"""
import pytest

from heff.conventions import Conventions
from heff.params import (DEBYE_TO_MHZ_PER_V_CM, GV_PER_CM_TO_MHZ_PER_E_CM,
                         MU_B, MU_N, Param, ParamSet, thf_v1, thf_v2)


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


def test_untagged_d_mf_raises():
    """[HAM] S2.7: an untagged d_mf would silently skip the origin check that
    catches a 20 % Stark error, whatever conventions.dipole_origin says."""
    with pytest.raises(ValueError, match="origin tag"):
        ParamSet({"d_mf": Param(3.37, "D")}, Conventions())
    assert thf_v1().value("d_mf") == pytest.approx(1.6964978, abs=1e-7)


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


def test_missing_symbol_returns_the_default_not_a_keyerror():
    assert thf_v1().value("no_such_knob", default=0.0) == 0.0
    with pytest.raises(KeyError):
        thf_v1().value("no_such_knob")


def test_table_names_every_symbol_with_its_status():
    text = thf_v1().table()
    for sym in ("B0", "D0", "omega_ef", "A_par", "d_mf", "G_par", "g_N", "c_I"):
        assert sym in text
    assert "estimate" in text and "measured" in text


def test_alpha_unit_is_registered_with_factor_one():
    assert Param(1.0, "(MHz/(V/cm))^2/MHz").canonical == pytest.approx(1.0)


def test_thf_v2_232_agrees_with_thf_v1_on_every_v1_symbol():
    """The parameter half of the master reduction: thf_v2('232') must be
    thf_v1() plus the Th knobs at zero plus the two-photon alphas -- never a
    parallel re-typed table that can drift from it."""
    v1 = thf_v1()
    v2 = thf_v2("232")
    for sym in v1.params:
        assert v2.value(sym) == pytest.approx(v1.value(sym))
        assert v2.params[sym].unit == v1.params[sym].unit
        assert v2.params[sym].status == v1.params[sym].status


def test_thf_v2_229_carries_the_documented_values_and_statuses():
    ps = thf_v2("229")
    expect = {
        "A_par_Th": (-1510, "MHz", "ab-initio"),
        "g_N_Th": (0.1464, "", "derived"),
        "eQq0_Th": (-2600, "MHz", "placeholder"),
        "eQq2_Th": (300, "MHz", "placeholder"),
        "c_I_Th": (0.0, "kHz", "held-fixed"),
        "Q_Th": (3.11, "e*b", "measured"),
    }
    for sym, (value, unit, status) in expect.items():
        p = ps.params[sym]
        assert p.value == pytest.approx(value)
        assert p.unit == unit
        assert p.status == status
        assert p.source
    assert "sign UNVERIFIED, gap G4" in ps.params["A_par_Th"].note
    assert "~1 kHz to ~1 MHz" in ps.params["c_I_Th"].note
    assert "UNVERIFIED normalisation bridge" in ps.params["eQq2_Th"].note
    assert "sqrt" in ps.params["eQq2_Th"].note.lower()


def test_thf_v2_227_A_par_is_a_labelled_placeholder():
    """The placeholder is the Schmidt single-particle value, not the
    mu(227Th) = mu(229Th) assumption -- the note must say so."""
    ps = thf_v2("227")
    a_par = ps.params["A_par_Th"]
    assert a_par.status == "placeholder"
    assert a_par.value == pytest.approx(39821, abs=1)
    assert "Schmidt" in a_par.note
    assert "docs/lit/lookup-227th-nuclear-moment.md" in a_par.note
    g_n = ps.params["g_N_Th"]
    assert g_n.status == "placeholder"
    assert g_n.value == pytest.approx(-3.826)
    assert "Schmidt" in g_n.note
    # eQq0_Th/eQq2_Th are structurally absent for I_Th = 1/2, not zero-valued
    assert "eQq0_Th" not in ps.params
    assert "eQq2_Th" not in ps.params


def test_thf_v2_229_a_par_th_sign_keyword_selects_the_trusted_calculation():
    """The sign of A_par_Th is a parameter choice made through thf_v2's
    keyword, not a conventions.py fork -- 229Th's Param is signed, and the
    227Th placeholder ignores the keyword entirely."""
    assert thf_v2("229").params["A_par_Th"].value == pytest.approx(-1510)
    assert thf_v2("229", a_par_th_sign="positive").params["A_par_Th"].value == \
        pytest.approx(1510)
    assert thf_v2("227", a_par_th_sign="positive").params["A_par_Th"].value == \
        pytest.approx(39821, abs=1)
    with pytest.raises(ValueError, match="a_par_th_sign"):
        thf_v2("229", a_par_th_sign="sideways")


def test_odd_thorium_transfers_do_not_claim_target_isotope_measurements():
    # A source-isotope error bar cannot be used as a transfer error bar.
    for isotope in ('229', '227'):
        ps = thf_v2(isotope)
        for name in ('B0', 'D0', 'omega_ef', 'A_par', 'd_mf', 'G_par'):
            p = ps.params[name]
            assert p.status == 'estimate'
            assert p.uncertainty is None
            assert p.value == thf_v1().params[name].value
            assert p.isotopologue == isotope + 'Th19F+'


def test_interval_uncertainty_propagates_to_B0_in_MHz():
    # Ng 2022 reports a 29.09733(4) GHz interval; dividing by four also divides sigma.
    assert thf_v1().params['B0'].uncertainty == pytest.approx(.00004 * 1000 / 4)


def test_unvalidated_quadrupole_transfer_has_no_calibrated_error_bar():
    for name in ('eQq0_Th', 'eQq2_Th'):
        assert thf_v2('229').params[name].uncertainty is None


def test_thorium_hyperfine_spread_is_not_a_complete_uncertainty():
    assert thf_v2("229").params["A_par_Th"].uncertainty is None
