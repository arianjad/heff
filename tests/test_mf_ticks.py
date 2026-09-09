"""M_F axis tick formatting: physicists write M_F as an integer or a reduced
half-integer fraction, never as a bare float (Scope A of
docs/superpowers/specs/2026-09-09-mf-resolved-level-plots.md)."""
import pytest

from scripts._mf_ticks import format_mF


def test_integers_render_with_no_decimal_point():
    """M_F = 0, +/-1 are whole quanta; the label must not carry a trailing
    ".0" the way str(float) would produce."""
    assert format_mF(0.0) == "0"
    assert format_mF(1.0) == "1"
    assert format_mF(-1.0) == "-1"


def test_half_integers_render_as_fractions():
    """Odd multiples of 1/2 (F half-integer, so M_F is too) print as a
    fraction, never as a decimal like "0.5"."""
    assert format_mF(0.5) == "1/2"
    assert format_mF(-0.5) == "-1/2"
    assert format_mF(1.5) == "3/2"
    assert format_mF(-3.5) == "-7/2"


def test_half_integer_fraction_is_always_odd_over_two():
    """The reduced fraction's numerator is odd and its denominator is fixed
    at 2 -- M_F is always a multiple of 1/2, so there is no other reduced
    form (never a mixed number "2 1/2", never a wrongly-reduced "5/4")."""
    assert format_mF(2.5) == "5/2"
    assert format_mF(-2.5) == "-5/2"


def test_non_half_integer_raises_rather_than_mangling():
    """M_F is always a multiple of 1/2 (integer for integer F, half-integer
    for half-integer F). A value off that lattice means an upstream bug
    (wrong quantum number, unit mix-up); silently rounding it would draw a
    tick that misreports the state, so this raises instead of guessing."""
    with pytest.raises(ValueError, match="mF"):
        format_mF(0.3)


def test_float_dirt_from_arithmetic_still_snaps_to_the_right_label():
    """A M_F value arriving via eigh/np.linalg arithmetic picks up ~1e-12
    dirt (e.g. 0.5 + 1e-13), not the exact literal. Tolerance policy: snap
    to the nearest multiple of 1/2 within 1e-6 (far above float64 roundoff
    in this pipeline, far below the 0.5 spacing between valid values), so
    such dirt still reads as the intended physical value."""
    assert format_mF(0.5 + 1e-10) == "1/2"
    assert format_mF(-1.5 - 1e-10) == "-3/2"


def test_mf_formatter_matches_format_mF_on_a_representative_tick_list():
    """The matplotlib formatter is a thin wrapper: axis ticks must read the
    same as format_mF for a mixed integer/half-integer J=3 tick set, or the
    figure and the plain-Python label would silently disagree."""
    from scripts._mf_ticks import mf_formatter

    ticks = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
    fmt = mf_formatter()
    assert [fmt(t, i) for i, t in enumerate(ticks)] == [format_mF(t) for t in ticks]
