"""Render M_F quantum numbers as a physicist writes them: integers with no
decimal point, half-integers as a reduced odd/2 fraction. Matplotlib is
imported lazily inside the formatter factory (heff/__init__.py docstring:
importing the package -- and this module -- must not pull in Matplotlib).
"""


def format_mF(value):
    """Format one M_F value, e.g. 1.0 -> "1", -0.5 -> "-1/2", 1.5 -> "3/2".

    Raises ValueError if value is not a multiple of 1/2 -- M_F never takes
    any other value, so anything else is an upstream bug, not a display
    corner case.
    """
    twice = round(value * 2)
    if abs(twice - value * 2) > 1e-6:
        raise ValueError(f"mF={value!r} is not a multiple of 1/2")
    if twice % 2 == 0:
        return str(twice // 2)
    return f"{twice}/2"


def mf_formatter():
    """Build a matplotlib tick formatter for an M_F axis, e.g.
    ax.xaxis.set_major_formatter(mf_formatter()). Matplotlib is imported
    here, not at module scope, so importing this module stays light.
    """
    from matplotlib.ticker import FuncFormatter

    return FuncFormatter(lambda value, pos: format_mF(value))
