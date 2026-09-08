# Source: ``Molecule-Structure/Source Code/formalism.py`` @ ``9eec91a``.
"""Quartic N²↔R² conversion (B&C Table 7.2); untagged dictionaries pass through."""
import warnings

# cm⁻¹ to MHz.
DEFAULT_C_CM = 29979.2458

# Verified B&C Table 7.2 generic-X pairs; this map bounds conversion.
CENTRIFUGAL_PARTNERS = {
    'p+2q':     'p2q_D',
    'Gamma_SR': 'Gamma_D',
    'q_lD':     'q_lD_D',
    'ASO':      'A_D',
}

# Sextic (H-order) keys: unsupported (quartic truncation; codebase has no H).
_SEXTIC_KEYS = {'H', 'p2q_H', 'Gamma_H', 'q_lD_H'}


def convert_formalism(params: dict, c_cm: float = DEFAULT_C_CM) -> dict:
    """Convert a copied, tagged parameter dictionary; ``Lambda`` is preserved."""
    out = dict(params)
    formalism = out.get('formalism')
    if formalism is None:
        return out
    if formalism not in ('N2', 'R2'):
        raise ValueError(
            f"Unknown 'formalism' {formalism!r}; expected 'N2' or 'R2'.")

    lam = out.get('Lambda')
    if lam is None:
        raise ValueError(
            "a set declaring 'formalism' requires 'Lambda' (|Λ|); set "
            "'Lambda': 0 for Σ. Λ is never inferred (spec §3, §6, §7).")

    sextic = _SEXTIC_KEYS.intersection(out)
    if sextic:
        raise ValueError(
            f"quartic-truncated converter; sextic keys {sorted(sextic)} "
            f"not supported.")

    # B&C Table 7.2; the origin's MHz terms require division by ``c_cm``.
    sgn = -1 if formalism == 'N2' else 1
    out['formalism'] = 'R2' if formalism == 'N2' else 'N2'

    lam2 = int(lam) ** 2
    if lam2 == 0:
        return out

    be_in = out.get('Be')
    d = out.get('D')
    if be_in is not None and d is not None:
        out['Be'] = be_in + sgn * 2 * lam2 * d

    for x, xd in CENTRIFUGAL_PARTNERS.items():
        if x in out and xd in out:
            out[x] = out[x] - sgn * lam2 * out[xd]
        elif xd in out and x not in out:
            warnings.warn(
                f"'{xd}' present without partner '{x}'; no conversion applied "
                f"(likely data error).")

    if 'Origin' in out and be_in is not None:
        d_over_c = (d / c_cm) if d is not None else 0.0
        out['Origin'] = (out['Origin']
                         - sgn * lam2 * be_in / c_cm
                         - (lam2 ** 2) * d_over_c)

    return out
