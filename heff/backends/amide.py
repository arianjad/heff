"""C2v amides: equivalent proton pair, nitrogen, and optional metal spin.

Core operators retain C2V-Molecules' conventions; an outer metal spin is
recoupled with B&C (5.172)-(5.176), printed p.173/PDF p.205. Metal contact and
Zeeman are implemented; metal anisotropic hyperfine/quadrupole are not implied.
"""
from dataclasses import dataclass
from functools import partial
from math import sqrt

from .. import elements_amide as core
from ..amide_basis import AmideBasisSpec, enumerate_amide_kets
from ..backend_registry import Backend
from ..elements_c import _ph
from ..params import MU_B, MU_N
from ..terms import Term
from ..wigner import w3j, w6j


@dataclass(frozen=True)
class AmideConventions:
    body_axis: str = "a"
    exchange: str = "two_equivalent_protons"
    hyperfine: str = "sears_hirota"
    version: str = "amide-v1"

    def __post_init__(self):
        if (self.body_axis, self.exchange, self.hyperfine, self.version) != (
                "a", "two_equivalent_protons", "sears_hirota", "amide-v1"):
            raise ValueError("unsupported amide convention; use the documented native conventions")

    def stamp(self):
        return vars(self).copy()


@dataclass(frozen=True)
class AmideContext:
    spec: AmideBasisSpec
    conventions: AmideConventions

    @property
    def frame(self):
        return self.spec.frame


@dataclass(frozen=True)
class AmideRules:
    rank: int
    dmF: tuple = (0,)

    def allows(self, bra, ket):
        return (bra['mF'] == ket['mF'] and bra['I_T'] == ket['I_T']
                and abs(bra['F'] - ket['F']) <= self.rank)


def make_spec(basis, electronic, spins):
    allowed = {'N_min', 'N_max', 'K', 'M', 'frame'}
    if set(basis) != allowed:
        raise ValueError(f"amide basis requires exactly {sorted(allowed)}")
    if len(spins) not in (2, 3):
        raise ValueError("amide spins must contain nitrogen, an equivalent proton pair, and optionally metal")
    for record, parent in zip(spins, ('J', 'F_N', 'F_core')):
        if record['couple_to'] != parent:
            raise ValueError(f"spin {record['label']!r} must couple_to {parent!r}")
    if spins[1].get('equivalent_count') != 2 or spins[1]['I'] != 0.5:
        raise ValueError("amide backend requires two equivalent spin-1/2 protons")
    if spins[0]['I'] not in (0.5, 1.0):
        raise ValueError("nitrogen I must be 1/2 or 1 for the implemented amide operators")
    if basis['frame'] != 'lab':
        raise ValueError("amide backend implements static fields in the lab frame")
    return AmideBasisSpec(
        S=electronic['S'], I_N=spins[0]['I'], i_H=spins[1]['I'],
        I_M=spins[2]['I'] if len(spins) == 3 else 0.0,
        N_range=(basis['N_min'], basis['N_max']), K_values=tuple(basis['K']),
        vibronic_sign=electronic['vibronic_sign'], M=basis['M'], frame=basis['frame'])


_CORE_KEYS = ('N', 'K', 'J', 'F_N', 'I_T', 'F_core')
_COEFFICIENTS = ('A', 'B', 'C', 'eps_xx', 'eps_yy', 'eps_zz', 'a_H', 'TH_aa',
                 'TH_bb', 'a_N', 'TN_aa', 'TN_bb', 'Q_aa', 'Q_bb', 'g_l')


def _args(bra, ket, mb, mk):
    return tuple(value for key in _CORE_KEYS for value in (bra[key], ket[key])) + (mb, mk)


def _core_value(fn, symbol, bra, ket, ctx, mb, mk):
    args = _args(bra, ket, mb, mk)
    if symbol is not None:
        if symbol in ('Q_aa', 'Q_bb') and ctx.spec.I_N < 1:
            return 0.0
        params = dict.fromkeys(_COEFFICIENTS, 0.0)
        params.update(S=ctx.spec.S, I_N=ctx.spec.I_N)
        if symbol.startswith('g_l_'):
            params['g_l'] = {k: {q: 0.0 for q in range(-k, k+1)} for k in range(3)}
            k, q = map(int, symbol[-2:])
            params['g_l'][k][q] = 1.0
            params['g_l'][k][-q] = 1.0
        else:
            params[symbol] = 1.0
        return fn(*args, params)
    if fn is core.hydrogen_zeeman_hamiltonian_element:
        return fn(*args)
    if fn is core.nitrogen_zeeman_hamiltonian_element:
        return fn(*args, ctx.spec.I_N)
    return fn(*args, ctx.spec.S, ctx.spec.I_N)


def _core_reduced(fn, symbol, bra, ket, ctx):
    """Recover a rank-1 core RME from a nonzero angular factor, B&C (5.172)."""
    fb, fk = float(bra['F_core']), float(ket['F_core'])
    if abs(fb - fk) > 1 or fb + fk < 1:
        return 0.0
    m = min(fb, fk)
    angular = _ph(fb - m) * w3j(fb, 1, fk, -m, 0, m)
    return _core_value(fn, symbol, bra, ket, ctx, m, m) / angular


def _lift(fn, symbol, rank, factor, bra, ket, ctx):
    if bra['mF'] != ket['mF']:
        return 0.0
    if ctx.spec.I_M == 0:
        return factor * _core_value(fn, symbol, bra, ket, ctx, bra['mF'], ket['mF'])
    if rank == 0:
        if bra['F'] != ket['F']:
            return 0.0
        return factor * _core_value(fn, symbol, bra, ket, ctx,
                                     bra['F_core'], ket['F_core'])
    gb, gk, fb, fk = (float(bra['F_core']), float(ket['F_core']),
                      float(bra['F']), float(ket['F']))
    m, I = float(ket['mF']), ctx.spec.I_M
    reduced = _core_reduced(fn, symbol, bra, ket, ctx)
    return (factor * _ph(fb-m) * w3j(fb, 1, fk, -m, 0, m)
            * _ph(fk+gb+1+I) * sqrt((2*fb+1)*(2*fk+1))
            * w6j(gk, fk, I, fb, gb, 1) * reduced)


def _metal_contact(bra, ket, ctx):
    if bra['F'] != ket['F'] or bra['mF'] != ket['mF'] or ctx.spec.I_M == 0:
        return 0.0
    I, f = ctx.spec.I_M, float(ket['F'])
    gb, gk = float(bra['F_core']), float(ket['F_core'])
    reduced = _core_reduced(core.spin_iso_zeeman_hamiltonian_element, None, bra, ket, ctx)
    return (_ph(gk+f+I) * w6j(I, gk, f, gb, I, 1) * reduced
            * sqrt(I*(I+1)*(2*I+1)))


def _metal_zeeman(bra, ket, ctx):
    if any(bra[key] != ket[key] for key in _CORE_KEYS) or bra['mF'] != ket['mF']:
        return 0.0
    I, g = ctx.spec.I_M, float(ket['F_core'])
    fb, fk, m = float(bra['F']), float(ket['F']), float(ket['mF'])
    return (-MU_N * _ph(fb-m) * w3j(fb, 1, fk, -m, 0, m)
            * _ph(fb+g+I+1) * sqrt((2*fb+1)*(2*fk+1)*I*(I+1)*(2*I+1))
            * w6j(I, fk, g, fb, I, 1))


def _registry():
    terms, units = {}, {}
    source = "C2V-Molecules atm_core/physics.py @ 3b77b021; B&C (5.172)-(5.176), printed p.173/PDF p.205"
    def add(name, params, rank, fn):
        terms[name] = Term(name, params, ('amide_c2v',), AmideRules(rank), True, True, source, fn)
    for fn, fields in (
        (core.rotational_hamiltonian_term, [('rotation_A','A'), ('rotation_B','B'), ('rotation_C','C')]),
        (core.spin_rotation, [(f'spin_rotation_{s}', f'eps_{s}') for s in ('xx','yy','zz')]),
        (core.hydrogen_hyperfine, [('hyperfine_H_contact','a_H'), ('hyperfine_H_aa','TH_aa'), ('hyperfine_H_bb','TH_bb')]),
        (core.nitrogen_hyperfine, [('hyperfine_N_contact','a_N'), ('hyperfine_N_aa','TN_aa'), ('hyperfine_N_bb','TN_bb')]),
        (core.quadrupole, [('quadrupole_N_aa','Q_aa'), ('quadrupole_N_bb','Q_bb')]),
    ):
        for name, symbol in fields:
            add(name, (symbol,), 0, partial(_lift, fn, symbol, 0, 1.0))
            units[symbol] = 'MHz'
    for name, symbol, knob, fn, field_symbol, factor in (
        ('stark_z','d_0','E_z',core.stark_hamiltonian_element,None,-1),
        ('zeeman_electron','g_s','B_z',core.spin_iso_zeeman_hamiltonian_element,None,MU_B),
        *((f'zeeman_anisotropic_{c}', f'g_l_{c}', 'B_z',
           core.spin_aniso_zeeman_hamiltonian_element, f'g_l_{c}', -MU_B)
          for c in ('00', '20', '22')),
        ('zeeman_protons','g_H','B_z',core.hydrogen_zeeman_hamiltonian_element,None,-MU_N),
        ('zeeman_nitrogen','g_N','B_z',core.nitrogen_zeeman_hamiltonian_element,None,-MU_N),
    ):
        add(name, (symbol,knob), 1, partial(_lift, fn, field_symbol, 1, factor))
        units[symbol] = 'MHz/(V/cm)' if symbol == 'd_0' else ''
    add('hyperfine_M_contact', ('a_M',), 0, _metal_contact)
    add('zeeman_metal', ('g_M','B_z'), 1, _metal_zeeman)
    units.update(a_M='MHz', g_M='')
    return terms, units


REGISTRY, UNITS = _registry()
AMIDE_BACKEND = Backend(
    id='amide_c2v', case='amide_c2v', registry=REGISTRY,
    default_terms=tuple(REGISTRY), canonical_units=UNITS,
    runtime_knobs=frozenset({'E_z','B_z'}), optional_zero_parameters=frozenset(),
    make_spec=make_spec, default_conventions=AmideConventions,
    enumerate_kets=enumerate_amide_kets,
    make_context=lambda spec, params: AmideContext(spec, params.conventions))
