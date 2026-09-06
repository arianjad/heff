"""Generate the ThF+ isotopologue and two-photon notebook.

The notebook is generated from a script so it is reviewable as a diff and
regenerable after an API change -- same shape as notebooks/build_tutorial.py:
a CELLS list of md(...)/code(...) tuples, a build() that writes the .ipynb via
nbformat, and `if __name__ == "__main__": build()`. Nothing is imported from
build_tutorial.py; the helper patterns (write_state, dominant, select_q,
display_levels) are re-implemented here for the two-spin basis.

Run:  conda run -n heff python notebooks/build_isotopologues.py
Then: conda run -n heff jupyter execute --inplace notebooks/ThF_plus_Isotopologues.ipynb
"""
import pathlib

import nbformat as nbf

OUT = pathlib.Path(__file__).with_name("ThF_plus_Isotopologues.ipynb")


def md(text):
    return ("md", text.strip("\n"))


def code(text):
    return ("code", text.strip("\n"))


CELLS = [

    # ------------------------------------------------------------- 1 -----
    md("""
# ThF⁺ isotopologues and the two-photon operator, with `heff`

The X ³Δ₁ ground vibronic state of ²³²Th¹⁹F⁺, ²²⁹Th¹⁹F⁺ and ²²⁷Th¹⁹F⁺, J = 1–4.
²³²Th is spin-0 (no Th hyperfine), so this notebook's ²³² panels use `heff`'s
v1 one-spin basis `|J, Ω, F, m_F⟩` (`KET_C`, half-integer `F`/`m_F`, F = J ± ½
from ¹⁹F alone). ²²⁹Th (I = 5/2) and ²²⁷Th (I = 1/2) each add a second coupled
nuclear spin, giving the two-spin basis `|((J I_Th) F₁, I_F) F, m_F⟩`
(`KET_C2`, integer `F`/`m_F`).

**Conventions in force**, stamped from `pset.conventions.stamp()` — printed
below. The v1 tutorial (`ThF_plus_X3Delta1_Tutorial.ipynb`) already fixes
`n_hat`, `ef_rule`, `zeeman_sign`, `zeeman_energy`, `dg_def`, `edm_factor`,
`dipole_origin`, `formalism`. This notebook adds the four v2 conventions:
`quadrupole_convention` (B&C's q₀-is-negative-EFG sign, single-valued),
`eqq2_norm` (B&C (9.52)'s q = ±2 normalisation — `'bc_9p52_q2'`, the only
implemented value; `'petrov2018_eq23'` raises, OPEN-17), `two_photon_norm`
(`'bc_5p142_reduced'`, unit reduced one-photon elements per channel), and the
`version` tag itself, `'thf-v2'`.

**Standing caveat, read this before any ²²⁹/²²⁷ number below**: **no Th
hyperfine constant has ever been measured for any ThF⁺ isotopologue**
([TH] §2.1, gap G5 — every Th hyperfine constant used here is ab initio, an
estimate, or a Schmidt-moment placeholder). Every figure that carries a
²²⁹ or ²²⁷ number states which of those three statuses is behind it, because
the three read very differently: an *ab initio* number (`A_par_Th` for ²²⁹)
comes from an actual electronic-structure calculation of ThF⁺; an *estimate*
(`eQq0_Th`, `eQq2_Th`) is transferred from the isoelectronic ¹⁷⁷HfF⁺ ion by a
scaling argument; a *placeholder* (all of ²²⁷Th's Th constants) is a
single-particle nuclear-shell-model number standing in for a moment nobody has
measured or calculated for that nucleus.

Physics sources cited throughout: [HAM] = `docs/thf-plus-x3delta1-effective-
hamiltonian.md`, [TH] = `docs/digest-literature-th-hyperfine.md`,
[2γ] = `docs/digest-literature-two-photon.md`, [SPEC-v2] = `docs/superpowers/
specs/2026-09-05-heff-v2-isotopologues-two-photon.md`. Nothing in this
notebook asserts a physics result — `tests/` does that; here we display, and
every non-trivial number is computed by the code on screen, not transcribed.
"""),

    # ------------------------------------------------------------- 2 -----
    md("""
## 1. The basis

`|J, Ω, F₁, F, m_F⟩` with `F₁ = J + I_Th` (the inner coupling) and
`F = F₁ + I_F` (¹⁹F the outer spin) — inner-spin-first, `[SPEC-v2] §2.1, [HAM]
§9` preamble. ²³²Th (I = 0) collapses this to the v1 `F = J + I_F` basis
exactly (`KET_C`; V16 proves every `REGISTRY_C2` term equals its v1 twin at
`I_Th = 0`), so `thf_spec('232')` returns the same object as `thf_spec()`.

**Why coupled, not decoupled.** Petrov 2018's alternative is the fully
decoupled product basis `|J, m_J⟩|I_Th, m_1⟩|I_F, m_2⟩` ([SPEC-v2] §2.1). The
coupled `F₁` scheme is used here because it is what B&C's spectator theorems
((5.172), (5.174), (5.186)) are written for, and because `F₁` turns out to be
an (almost) good quantum number — [TH] §4.5 finds `ΔF₁ = ±1` mixing is a
25 kHz-level effect against 190–2670 MHz `F₁` spacings, so the coupled basis
is not just formally convenient, it tracks the real near-conserved quantity.

The dimension closed form ([SPEC-v2] §2.2): `(2I_Th+1)(2I_F+1) Σ_{J=1}^{J_max}
2(2J+1)` for the two-spin case; `(2I_F+1) Σ 2(2J+1)` for ²³² (`I_Th = 0`).
"""),

    code("""
import numpy as np

import heff
from heff import (block_by_mF, build_term_matrices, ctx_from, enumerate_kets,
                  g_factors, hamiltonian, label_lines, line_strengths,
                  sweep, terms_for_case, thf_spec, thf_v2)
from heff.elements_c2 import REGISTRY_C2, j_convergence
from heff.spec import ElecState, Spin, StateSpec
from heff.twophoton import dyad_weights, two_photon_line_strengths

import matplotlib.pyplot as plt
%matplotlib inline
plt.rcParams['figure.dpi'] = 90
np.set_printoptions(precision=6, suppress=True, linewidth=120)

ISOTOPOLOGUES = ('232', '229', '227')
J_MAX = 4

pset229 = thf_v2('229')
print('conventions.stamp():')
for key, val in pset229.conventions.stamp().items():
    print(f"  {key:<22s} {val}")


def load(iso, J_max=J_MAX, mF=None):
    \"\"\"One isotopologue's field-free basis, one signed-m_F block, the term
    matrices for it, and the Blocking of the FULL basis (so a caller can look
    up dimensions of other blocks too). '229'/'227' build on KET_C2 through
    REGISTRY_C2; '232' builds on the v1 KET_C basis through the default
    registry -- the same helper handles both dtypes (controller ruling).
    \"\"\"
    spec = thf_spec(iso, J_max=J_max)
    kets = enumerate_kets(spec)
    blocks = block_by_mF(kets)
    if mF is None:
        mF = 0.5 if iso == '232' else 0.0
    sub = kets[blocks.index[mF]]
    pset = thf_v2(iso)
    ctx = ctx_from(spec, pset)
    if iso == '232':
        tm = build_term_matrices(sub, ctx)
    else:
        tm = build_term_matrices(sub, ctx, case='c2', registry=REGISTRY_C2)
    return sub, blocks, ctx, tm


print()
for iso in ISOTOPOLOGUES:
    kets, blocks, ctx, tm = load(iso)
    I_Th = ctx.spins[0].I if ctx.spins else 0.0
    tag = 'v1 KET_C basis' if iso == '232' else 'KET_C2 basis'
    print(f"{iso}ThF+  total dim = {blocks.n_total}  I_Th = {I_Th}  ({tag})")
    print(f"  {len(blocks.labels)} signed-m_F blocks, sizes "
          f"{ {lab: len(ix) for lab, ix in sorted(blocks.index.items())} }")
"""),

    # ------------------------------------------------------------- 3 -----
    md("""
## 2. The term catalogue and its matrix elements

`terms_for_case('c2', registry=REGISTRY_C2)` for ²²⁹/²²⁷; the default
`terms_for_case('c')` for ²³². Every term below is registered with its
declared selection rules (`Rules`, checked by gate A5) and a citation string
(`term.cite`); the printout shows both, but the physics is the formula, so it
is copied out here from [HAM] with its B&C equation numbers.

**Case-(c) v1 terms** ([HAM] §2, `heff/elements_c.py`):

- **Rotation + centrifugal** (§2.1): `H_rot = B₀ J(J+1) − D₀ [J(J+1)]²`,
  diagonal.
- **Ω-doubling** (§2.3, Ng Eq. C.3): off-diagonal element between `Ω = ±1` at
  fixed `J`, `−ω_ef J(J+1)/4`; splitting `ω_ef J(J+1)/2`.
- **¹⁹F hyperfine, ΔJ = 0** (§2.4, B&C (9.50)): diagonal,
  `A∥ [F(F+1) − I(I+1) − J(J+1)] / [2J(J+1)]`.
- **¹⁹F hyperfine, ΔJ = ±1** (§2.5, B&C (9.51)): off-diagonal in `J`, same
  `A∥`, dropped by Ng but above the kHz floor.
- **Nuclear spin–rotation c_I** (§2.6, B&C (8.7)/(8.20)):
  `c_I [F(F+1) − I(I+1) − J(J+1)] / 2`.
- **Stark** (§2.7): `−n̂_sign · d_mf E_z` × the rank-1 axial geometry
  (`elements_c.dipole_geometry`), `Δm_F = 0`, parity-odd.
- **Zeeman `+G∥`** (§2.8, sign corrected from Ng Eq. C.6): `+G∥ μ_B Ω (J·n̂)(n̂·B)`.
- **Nuclear Zeeman, PT-odd EDM and scalar–pseudoscalar** (§2.8, §2.12): opt-in,
  zero by default.

**Case-(c2) v2 terms, added on top** ([HAM] §9, `heff/elements_c2.py`) — every
element below is `axial_geometry`'s master formula (the two-spectator chain,
B&C (5.172) + (5.174) twice + (5.186)) at the rank/component the physics needs,
or a v1 formula run through a field adapter:

```
<J',Om',F1',F',m'| T^k_p(molecule-frame, component q) |J,Om,F1,F,m>
  = (-1)^(F'-m') ( F'  k  F ; -m'  p  m )                                 <- (5.172)
  x (-1)^(F+F1'+k+I_F) sqrt((2F'+1)(2F+1)) { F1  F  I_F ; F' F1' k }      <- (5.174), I_F spectator
  x (-1)^(F1+J'+k+I_Th) sqrt((2F1'+1)(2F1+1)) { J F1 I_Th ; F1' J' k }    <- (5.174), I_Th spectator
  x (-1)^(J'-Om') sqrt((2J'+1)(2J+1)) ( J' k J ; -Om' q Om )              <- (5.186)
```

with `q = Om' − Om`, `Δm_F` forced to `p` by the first 3j. [HAM] §9.1: this
reduces analytically to the v1 element at `I_Th = 0` — an identity, not an
approximation.

- **Five v1 diagonal-in-`F₁` terms** (rotation, centrifugal, Ω-doubling,
  the two PT-odd terms) are the v1 formulas unchanged (§9.2, B&C (5.176):
  operators built from `J`/`Ω` alone are diagonal in `F₁`, `F`, `m_F`).
- **`hyperfine_A_par_Th`, `hyperfine_A_par_Th_dJ1`, `spin_rotation_cI_Th`**
  (§9.2 substitution table): the v1 `elements_c.hyperfine_A_par` /
  `hyperfine_A_par_dJ1` / `spin_rotation_cI` formulas with `(I, F) → (I_Th,
  F₁)` — B&C (9.50)/(9.51)/(8.20) again, same equations, `I_Th`/`F₁` in place
  of `I`/`F`.
- **`zeeman_nuclear_Th`** (§9.2.1): the one Th operator the substitution rule
  does *not* cover (a lab-frame rank-1 operator, not a scalar) — the
  two-spectator chain at `k = 1`.
- **`hyperfine_A_par_F`, `hyperfine_A_par_F_dJ1`, `spin_rotation_cI_F`**
  (§9.3): the ¹⁹F operators recoupled with `I_Th` now a spectator — `ΔF₁ =
  0, ±1` (¹⁹F is the *outer* spin, so B&C (5.176) does not apply to it).
- **`quadrupole_eQq0_Th`, `quadrupole_eQq2_Th`** (§9.4, B&C (9.52)/(9.53)):
  the Th electric quadrupole, `q = 0` diagonal and `q = ±2` mixing `Ω = ±1`.
  Structurally absent (not zero) for `I_Th < 1`, i.e. for ²²⁷Th and ²³²Th.
"""),

    code("""
print(f"{'case c (v1)':-^100}")
for t in terms_for_case('c'):
    print(f"{t.name:28s} param={str(t.param):26s} dJ={t.rules.dJ} dOm={t.rules.dOm} dF={t.rules.dF}")

print()
print(f"{'case c2 (v2)':-^100}")
for t in terms_for_case('c2', registry=REGISTRY_C2):
    dF1 = t.rules.dF1
    print(f"{t.name:28s} param={str(t.param):26s} dJ={t.rules.dJ} dOm={t.rules.dOm} "
          f"dF1={dF1} dF={t.rules.dF}")
"""),

    # ------------------------------------------------------------- 4 -----
    md("""
## 3. ²³²ThF⁺ field-free level diagram

The v1 basis, no Th hyperfine. Per J, referred to `B₀ J(J+1)`: the ¹⁹F
hyperfine (`(3/4)|A∥| = 15.07 MHz` at J = 1, [HAM] §2.4) and the Ω-doubling
(`ω_ef J(J+1)/2`, [HAM] §2.3). `A_par` is **measured** (Ng 2022 Table I);
`omega_ef` is **measured** (Ng 2022 Table I) — both carried unchanged by
`thf_v2('232')` from `thf_v1()`.

**e/f labels**: the upper Ω-doublet component has parity `(−1)^J` at every
J — equivalently, e lies `ω_ef J(J+1)/2` above f, uniformly in J (Brown 1975
convention). This is **OQ-A**, closed 2026-09-05 (`docs/open-questions.md`);
"true parity" here means the parity of the actual eigenstate — the e/f label
already has the J(J+1)-dependent rotation factored out, so "true parity
`(−1)^J`" and "e above f uniformly" are the same statement, not two different
orderings.
"""),

    code("""
kets232, blocks232, ctx232, tm232 = load('232')
pset232 = thf_v2('232')
H0_232 = hamiltonian(tm232, pset232, {'E_z': 0.0, 'B_z': 0.0})
w232, v232 = np.linalg.eigh(H0_232)
lab232 = label_lines(kets232, v232, ctx232.S, rule=pset232.conventions.ef_rule, ell=0.0, s=0.0)

OFF232 = 2.0 * pset232.value('B0')
print("232ThF+, m_F = +1/2 block, E - 2B0 (MHz):")
for s in range(len(w232)):
    l = lab232[s]
    if l['J'] > 2:
        continue
    print(f"  E[{s}] = {w232[s] - OFF232:12.5f}   J={l['J']:.0f} F={l['F']:.1f} {l['ef']}"
          f" (parity {l['parity']:+d})")

fig, ax = plt.subplots(figsize=(5, 4))
for s in range(len(w232)):
    l = lab232[s]
    if l['J'] > 2:
        continue
    e = w232[s] - OFF232 if l['J'] == 1 else w232[s] - 6 * pset232.value('B0')
    ax.hlines(e, l['J'] - 0.2, l['J'] + 0.2, color='C0' if l['ef'] == 'e' else 'C3')
ax.set_xticks([1, 2])
ax.set_xlabel('J')
ax.set_ylabel('E - B0 J(J+1)  (MHz)')
ax.set_title('232ThF+ field-free, m_F=+1/2 (A_par measured, omega_ef measured)')
plt.tight_layout()
plt.show()
print("blue = e, red = f -- upper component is e at every J (OQ-A, closed).")
"""),

    # ------------------------------------------------------------- 5 -----
    md("""
## 4. ²²⁹ and ²²⁷ThF⁺ level diagrams, side by side with ²³²

**Every ¹⁹F structure and the whole Ω-doubling now live inside one F₁
level.** The ²²⁹ J = 1 Th hyperfine spread is **~4.6 GHz — 63 % of B₀**
([TH] §4.1), so the two panels below need separate energy scales (a broken
axis would show mostly white space at this ratio): the ²²⁹/²²⁷ panel is
plotted on its own axis, referred to `B₀ J(J+1)` exactly as ²³² was, but at a
much larger vertical scale.

`A_par_Th` (²²⁹) is **ab-initio**, sign-unresolved (OPEN-16, cell 6 below).
`A_par_Th` (²²⁷) is a **Schmidt-moment placeholder** ([HAM] §9.6, OPEN-20) —
no measured or estimated μ(²²⁷Th) exists anywhere ([TH] §1.4). Both spreads
below are **computed by the code**, not the [TH] table numbers.
"""),

    code("""
def level_spread_J1(iso):
    kets, blocks, ctx, tm = load(iso)
    pset = thf_v2(iso)
    H0 = hamiltonian(tm, pset, {'E_z': 0.0, 'B_z': 0.0})
    w, v = np.linalg.eigh(H0)
    J = kets['J'][np.argmax(np.abs(v), axis=0)]
    w1 = w[J == 1]
    return kets, w, v, w1.min(), w1.max(), w1.max() - w1.min()


fig, axes = plt.subplots(1, 2, figsize=(10, 4))

_, w232b, _, lo232, hi232, spread232 = level_spread_J1('232')
OFF = 2 * pset232.value('B0')
axes[0].hlines((w232b[(np.arange(len(w232b)))] - OFF)[:4], 0.8, 1.2)
axes[0].set_title(f'232ThF+ J=1  spread = {spread232*1e3:.1f} kHz (measured A_par, omega_ef)')
axes[0].set_ylabel('E - 2B0 (MHz)')

colors = {'229': 'C1', '227': 'C2'}
for iso in ('229', '227'):
    kets, w, v, lo, hi, spread = level_spread_J1(iso)
    OFFi = 2 * thf_v2(iso).value('B0')
    J = kets['J'][np.argmax(np.abs(v), axis=0)]
    axes[1].hlines((w[J == 1] - OFFi), 0.8 if iso == '229' else 1.2,
                   1.2 if iso == '229' else 1.6, color=colors[iso], label=f'{iso}Th')
    status = thf_v2(iso).params['A_par_Th'].status
    print(f"{iso}ThF+ J=1 spread = {spread/1e3:.3f} GHz  "
          f"({spread/pset229.value('B0')*100:.1f}% of B0)   A_par_Th status={status!r}")
axes[1].legend()
axes[1].set_title('229/227ThF+ J=1 (Th hyperfine)')
axes[1].set_ylabel('E - 2B0 (MHz)')
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------- 6 -----
    md("""
## 5. The ¹⁹F doublet-ordering flip — a falsifiable prediction

Once `I_Th ≠ 0`, the ¹⁹F doublet splitting inside each `F₁` is no longer the
v1 `A∥^F (2J+1)/[2J(J+1)]` law — it depends on `F₁` and can even **flip
sign** ([TH] §4.5). Computed below directly from the `hyperfine_A_par_F` term
matrix (the diagonal projection, in units of `A∥^F`) on the default branch
`thf_v2('229')` (`A_par_Th = −1510 MHz`, the **signed, ab-initio** value
recommended by `docs/lit/lookup-apar-th-sign-convention.md`). The alternative
branch is `thf_v2('229', a_par_th_sign='positive')`, trusting Denis 2015
instead of Skripnikov & Titov 2015; **the sign of A∥(Th) is unresolved
(OPEN-16)**, and only one branch is drawn here since the ¹⁹F splitting PATTERN
below is a sign-independent geometric fact of the recoupling (only the F₁
*ordering* — which A∥_Th's sign controls — decides which physical F₁ manifold
sits where in energy).
"""),

    code("""
kets229, blocks229, ctx229, tm229 = load('229')
idxF = tm229.names.index('hyperfine_A_par_F')
mask1 = (kets229['J'] == 1) & (kets229['Om'] == 1)
rows = sorted(zip(kets229['F1'][mask1], kets229['F'][mask1], np.flatnonzero(mask1)))
diag = {(F1, F): tm229.mats[idxF][i, i] for F1, F, i in rows}

print("19F doublet splitting at J=1, in units of A_par^F (hyperfine_A_par_F term, diagonal):")
print(f"{'F1':>5}  {'E(F=F1+1/2)-E(F=F1-1/2)':>26}   (I_Th=0 reference: +0.7500)")
for F1 in sorted({F1 for F1, F in diag}):
    hi = diag[(F1, F1 + 0.5)]
    lo = diag[(F1, F1 - 0.5)]
    print(f"{F1:5.1f}  {hi - lo:26.4f}")

print()
print("A_par_Th status:", pset229.params['A_par_Th'].status,
      " value:", pset229.value('A_par_Th'), "MHz  (sign UNVERIFIED, OPEN-16)")
print("Alternative branch: thf_v2('229', a_par_th_sign='positive') trusts Denis 2015,",
      "+1510 MHz instead of the default -1510 MHz.")
"""),

    # ------------------------------------------------------------- 7 -----
    md("""
## 6. The eQq₂ versus ω_ef competition (no parameter sweep — Arian's ruling)

`eQq2_Th` mixes `Ω = +1 ↔ −1` at fixed `(J, F₁, F, m_F)`, in the same matrix
slot as `omega_doubling` ([HAM] §9.4.4, [TH] §4.4). Tabulated below at the
**default eQq2_Th = 300 MHz estimate** (`status='estimate'`, band **200–400
MHz**, OPEN-17: the bridge to Petrov 2018's normalisation carries an
unresolved `√2` and an unresolved sign) against two reference scales: the
Ω-doubling off-diagonal element `ω_ef J(J+1)/4 = 2.65 MHz` at J = 1, and the
fully-polarised Stark shift `γ_F m_F d_mf E = 50.9 MHz` at F = 3/2, m_F = 3/2,
E = 60 V/cm (the JILA operating field, [HAM] §2.7). **At the JILA field the
Stark and quadrupole scales are comparable, so the three compete rather than
one simply winning** ([TH] §4.4).
"""),

    code("""
from heff.assemble import coefficients

idx_om = tm229.names.index('omega_doubling')
idx_q2 = tm229.names.index('quadrupole_eQq2_Th')
c229 = coefficients(tm229, pset229, {'E_z': 0.0, 'B_z': 0.0})

maskJ1 = kets229['J'] == 1
mask_p1 = maskJ1 & (kets229['Om'] == 1)
mask_m1 = maskJ1 & (kets229['Om'] == -1)
rows = []
for i in np.flatnonzero(mask_p1):
    F1, F = kets229['F1'][i], kets229['F'][i]
    for j in np.flatnonzero(mask_m1):
        if kets229['F1'][j] == F1 and kets229['F'][j] == F:
            om = c229[idx_om] * tm229.mats[idx_om][i, j]
            q2 = c229[idx_q2] * tm229.mats[idx_q2][i, j]
            rows.append((F1, F, om, q2))

print(f"229ThF+, J=1, Omega=+1 <-> -1 mixing element (MHz):")
print(f"{'F1':>5} {'F':>5} {'omega_doubling':>16} {'eQq2_Th (300 MHz)':>20} {'ratio':>8}")
for F1, F, om, q2 in rows:
    print(f"{F1:5.1f} {F:5.1f} {om:16.4f} {q2:20.4f} {q2 / om:8.2f}")

d_mf = pset229.value('d_mf')
gam = lambda J, F, I=0.5: (J * (J + 1) + F * (F + 1) - I * (I + 1)) / (2 * F * (F + 1) * J * (J + 1))
stark_60 = gam(1, 1.5) * 1.5 * d_mf * 60.0

fig, ax = plt.subplots(figsize=(6, 4))
labels = [f"F1={F1:g}" for F1, F, om, q2 in rows]
ax.bar(np.arange(len(rows)) - 0.2, [abs(om) for *_, om, q2 in rows], width=0.4, label='omega_doubling')
ax.bar(np.arange(len(rows)) + 0.2, [abs(q2) for *_, om, q2 in rows], width=0.4, label='eQq2_Th (300 MHz est.)')
ax.axhline(stark_60, ls='--', c='k', label=f'Stark @ 60 V/cm = {stark_60:.1f} MHz')
ax.set_xticks(range(len(rows)))
ax.set_xticklabels(labels)
ax.set_ylabel('|Omega=+1 <-> -1 element| (MHz)')
ax.set_title('eQq2_Th [estimate, 200-400 MHz band, OPEN-17] vs omega_ef [measured] vs Stark')
ax.legend()
plt.tight_layout()
plt.show()

print(f"\\nomega_ef J(J+1)/4 at J=1 = {pset229.value('omega_ef') * 2 / 4:.3f} MHz  [omega_ef: measured]")
print(f"eQq2_Th is {min(abs(q2/om) for *_, om, q2 in rows):.1f}-{max(abs(q2/om) for *_, om, q2 in rows):.1f}x "
      "the omega-doubling element at its 300 MHz default -- 20x for two of the three F1 branches, "
      "smaller for the third.")
"""),

    # ------------------------------------------------------------- 8 -----
    md("""
## 7. J convergence

`j_convergence('229')` and `j_convergence('227')` at `J_max = 2, 4, 6`
(field-free; `m_F` resolved per-isotopologue by the helper — 0 for both,
since `F`/`m_F` are integer for both). The Th `ΔJ = ±1` hyperfine is a
**~2 GHz off-diagonal element** ([TH] §4.2) whose second-order shift is **not**
absorbable into `B₀`, because — unlike the analogous ¹⁹F term — it depends on
`F₁`. `J_max` is a `StateSpec` knob (OPEN-22): the notebook's own figures use
`J_max = 4` throughout; the table below is the evidence for that choice.
"""),

    code("""
for iso in ('229', '227'):
    report = j_convergence(iso, J_maxes=(2, 4, 6))
    print(f"{iso}ThF+  (m_F = {report['mF']}, field-free):")
    print(f"  {'J_max':>6} {'dimension':>10} {'max shift vs J_max=6 (MHz)':>28}")
    for J_max, row in report['by_J_max'].items():
        print(f"  {J_max:6d} {row['dimension']:10d} {row['max_shift_MHz']:28.6g}")
    print()

print("The J_max=1 -> 4 jump, at the full J=1 dimension (V23's own demonstration):")
for iso, mF in (('229', 0.0), ('227', 0.0)):
    kets1, blocks1, _, _ = load(iso, J_max=1)
    n1 = len(blocks1.index[mF])
    r = j_convergence(iso, J_maxes=(1, 4), mF=mF, n_levels=n1)
    print(f"  {iso}ThF+: {r['by_J_max'][1]['max_shift_MHz']/1e3:.3f} GHz "
          f"({'placeholder-driven, A_par_Th = +39.8 GHz' if iso == '227' else 'ab-initio A_par_Th'})")
"""),

    # ------------------------------------------------------------- 9 -----
    md("""
## 8. Zeeman maps and g-factors

Exact Hellmann–Feynman g-factors (`observe.g_factors`), one non-zero
signed-`m_F` block per isotopologue: `m_F = +1/2` for ²³² (matches the v1
tutorial's block), `m_F = +1` for ²²⁹/²²⁷ (smallest non-zero integer m_F).
**Predicted g at J > 1 carries an unquantified ~1 % error from the absorbed
rotational g-factor `g_r`** ([HAM] §2.10, OPEN-7) — not fitted for any
isotopologue, so the same caveat applies to every g below.
"""),

    code("""
gblocks = {}
for iso, mF in (('232', 0.5), ('229', 1.0), ('227', 1.0)):
    kets, blocks, ctx, tm = load(iso, mF=mF)
    pset = thf_v2(iso)
    gblocks[iso] = (kets, ctx, tm, pset, mF)

B = np.linspace(0.0, 5.0, 51)
fig, ax = plt.subplots(figsize=(6, 4))
for iso, (kets, ctx, tm, pset, mF) in gblocks.items():
    res0 = g_factors(tm, pset, {'E_z': 0.0, 'B_z': 0.0}, ctx=ctx)
    j1 = kets['J'] == 1
    g_j1 = res0['g'][j1]
    print(f"{iso}ThF+ (m_F={mF:g}): g at J=1 (n={j1.sum()} states) = "
          f"{np.unique(np.round(g_j1, 6))}")
    r = sweep(tm, pset, {'E_z': np.zeros_like(B), 'B_z': B})
    shift = r.evals[:, np.flatnonzero(j1)[0]] - r.evals[0, np.flatnonzero(j1)[0]]
    ax.plot(B, shift, label=f'{iso}ThF+, lowest J=1 state')
ax.set_xlabel('B_z (G)')
ax.set_ylabel('E(B) - E(0)  (MHz)')
ax.set_title('Zeeman shift, J=1, field-free m_F block per isotopologue')
ax.legend()
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------ 10 -----
    md("""
## 9. Stark maps

The three isotopologues' `m_F` blocks from §8 (`m_F = +1/2` for ²³², `+1` for
²²⁹/²²⁷), swept in `E_z`. The ²³² closed form `Ω m_F γ_F d_mf E`, evaluated
for the `J = 1, F = 1/2` component of that block, is drawn as the reference
curve on the ²³² panel — the ²²⁹/²²⁷ panels have no such closed form because
the Th hyperfine mixes many more `(F₁, F)` states at comparable energy, so
only the direct diagonalisation is shown for them.
"""),

    code("""
E = np.linspace(0.0, 60.0, 61)
fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True)
gam = lambda J, F, I=0.5: (J * (J + 1) + F * (F + 1) - I * (I + 1)) / (2 * F * (F + 1) * J * (J + 1))

for ax, iso in zip(axes, ('232', '229', '227')):
    kets, ctx, tm, pset, mF = gblocks[iso]
    j1 = np.flatnonzero(kets['J'] == 1)
    r = sweep(tm, pset, {'E_z': E, 'B_z': np.zeros_like(E)})
    for s in j1:
        ax.plot(E, r.evals[:, s] - r.evals[0, s], lw=1.1)
    ax.set_title(f'{iso}ThF+  m_F={mF:g}  ({len(j1)} J=1 states)')
    ax.set_xlabel('E_z (V/cm)')

axes[0].set_ylabel('E(E_z) - E(0)  (MHz)')
d_mf = pset232.value('d_mf')
closed = np.abs(1.0 * 0.5 * gam(1, 0.5) * d_mf * E)
axes[0].plot(E, closed, 'k--', lw=1, label='closed form |Om m_F gam_F(F=1/2) d_mf E|')
axes[0].legend(fontsize=8)
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------ 11 -----
    md("""
## 10. One-photon E1 spectra

Two adjacent `m_F` blocks per isotopologue, each diagonalised separately at
zero field, then `heff.spectra.line_strengths` (v1 geometry for ²³²,
`functools.partial(heff.elements_c2.axial_geometry, k=1, q=0.0)` for
²²⁹/²²⁷) between them, filtered by the resulting labels into **within J = 1**
(`ΔJ = 0`, same-Ω dipole transitions that only exist because both `Ω = ±1`
mix inside one `m_F` block) and **J = 1 ↔ 2**. Lines are labelled by
`(J, F₁, F, parity)` via `label_lines`.
"""),

    code("""
import functools
from heff.elements_c2 import axial_geometry

GEO = functools.partial(axial_geometry, k=1, q=0.0)


def e1_spectrum(iso, mF_a, mF_b):
    kets_a, _, ctx, tm_a = load(iso, mF=mF_a)
    kets_b, _, _, tm_b = load(iso, mF=mF_b)
    pset = thf_v2(iso)
    wa, va = np.linalg.eigh(hamiltonian(tm_a, pset, {'E_z': 0.0, 'B_z': 0.0}))
    wb, vb = np.linalg.eigh(hamiltonian(tm_b, pset, {'E_z': 0.0, 'B_z': 0.0}))
    geo = None if iso == '232' else GEO
    freqs, S = line_strengths(wa, va, kets_a, wb, vb, kets_b, ctx,
                              polarizations=(-1,), weights={-1: 1.0}, geometry=geo)
    lab_a = label_lines(kets_a, va, ctx.S, rule=pset.conventions.ef_rule, ell=0.0, s=0.0)
    lab_b = label_lines(kets_b, vb, ctx.S, rule=pset.conventions.ef_rule, ell=0.0, s=0.0)
    return freqs, S, lab_a, lab_b


PAIRS = {'232': (0.5, 1.5), '229': (0.0, 1.0), '227': (0.0, 1.0)}
fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharey='row')
for col, iso in enumerate(('232', '229', '227')):
    mFa, mFb = PAIRS[iso]
    freqs, S, lab_a, lab_b = e1_spectrum(iso, mFa, mFb)
    mask = S > 1e-6 * S.max()
    ia, ib = np.nonzero(mask)
    for row, (Jsel_a, Jsel_b, title) in enumerate(
            [(1, 1, 'within J=1'), (1, 2, 'J=1 -> 2')]):
        sel = [(i, j) for i, j in zip(ia, ib)
               if lab_a[i]['J'] == Jsel_a and lab_b[j]['J'] == Jsel_b]
        f0 = 0.0 if Jsel_a == Jsel_b else 4 * thf_v2(iso).value('B0')
        xs = [freqs[i, j] - f0 for i, j in sel]
        ys = [S[i, j] for i, j in sel]
        if xs:
            axes[row, col].stem(xs, ys, basefmt=' ')
        axes[row, col].set_title(f'{iso}ThF+ {title} ({len(sel)} lines)')
        print(f"{iso}ThF+ {title}: {len(sel)} lines above 1e-6 of the strongest")
    axes[1, col].set_xlabel('freq - ref (MHz)')
axes[0, 0].set_ylabel('strength (d_mf^2)')
axes[1, 0].set_ylabel('strength (d_mf^2)')
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------ 12 -----
    md("""
## 11. The two-photon operator

Adiabatic elimination of a far-detuned intermediate manifold gives, between
X-state levels, `T_eff = Σ_i (d·ε₂*) |i⟩⟨i| (d·ε₁) / Δ_i` ([HAM] §9.5, [2γ]
§3.1). With a common detuning, B&C Eq. (5.142) read backwards makes the
intermediate sum the reduced element of a single **rank-K polarisability**,
`α^K = T^K(d,d)/Δ` — one scalar per `(K, ΔΩ)` channel times parameter-free
geometry ([HAM] §9.5.1(1)). `heff` computes the geometry; `alphas` is the
caller's scalar per channel.

**Which channels.** Two rank-1 dipole operators couple to `K = 0, 1, 2`
(B&C (5.141)). Within X ³Δ₁, `|Ω| = 1` on both sides, so `ΔΩ ∈ {0, ±2}`, and
`ΔΩ = ±2` needs `K = 2` **and an `Ω = 0` intermediate state** — each E1 leg
carries `|q_i| ≤ 1`, and only an `Ω_i = 0` intermediate lets the two legs sum
to `ΔΩ = ±2` ([2γ] §3.3). **`K = 1` is not registered**: in *exact closure*
(a complete opposite-parity intermediate manifold) it is the antisymmetric
part of the dyad, proportional to `P_X [d_a, d_b] P_X`, which is identically
zero because the Cartesian components of the dipole commute — confirmed
numerically to `2.2 × 10⁻¹⁶` on a complete spherical-harmonic closure test
([HAM] §9.5.3, OPEN-21 resolved: `K ∈ {0, 2}` only, registered; `K = 1`
reappears only at `O(δ/Δ)` in the resolved sum or `O(1)` for a *restricted*
intermediate manifold — neither is in this operator).

**Parity.** The operator is **parity-even**: `P d P† = −d` twice. At zero
field, where `P` commutes with `H`, this means the two-photon operator
connects `e → e` and `f → f` and **never `e ↔ f`** ([2γ] §3.3). At non-zero
`E_rot` parity is not a good quantum number (the Stark term does not commute
with `P`), and the e/f restriction lifts — the regime JILA actually runs in.

**Δm_F reach.** `Δm_F ∈ {0, ±1, ±2}`, **never ±3** — the Wigner-Eckart
projection bounds `|Δm_F| ≤ K ≤ 2` ([HAM] §9.5.4). This is why the two-photon
operator is **not** the eEDM π/2 pulse (`Δm_F = ±3` within J = 1, `|m_F| =
3/2 → ∓3/2`): that transfer is done by an `E_rot` amplitude ramp, not by
2 × E1 ([2γ] §1 row 4).

**The validity condition — quoted verbatim** ([HAM] §9.5.5):

> The closure form of the two-photon operator requires a detuning large
> compared with the intermediate rotational structure, `Δ ≫ B_i ≈ 7 GHz` for
> ThF⁺ (`B_i` = the **intermediate** electronic state's rotational constant,
> `B_e ≈ 0.23 cm⁻¹`), whereas the JILA experiments run at 0.16–1.5 GHz — so
> this is the right *operator shape* and the wrong *limit for the current
> experiment*, and the α's could later be generated by a resolved sum once
> the intermediate ladder and its 0⁺/0⁻ labels are settled.

The `alpha_K*_dOm*` scalars below are **placeholders** (`status='placeholder'`,
value 1.0 in their unit): no ThF⁺ two-photon polarisability exists in any
source read ([2γ] gap 1). So the spectra in §12–13 show **geometry** — the
relative strengths the rank-K selection rules and the polarisation dyad
allow — not physical rates.
"""),

    # ------------------------------------------------------------ 13 -----
    md("""
## 12. Two-photon spectra within J = 1 and J = 1 ↔ 2

Same two `m_F` blocks as §10, now driven by `two_photon_line_strengths` at
`alphas = 1` for every registered channel (`alpha_K0_dOm0`, `alpha_K2_dOm0`,
`alpha_K2_dOm2` — all `placeholder`). Three polarisation pairs, **Raman
reading** (`dyad_weights` conjugates `ε₂`): `(σ⁺, σ⁺) → Δm_F = 0`,
`(σ⁺, σ⁻) → Δm_F = +2`, `(σ⁻, σ⁻) → Δm_F = −2`.

**Ng's opposite-helicity case, highlighted** — Ng thesis p. 102, verbatim:
"An alternative to the π-polarized microwaves is to use a two-photon Raman
process, using photons of **opposite helicities**. We would need to do some
spectroscopy to make this happen." His stated target
`|J=1,F=3/2,m_F=+3/2⟩ → |m_F=+1/2⟩` is `Δm_F = −1` and is **genuinely out of
reach** with σ± alone. `m_F = +3/2 → −1/2` is `Δm_F = −2` and **is**
reachable — [HAM] §9.5.4's last paragraph, exactly: "with a same-helicity
σ⁻σ⁻ pair in the ladder reading, or with ε₁ = σ⁻, ε₂ = σ⁺ in the Raman
reading of §9.5.1(3)" — `dyad_weights` implements the **Raman** reading, so
the `(σ⁻, σ⁺)` pair below is the one that reaches it (`[2γ] §3.3`, marked
there as **derived**, not something JILA has stated it intends).
"""),

    code("""
SQ2 = np.sqrt(2.0)
SIGMA_P = np.array([-1.0, -1.0j, 0.0]) / SQ2
SIGMA_M = np.array([1.0, -1.0j, 0.0]) / SQ2
ALPHAS = {'alpha_K0_dOm0': 1.0, 'alpha_K2_dOm0': 1.0, 'alpha_K2_dOm2': 1.0}

# label -> (eps1, eps2, Delta m_F).  Raman reading (dyad_weights conjugates
# eps2): (sig+,sig+) -> dmF=0; (sig+,sig-) -> dmF=+2; (sig-,sig+) -> dmF=-2
# ([HAM] S9.5.1(3) row (d)). two_photon_spectrum's first block argument is
# the BRA side and the second is the KET side (two_photon_matrix(kets_a,
# kets_b, ...) evaluates two_photon_geometry(kets_a[i], kets_b[j], ...),
# whose bra/ket order matches axial_geometry: Delta m_F = P = bra_mF -
# ket_mF), so the bra-side block sits at ket_mF + Delta m_F. Picking the
# wrong pair gives an all-zero spectrum, not an error, so this mapping is
# the load-bearing part.
TWOP_PAIRS = {'(sig+,sig+) dmF=0': (SIGMA_P, SIGMA_P, 0),
             '(sig+,sig-) dmF=+2': (SIGMA_P, SIGMA_M, 2),
             '(sig-,sig+) dmF=-2 [Ng p.102]': (SIGMA_M, SIGMA_P, -2)}


def _v2_shaped_232(J_max=J_MAX):
    \"\"\"232ThF+ padded into the KET_C2 dtype at I_Th=0 -- needed only here,
    because the two-photon operator lives in REGISTRY_C2/axial_geometry,
    which reads ctx.spins[0] and ket['F1'] unconditionally. V16 proves this
    basis is ket-for-ket identical to the v1 basis at I_Th=0 (tests/
    test_elements_c2_reduction.py), so this is not a second code path.
    \"\"\"
    v1spec = thf_spec('232', J_max=J_max)
    return StateSpec(case='c', electronic=v1spec.electronic, I=0.5,
                     J_range=v1spec.J_range, v=0, M='blocks', frame='rotating',
                     spins=(Spin(label='232Th', I=0.0, couple_to='J'),
                            Spin(label='19F', I=0.5, couple_to='F1')))


def two_photon_spectrum(iso, mF_a, mF_b, eps1, eps2):
    \"\"\"Diagonalise the mF_a (BRA side) and mF_b (KET side) blocks
    separately (field-free) and return the two-photon line list between
    them, geometry only (alphas = 1 for every registered channel -- the
    alphas are placeholders, so the strengths below are geometry, not rate).
    The physical Delta m_F of the transition is mF_a - mF_b.
    \"\"\"
    if iso == '232':
        spec2 = _v2_shaped_232()
        kets2 = enumerate_kets(spec2)
        blocks2 = block_by_mF(kets2)
        ctx = ctx_from(spec2, thf_v2('232'))
        kets_a, kets_b = kets2[blocks2.index[mF_a]], kets2[blocks2.index[mF_b]]
        tm_a = build_term_matrices(kets_a, ctx, case='c2', registry=REGISTRY_C2)
        tm_b = build_term_matrices(kets_b, ctx, case='c2', registry=REGISTRY_C2)
    else:
        kets_a, _, ctx, tm_a = load(iso, mF=mF_a)
        kets_b, _, _, tm_b = load(iso, mF=mF_b)
    pset = thf_v2(iso)
    wa, va = np.linalg.eigh(hamiltonian(tm_a, pset, {'E_z': 0.0, 'B_z': 0.0}))
    wb, vb = np.linalg.eigh(hamiltonian(tm_b, pset, {'E_z': 0.0, 'B_z': 0.0}))
    freqs, S = two_photon_line_strengths(wa, va, kets_a, wb, vb, kets_b, ctx,
                                         eps1=eps1, eps2=eps2, alphas=ALPHAS)
    lab_a = label_lines(kets_a, va, ctx.S, rule=pset.conventions.ef_rule, ell=0.0, s=0.0)
    lab_b = label_lines(kets_b, vb, ctx.S, rule=pset.conventions.ef_rule, ell=0.0, s=0.0)
    return freqs, S, lab_a, lab_b


mF0_229 = 0.0
fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharey='row')
for col, (label, (eps1, eps2, dmF)) in enumerate(TWOP_PAIRS.items()):
    freqs, S, lab_a, lab_b = two_photon_spectrum('229', mF0_229 + dmF, mF0_229, eps1, eps2)
    mask = S > 1e-9 * max(S.max(), 1e-30)
    ia, ib = np.nonzero(mask)
    for row, (Ja, Jb, title) in enumerate([(1, 1, 'within J=1'), (1, 2, 'J=1->2')]):
        sel = [(i, j) for i, j in zip(ia, ib) if lab_a[i]['J'] == Ja and lab_b[j]['J'] == Jb]
        f0 = 0.0 if Ja == Jb else 4 * pset229.value('B0')
        if sel:
            axes[row, col].stem([freqs[i, j] - f0 for i, j in sel],
                                [S[i, j] for i, j in sel], basefmt=' ')
        axes[row, col].set_title(f'{label}\\n{title} ({len(sel)} lines)', fontsize=8)
        print(f"229ThF+ {label} {title}: {len(sel)} lines")
    axes[1, col].set_xlabel('freq - ref (MHz)')
axes[0, 0].set_ylabel('strength (alpha^2)')
axes[1, 0].set_ylabel('strength (alpha^2)')
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------ 14 -----
    md("""
## 13. Two-photon spectra, all three isotopologues

Same `(σ⁺, σ⁻)`, `Δm_F = +2`, within-J = 1 panel, for ²³², ²²⁹ and ²²⁷.
²³² needs the I_Th = 0 `KET_C2`-shaped basis (`_v2_shaped_232` above) since
`axial_geometry`/`two_photon_geometry` require a two-spin `ctx`.
"""),

    code("""
fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
for ax, iso in zip(axes, ('232', '229', '227')):
    # 232's J=1 manifold only reaches |F| <= 1.5 (I_Th=0), so the ket-side
    # m_F has to sit at -0.5 (bra at +1.5) for both sides to have J=1 states
    # at Delta m_F=+2; 229/227 reach F=4 and F=2 respectively at J=1, so m_F=0
    # (bra at +2) already works.
    ket_mF = -0.5 if iso == '232' else 0.0
    freqs, S, lab_a, lab_b = two_photon_spectrum(iso, ket_mF + 2, ket_mF, SIGMA_P, SIGMA_M)
    mask = S > 1e-9 * max(S.max(), 1e-30)
    ia, ib = np.nonzero(mask)
    sel = [(i, j) for i, j in zip(ia, ib) if lab_a[i]['J'] == 1 and lab_b[j]['J'] == 1]
    if sel:
        ax.stem([freqs[i, j] for i, j in sel], [S[i, j] for i, j in sel], basefmt=' ')
    ax.set_title(f'{iso}ThF+  ({len(sel)} lines)')
    ax.set_xlabel('freq (MHz)')
    print(f"{iso}ThF+ within J=1, (sig+,sig-) dmF=+2: {len(sel)} lines")
axes[0].set_ylabel('strength (alpha^2, placeholder units)')
plt.tight_layout()
plt.show()
"""),

    # ------------------------------------------------------------ 15 -----
    md("""
## 14. What this model contains and does not

**In** (the v2 term list, on top of v1's nine terms): `hyperfine_A_par_Th`,
`hyperfine_A_par_Th_dJ1`, `spin_rotation_cI_Th`, `zeeman_nuclear_Th`,
`hyperfine_A_par_F`, `hyperfine_A_par_F_dJ1`, `spin_rotation_cI_F`,
`quadrupole_eQq0_Th`, `quadrupole_eQq2_Th`; the rank-K two-photon operator
(`REGISTRY_2G`, never summed into a Hamiltonian).

**Out** — sizes from [HAM] §5, unchanged by v2: `Ω = ±2, ±3` (³Δ₂, ³Δ₃);
`e_Δ` (hyperfine-dependent Ω-doubling, OPEN-8); parity-dependent Zeeman
(OPEN-10); the rotational `g_r` (OPEN-7); the rotating-frame term
(`ħω_rot F_x`, breaks m_F blocking).

**Statuses that matter, one line each**:

- `A∥(Th)` for ²²⁹Th: **ab-initio**, sign unresolved between two
  calculations that agree in magnitude — **OPEN-16**.
- `eQq₂_Th` for ²²⁹Th: **estimate**, 200–400 MHz band, inherits an
  unresolved normalisation factor between B&C and Petrov 2018 — **OPEN-17**.
- `eQq₀_Th` for ²²⁹Th: **estimate**, no published ThF⁺ or ThO quadrupole
  constant exists at all — **OPEN-18**.
- `c_I(Th)`: held at 0, unconstrained over three decades (~1 kHz–1 MHz) —
  **OPEN-19**.
- Every ²²⁷Th constant: **placeholder**, a Schmidt single-particle moment for
  a *tentative* (1/2⁺) spin assignment with no measured or estimated moment
  anywhere — **OPEN-20**.
- The two-photon `K = 1` channel: **resolved**, not registered — exact
  closure makes it identically zero — **OPEN-21** (this notebook, §11).
- `J_max`: a `StateSpec` knob, this notebook uses 4 throughout with a
  convergence table (§7) — **OPEN-22**.
- The two-photon closure form: valid at `Δ ≫ 7 GHz`, JILA runs at
  0.16–1.5 GHz — **OPEN-23**.

Where the physics is gated rather than displayed: `tests/`. `docs/open-
questions.md` carries OPEN-16 through OPEN-23 in full, each with its
citation. Design: [SPEC-v2] `docs/superpowers/specs/2026-09-05-heff-v2-
isotopologues-two-photon.md`.
"""),
]


def build():
    nb = nbf.v4.new_notebook()
    nb.cells = [nbf.v4.new_markdown_cell(src) if kind == "md"
                else nbf.v4.new_code_cell(src) for kind, src in CELLS]
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, OUT)
    n_md = sum(1 for kind, _ in CELLS if kind == "md")
    print(f"wrote {OUT} with {len(nb.cells)} cells ({n_md} markdown, {len(CELLS) - n_md} code)")


if __name__ == "__main__":
    build()
