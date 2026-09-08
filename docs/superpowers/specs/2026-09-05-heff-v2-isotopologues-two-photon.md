# Reference design — Th isotopologues and the effective two-photon operator

This document records the two-nuclear-spin and effective two-photon contracts
adopted for `heff` on 2026-09-05. It extends the
[v1 reference](2026-09-05-heff-design.md): term-matrix catalogue
`H = Σ_k c_k M_k`,
declared selection rules as data, a conventions block stamped on every result, blocking
as a first-class object, typed `Param` with unit/status/source, no hash, snapshot or
pinned-spectrum gates, every gate naming the failure mode it uniquely catches with both
PASS and FAIL reachable, Hamiltonian-agnostic kernel gates.

The scientific sources are the [Hamiltonian reference](../../thf-plus-x3delta1-effective-hamiltonian.md)
(**[HAM]**), [Th hyperfine digest](../../digest-literature-th-hyperfine.md)
(**[TH]**), [two-photon digest](../../digest-literature-two-photon.md)
(**[2γ]**), and [²²⁷Th nuclear-moment lookup](../../lit/lookup-227th-nuclear-moment.md).

---

## 1. Problem restatement, goals, non-goals

**Problem.** The original model computes ²³²Th¹⁹F⁺ X ³Δ₁ with one nuclear spin.
The isotope and two-photon extension requires two structures it cannot express:

1. **²²⁹ThF⁺ (I_Th = 5/2) and ²²⁷ThF⁺ (I_Th = ½, tentative)** — a second nuclear spin,
   a Th magnetic hyperfine constant ≈ 1.5 GHz (21 % of B₀), and, for I_Th = 5/2, an
   electric quadrupole with a ΔΩ = ±2 component estimated at 13–30× the Ω-doubling
   operator ([TH] §4.4). The v1 ket dtype has no place to put F₁ and the v1 element set
   has no rank-2 nuclear operator and no ΔΩ = ±2 operator other than `omega_doubling`.
2. **An effective two-photon (2 × E1) operator within X**, distinct from an
   electric-quadrupole transition. The original `spectra.py` is one-photon only and
   `elements_c.dipole_geometry` hard-codes the molecule-frame component `q = 0`.

**Goal — the v2 slice.** ²³²/²²⁹/²²⁷Th¹⁹F⁺ X ³Δ₁, v = 0, in the coupled case-(c) basis
`|J, Ω, F₁, F, m_F⟩` with `F₁ = J + I_Th`, `F = F₁ + I_F`, blocked by signed m_F; a Th
hyperfine, quadrupole and spin–rotation term set; every v1 term re-coupled through the
extra spectator spin; a parameter-free rank-K two-photon operator within X; and a new
notebook that narrates the structure and the one- and two-photon spectra for the three
isotopologues side by side.

**The master requirement, and it drives every decision below.** With **I_Th = 0** the v2
basis must reduce **exactly** to the v1 basis, and every v2 matrix element must equal its
v1 counterpart to machine precision. This is the one gate that can catch a wrong 6j column
order or a wrong spectator phase in any of the seven re-coupled terms, and both outcomes
are reachable (a deliberately transposed 6j fails it; the shipped code passes it).

**Non-goals.** No new coupling case (still case (c), one electronic state — T2/T3 of
[HAM] §1.4 stay out). No resolved-intermediate two-photon spectrum as a *product*
(§3 keeps it as a *test fixture* only, because the ThF⁺ intermediate ladder is contested,
[2γ] §2.2, gap 3). No MQM operator (it needs I ≥ 1 and W_M, and the observable is not in
this slice; [HAM] §2.12). No fitting. No parity-dependent Zeeman, no `e_Δ` — those are
[HAM] OPEN-8 and OPEN-10 and are unchanged by isotopologues.

---

## 2. Architecture decision A — the two-nuclear-spin basis

### 2.1 The choice, and the two alternatives

| | **A1 coupled `|((J I_Th)F₁, I_F)F, m_F⟩`** | A2 decoupled `|J Ω m_J⟩|I_Th m₁⟩|I_F m₂⟩` | A3 hybrid: couple I_F to J, I_Th free |
|---|---|---|---|
| dimension | (2I_Th+1)(2I_F+1) Σ_J 2(2J+1) | identical | identical |
| Th hyperfine | the *existing* v1 closed forms with `F → F₁`, by B&C (5.176) ([TH] §3.2) | Clebsch sums, no 6j, but no closed form to gate against | worst of both |
| ¹⁹F hyperfine | one new recoupler (B&C 5.173 + 5.174), ΔF₁ = 0, ±1 ([TH] §3.3) | Clebsch sums | closed form |
| reduces to v1 at I_Th = 0 | **yes, ket-for-ket and element-for-element** | **no** — a different basis entirely; every v1 gate would have to be rewritten | partially |
| labels | F₁ is a near-exact quantum number: ΔF₁ = ±1 ¹⁹F elements are 1.3–3.8 MHz against F₁ spacings of 190–2670 MHz ⇒ second-order shifts ≲ 25 kHz ([TH] §4.5) | m_F only | — |
| precedent | Petrov 2018 *labels* every ¹⁷⁷/¹⁷⁹HfF⁺ result by F₁ = J + I_Hf ([TH] §2.3) | Petrov 2018 *computes* in this basis (their Eq. 15) | none |

**Decision: A1, the fully coupled basis, ordered inner-spin-first.**

The deciding argument is not elegance, it is the master gate. A2 is a perfectly good basis
— it is what Petrov actually diagonalises — but it cannot reduce to v1, so adopting it
would throw away 164 passing tests and every closed form in [HAM] §6 as a check on the new
code. A1 keeps them all: at I_Th = 0 the F₁ column is identically J, the recouplers
collapse to identity 6j's, and `hyperfine_A_par`, `stark_z`, `zeeman_Gpar`,
`zeeman_nuclear`, `omega_doubling`, `spin_rotation_cI` and both PT-odd terms must come back
bit-for-bit. A2's own advantage — no 6j's — is worth nothing here, because the term-matrix
architecture pays the Wigner cost **once per basis**, not once per parameter set (v1 spec
§3.3); at dim 576 that is seconds.

Petrov's choice is not evidence against A1: he needs a decoupled basis because he carries
four electronic states and off-diagonal electronic operators, and he still reports
everything by F₁ because that is the physical label.

*What could go wrong.* (i) F₁ stops being a good label if `eQq₂` or a strong field mixes
F₁ manifolds — it does not, at the estimated sizes ([TH] §4.5), but the notebook must
report the F₁ purity of each eigenvector rather than assert the label. (ii) Coupling order
matters: `((J I_Th) F₁, I_F) F` and `((J I_F) F', I_Th) F` are different bases related by a
9j, and the spectator phases differ between them (v1 spec §2.2, thesis item 14). The order
is fixed once, in the enumerator, and stamped in the manifest.

### 2.2 Ket dtype: one new field, not a general list of coupled spins

**Decision: a second dtype `KET_C2 = [('J','f8'),('Om','f8'),('F1','f8'),('F','f8'),('mF','f8')]`,
with `KET_C` left byte-identical.**

Three expressible options:

- **(a) add `F1` to `KET_C`.** Rejected: it changes the v1 dtype, and
  `tests/test_spec.py::test_ket_dtype_is_named_float64_fields` asserts
  `KET_C.names == ("J","Om","F","mF")`. A v2 that starts by editing a v1 gate has given up
  the master gate's independence.
- **(b) a general list of coupled spins — an object ket, or a variable-width dtype.**
  Rejected on the v1 spec's own grounds (§3.1): a structured float64 array is what keeps
  selection-rule masks vectorised and the cache key canonical, and a variable-width dtype
  makes both hard. It also buys nothing until a *third* spin exists, and no ThF⁺ isotopologue
  has one. This is the YAGNI line: pay one field now, pay a redesign if a polyatomic with
  three quadrupolar nuclei ever arrives (named trigger).
- **(c) two dtypes, one enumerator that dispatches. — chosen.** `enumerate_kets(spec)`
  returns `KET_C` when the spec carries one spin and `KET_C2` when it carries two. Element
  functions read named fields, so the dtype is invisible to them except where F₁ is
  genuinely needed.

**The spin chain is declared on the spec, not inferred.** `Spin(label, I, couple_to)` —
the frozen dataclass the v1 spec §3.1 already sketched — with
`spins=(Spin('229Th', 2.5, 'J'), Spin('19F', 0.5, 'F1'))`. `couple_to` is data, and the
enumerator is the only code that reads it. Ordering is inner-first, and a spec whose
`couple_to` chain is not a chain raises (structure, so raising is correct).

**`thf_spec()` stays backward-compatible by staying the same function with the same
default.** `thf_spec()` with no arguments returns today's `StateSpec` object unchanged
(`I=0.5`, `spins=()`), so `enumerate_kets(thf_spec())` still returns the same 96 `KET_C`
kets in the same order. `thf_spec('229')`, `thf_spec('227')`, `thf_spec('232')` return the
isotopologue specs; `thf_spec('232')` is required to compare equal to `thf_spec()`.
`J_max` is a keyword (`thf_spec('229', J_max=6)`), default 4, because the Th ΔJ = ±1
hyperfine is a ~2 GHz off-diagonal element ([TH] §4.2) and the truncation is a physics
choice, not a constant.

**Dimension is a derived closed form, and therefore a gate:**
`dim = (2I_Th+1)(2I_F+1) Σ_{J=J_min}^{J_max} 2(2J+1)`, because the F₁/F recoupling is a
change of basis inside fixed (J, Ω, m_total). For J = 1–4: 48 × 2 = 96 (v1 ✓), × 6 = **576**
for ²²⁹ThF⁺, × 2 = **192** for ²²⁷ThF⁺ `[derived]`. Per-signed-m_F block dimensions are
**computed by the enumerator and never asserted from a document**, exactly as in v1.

### 2.3 Blocking and the `Blocking`/`merge` API survive untouched

`block_by_mF` reads only `kets["mF"]`; `Blocking` carries `{label: index_array}` and
`merge`; `blocks_for` switches on `spec.M`. None of them mentions F, F₁ or a dtype. **No
change.** The `build_term_matrices` guard that refuses to build a Δm_F ≠ 0 term on a
single-m_F block is likewise unchanged and becomes load-bearing for the first time in v2:
the two-photon operator declares `dmF=(0,±1,±2)` and must be built on merged blocks.

*What could go wrong.* A larger basis makes the per-block build cost visible: 576 states over
~21 signed-m_F blocks with ~15 terms is roughly 10× v1's Wigner count. The mitigation is the
one already in the architecture — build once, resum forever — plus the declared-rules
sparsity mask that `build_term_matrices` already applies. If a build exceeds ~30 s the plan's
pre-flight number is what makes that visible; nothing here gates on it.

### 2.4 The element layer: one generic kernel, four consumers

This is where v2 either stays small or doubles the element file. The lazy reading of [TH]
§3.2–3.3 is that **there is one geometric object underneath five of the new terms**: a
molecule-frame rank-k operator with a definite axial component q, threaded out through
however many spectator spins the basis has.

```
axial_geometry(bra, ket, ctx, *, k, q, p)
  = [Wigner–Eckart 3j on F, lab component p]            B&C (5.172)
  × [spectator reduction I_F out of F]                  B&C (5.174), j1=F1, j2=I_F
  × [spectator reduction I_Th out of F1]                B&C (5.174), j1=J,  j2=I_Th
  × [rotation-matrix reduced element, (J' k J; -Ω' q Ω)] B&C (5.186)
```

with `q = Ω_bra − Ω_ket` read off the kets, guarded by `|q| ≤ k`. Consumers:

| consumer | k | q | p |
|---|---|---|---|
| Stark ([HAM] §2.7, Ng Eq. C.5) | 1 | 0 | 0 (and ±1 for E1 spectra) |
| Zeeman G∥ ([HAM] §2.8) | 1 | 0 | 0 |
| E1 transition dipole (`spectra`) | 1 | 0 | 0, ±1 |
| two-photon rank-K (§3) | K ∈ {0,1,2} | 0, ±2 | −K…K |

and, at one spin and `k=1, q=0`, it must reproduce `elements_c.dipole_geometry` exactly —
which is the second non-vacuous gate on the kernel, independent of the master gate.

The **axial hyperfine** is *not* this object; it is the scalar product of a rank-1
molecule-frame axial operator (the `m3` factor above at k=1, q=0, **not** diagonal in J)
with `T¹(I_x)`. So:

- **Th (inner spin).** By B&C (5.176) — a scalar built entirely from operators acting
  inside F₁ is diagonal in F, m_F and independent of them ([TH] §3.2) — the matrix element
  is *literally* the v1 formula with `I → I_Th, F → F₁`. **Reuse the v1 functions through a
  field adapter; write no new algebra.** Same for `spin_rotation_cI` (B&C 8.7/8.20) and for
  the ΔJ = ±1 form (B&C 9.51).
- **¹⁹F (outer spin).** Needs the genuine recoupler: B&C (5.173) across the F₁/I_F pair,
  then B&C (5.174) with I_Th as the spectator ([TH] §3.3). ΔF₁ = 0, ±1. This kernel also
  supplies the ΔJ = ±1 ¹⁹F term, which [TH] §3.3 explicitly leaves "derivable-not-derived"
  — it cannot be deferred, because at I_Th = 0 it must reproduce
  `hyperfine_A_par_dJ1`.

The **quadrupole** gets its own kernel from B&C (9.52)/(9.53) with `I → I_Th, F → F₁`,
because it is a scalar product of two rank-2 tensors with a `(I 2 I; −I 0 I)^{-1}` nuclear
normalisation that no other term shares. Its q = 0 and q = ±2 components are two registered
terms with different declared rules, sharing one body.

**Registry: a second registry, not new machinery.** `terms.term()` already takes
`registry=REGISTRY` as a keyword, so `heff/elements_c2.py` declares
`REGISTRY_C2 = {}` and registers `cases=('c2',)` terms into it;
`build_term_matrices(kets, ctx, case='c2', registry=REGISTRY_C2)` already accepts both.
**The v1 `REGISTRY` is not touched**, so every v1 test that iterates it keeps seeing exactly
eleven terms. `Rules` gains one optional field `dF1: tuple | None = None`, ignored when
`None` or when the ket dtype has no `F1` — backward-compatible, and it is what lets gate A5
run unchanged on the v2 terms.

**Two small generalisations to existing modules, each one line of intent:**

- `conventions.parity_operator` keys its ±Ω pairing on the literal four v1 field names. It
  becomes "key on every dtype field, with `Om` negated" — dtype-agnostic, identical output
  on `KET_C`.
- `spectra.dipole_matrix` hard-calls `dipole_geometry(…, ctx.I, p)`. It gains a
  `geometry=dipole_geometry` keyword; v2 passes the recoupled kernel. No duplication of
  `line_strengths`, `_strengths_from_matrices` or `label_lines`.
- `terms.Ctx` gains `spins: tuple = ()` (default empty ⇒ v1 behaviour), and `ctx_from`
  fills it from `spec.spins`. `Ctx.I` keeps meaning "the ¹⁹F spin" so v1 elements are
  untouched.

---

## 3. Architecture decision B — the effective two-photon operator

### 3.1 The choice, and why not the resolved sum

[2γ] §4.2 lays out two shapes and declines to pick. **Arian picked B: a rank-K
polarisability within X** ([inventory] §5 plus his instruction for this design). The spec
records the trade, because the notebook narrative has to state it honestly.

| | **B — rank-K effective operator within X** (chosen) | A — resolved sum over an explicit intermediate manifold |
|---|---|---|
| fits `H = Σ c_k M_k` | yes: parameter-free geometry × one coefficient per channel, so it gets the sparsity mask, the A5 gate, the citation field | no: the amplitude depends on Δ non-linearly ([2γ] §4.2) |
| new inputs required | the α^K_{ΔΩ}, for which **no ThF⁺ value exists** | the ThF⁺ intermediate ladder, which is **contested**: Denis puts Ω = 0 at 6344–6747 cm⁻¹, Petrov puts ³Π₀∓ at 3044/3395 cm⁻¹, nothing is measured between 3150 and 10 472 cm⁻¹ ([2γ] §2.2, gap 3); no 0⁺/0⁻ labels (gap 2); no excited-state hyperfine (gap 5); no lifetimes (gap 4) |
| validity | the unresolved/closure limit, Δ ≫ intermediate rotational spacing ≈ **2B ≈ 7 GHz** ([2γ] §3.4, from B ≈ 0.23 cm⁻¹) — *not* the JILA 0.16–1.5 GHz regime | exact at any detuning |
| what it hides | intermediate-symmetry interference in the ΔΩ = 0 Q branch, which Bray & Hochstrasser identify as the physically interesting content ([2γ] §3.0) | nothing |

**The honest statement, which the notebook must carry:** B is the *shape* the package
should have, and A is the *regime the experiment runs in*. [2γ] §4.2's closing line is the
reconciliation — "A can generate B's parameters, B cannot generate A's spectra". So v2
ships B as the operator, ships A **only as a test fixture** (§3.3), and the notebook states
the validity condition (Δ ≫ ~7 GHz) and that the α's could later be generated by a resolved
sum once the ladder is settled.

### 3.2 Structure, and where it lives

The operator is `T_eff = Σ_K Σ_P (−1)^P (ε₁⊗ε₂)^K_{−P} α^K_P`, with the molecular rank-K
tensor built from two rank-1 dipoles on the same (electronic–rotational) part via
**B&C (5.141)** (the coupled tensor) and **B&C (5.142)** (its reduced element, whose
`Σ_{η″ j″}` sum *is* the intermediate sum and whose 6j `{1 1 K; j′ j j″}` is what closure
collapses) ([2γ] §3.2). The spectator theorem over the nuclear spins applies to the composite
rank-K operator exactly as it does to the one-photon dipole, replacing `{J F I; F′ J′ 1}` by
`{J F I; F′ J′ K}` — and, in v2's two-spin basis, by the same two nested B&C (5.174)
reductions as §2.4's kernel. So **the two-photon geometry is `axial_geometry` at k = K**,
not a new kernel.

**Channels.** Each E1 leg carries `|q_i| ≤ 1`, so ΔΩ = q₁ + q₂ ∈ {0, ±2} within X (±1 is
excluded because |Ω| = 1 on both ends), and Ω = +1 → −1 requires an Ω = 0 intermediate
([2γ] §1, §3.3, verified numerically there against an Ω = 1 control that returns the other
answer). In the rank-K form the molecule-frame component of `α^K` is `q = Ω_bra − Ω_ket`,
and a rank-K spherical tensor has components only for `|q| ≤ K`. Therefore:

> **The ΔΩ = ±2 channel exists only at K = 2.** `[derived from |q| ≤ K; to be confirmed
> numerically in the derivation task]`

and the ΔΩ = 0 channel exists at K = 0, 1, 2. **K = 1 is a live question, not a channel to
assume:** the antisymmetric part of `d ⊗ d` is `∝ d × d`, which vanishes for a vector
operator whose Cartesian components commute — i.e. K = 1 is identically zero in exact
closure with a single common denominator, and survives only through the detuning asymmetry
between the two orderings. `[inference from B&C 5.141 with k₁ = k₂ = 1 and A₁ = B₁ = d;
NOT verified — this is a deliverable of the derivation task, and the plan does not code a
K = 1 channel until it is settled.]`

**Selection rules that follow, and are declarable as data:** |ΔF| ≤ K, Δm_F = P with
|P| ≤ K, so Δm_F ∈ {0, ±1, ±2} and **never ±3** ([2γ] §3.2) — which is why the Δm_F = 3
Ramsey π/2 pulse is done by an E_rot ramp and not by 2 × E1 ([2γ] §1 row 4), a sentence the
notebook owes the reader. Under the JILA constraint that only σ± exist in the E_rot plane
(Ng thesis p. 102 fn. 4), the reach narrows to Δm_F ∈ {0, ±2} ([2γ] §3.3).

**Parity: the operator is parity-EVEN**, so at zero field it connects e→e and f→f and never
e↔f — proved and numerically verified in [2γ] §0 item 6 and §3.3
(`max|[T,P]| ≤ 5.6 × 10⁻¹⁷`). At non-zero E_rot parity is not good (v1 gate V8) and the
restriction lifts. The notebook must say this instead of advertising a "parity flip".

**Where it lives: `heff/twophoton.py`, registered in a third registry, built like
`spectra`, not like `assemble`.** Three options were weighed:

- **in `elements_c2.py`, `@term`-registered into the Hamiltonian registry** — rejected:
  a transition operator is not a Hamiltonian term, and `build_term_matrices` would sum it
  into every assembled `H`. That is a physics bug waiting behind an architectural convenience.
- **entirely inside `spectra.py` as bare functions** — rejected: it loses `Rules`, the
  citation field, and gate A5, which is exactly the metadata that made v1's dead-operator
  class of bug automatable.
- **its own module and its own registry `REGISTRY_2G`, matrices built by a
  `spectra`-shaped function — chosen.** The `@term` decorator supplies the declared rules,
  the citation and the A5 gate for free (`term(..., registry=REGISTRY_2G)` already exists);
  `two_photon_matrix(kets_a, kets_b, ctx, *, K, dOmega, P)` supplies the matrix between two
  blocks, mirroring `spectra.dipole_matrix`; and `two_photon_line_strengths(...)` mirrors
  `spectra.line_strengths` — **summing amplitudes over (K, P) and then squaring**, which is
  v1's gate-B8 convention and, for two photons, is the convention Cossel's measured
  destructive interference between two σ pathways demands ([2γ] §3.0, his Fig. 6.18).

The polarisation enters only as the coupled dyad weights `(ε₁⊗ε₂)^K_P`, computed from two
Jones vectors by one Clebsch–Gordan contraction — a small user-facing API (`eps1`, `eps2`),
and the one place where Ng's opposite-helicity proposal (σ⁺σ⁻) and Cossel's measured
`W ∝ sin²θ` linear-polarisation dependence both live.

### 3.3 The closure test — the one test that catches a wrong 6j or phase in B&C (5.142)

Everything above rests on one algebraic step: that the `Σ_{η″ j″}` in B&C (5.142) collapses
onto a rank-K reduced element when the denominator is common. That step is a 6j and a
phase, and a wrong one is invisible in every symmetry gate (parity, sum rules and Δm_F reach
all survive a wrong `{1 1 K; j′ j j″}`).

**So the plan builds the resolved sum as a fixture and compares.** Construct a *fake*
intermediate manifold — a synthetic Ω_i ∈ {0, ±1} state at an arbitrary term energy with
arbitrary reduced dipoles — sum the resolved two-leg amplitudes over
(Ω_i, J′, F₁′, F′, m′_F, q₁, q₂) with a **common** detuning, and require it to equal the
rank-K form contracted with the same polarisations, to machine precision, for every
(initial, final) pair. The fixture is fake on purpose: the ThF⁺ ladder is contested
([2γ] §2.2), so a *physical* manifold would make the test depend on an unsettled input,
while the algebra it checks does not. **FAIL is reachable in one line** — transpose two
columns of the `{1 1 K; j′ j j″}` 6j and the comparison breaks while every other two-photon
gate still passes.

---

## 4. Conventions block additions

Three fields, each a real fork or a real trap, each stamped on every result. Fields are added
to `Conventions` and to `_ALLOWED`, so `stamp()` gains three keys; no test asserts the stamp
as a whole (measured, §0), and `version` moves to `'thf-v2'` only when a v2 spec is in play.

| field | default | alternatives | why it exists |
|---|---|---|---|
| `quadrupole_convention` | `'bc_q0_is_negative_efg'` | (one value for now) | B&C state verbatim, on the page that carries (9.53), that "q₀ is the negative of the electric field gradient". [TH] §3.4 verified this numerically: B&C (9.53) at Ω = 0 reproduces the textbook Casimir function with a **uniform ratio of exactly −1** across (J, I) = (1,1), (2,1), (2,3/2), (3,5/2) and every F. Code it explicitly or the sign of every quadrupole splitting flips. Single-valued like v1's `zeeman_energy`: it is a trap to record, not a fork to choose. |
| `eqq2_norm` | `'bc_9p52_q2'` | `'petrov2018_eq23'` | The **UNVERIFIED normalisation bridge** ([TH] §3.2). B&C define `eq₀Q` by the q = 0 specialisation of (9.52); Petrov 2018 define `eQq₀`/`eQq₂` by their Eqs. (22)–(23) with a √6 and a Y₂₂ that have not been shown to map onto B&C's T²_{±2} with the same ¼ prefactor. The package computes in B&C's normalisation; the converter exists so a *published* Petrov-style number can be entered without hand-editing a factor, in exactly the pattern `n_hat` already uses. The derivation task resolves the factor and states it with citation **before** any eQq₂ is coded. **OPEN-17.** |
| `two_photon_norm` | `'bc_5p142_reduced'` | — | Pins what `α^K_{ΔΩ}` multiplies: the geometry returned by `two_photon_matrix` is dimensionless and normalised to B&C (5.142)'s reduced element with unit one-photon reduced elements, so a strength comes out in units of α². Recorded, not chosen — but recorded, because an α transplanted from a Placzek-convention source would otherwise be silently off by a rank-dependent factor. |

Unchanged: `n_hat`, `ef_rule`, `zeeman_sign`, `zeeman_energy`, `dg_def`, `edm_factor`,
`dipole_origin`, `formalism`. The e/f rule stays `brown1975`; OQ-A is closed (the upper
Ω-doublet component has true parity (−1)^J at every J, [inventory] §5) and the code is
unchanged.

---

## 5. Parameter table

`params.STATUSES` gains **`'placeholder'`** — "a value chosen only to make the geometry
visible; it carries no physical claim". It is distinct from `'estimate'` (which has a stated
method and an anchor) and from `'unspecified'` (which understates the case when the value is
1.0 by fiat). `_TO_MHZ` gains the unit string `'(MHz/(V/cm))^2/MHz'` with factor 1.0, the
same way `'MHz/(V/cm)'` and `'MHz/(e cm)'` already work.

### 5.1 Shared and ²³²Th¹⁹F⁺

`thf_v1()` is unchanged and keeps its name. `thf_v2('232')` returns a ParamSet whose values
are `thf_v1()`'s, plus the Th knobs at zero, plus the two-photon α's — and a gate requires
`thf_v2('232')` and `thf_v1()` to agree on every v1 symbol.

For 229/227, the shared measured/derived constants are unscaled transfers from
232ThF+, not measurements of the selected isotope. They retain their numerical
values but carry `status="estimate"` and no target-isotope uncertainty; the source
measurement uncertainty is retained in the note. In particular B0 should change
with reduced mass (B&C Eq.7.199), so these transfers are not precision isotope
predictions. See the 2026-09-08 parameter estimate audit.

### 5.2 ²²⁹Th¹⁹F⁺ (I_Th = 5/2)

| symbol | value | unit | uncert. | status | source / note |
|---|---|---|---|---|---|
| `A_par_Th` | −1510 (signed; `thf_v2(..., a_par_th_sign=)` selects the trusted calculation) | MHz | unquantified combined error | `ab-initio` | Skripnikov & Titov 2015 Table II FINAL(ThF⁺) −4163 (μ/μ_N) MHz and Denis 2015 +1833 MHz, both rescaled to μ = 0.366(6) μ_N ⇒ −1524 / +1491 MHz; mean of the two rescalings ([TH] §2.2). The historical 60 MHz spread is not a complete bound; the source gives a separate 7% theory-error scale. Note: **sign UNVERIFIED, gap G4** |

R14 (fix round 1): `A_par_Th`'s value is SIGNED, not a magnitude paired with a `conventions.a_par_th_sign` fork -- `params.thf_v2(isotopologue, *, a_par_th_sign="negative")` picks which ab initio calculation (Skripnikov & Titov 2015 vs Denis 2015) the sign comes from; this is a parameter choice, not a convention, so it does not appear in §4's table.
| `g_N_Th` | 0.1464 | — | 0.0024 | `derived` | μ(²²⁹Th)/I = 0.366(6)/(5/2) ([TH] §1.2, Porsev 2021 arXiv:2107.14723). Note: the 1974 value 0.46(4) still in ENSDF is superseded and must never be used to rescale a published A∥ |
| `eQq0_Th` | −2600 | MHz | unquantified | `placeholder` | HfF+ magnitude anchor only. Electronic EFG transfer bracket and signed Petrov-to-B&C conversion are not validated. See the 2026-09-08 quadrupole estimate audit; no calibrated ±1000 MHz error bar. |
| `eQq2_Th` | +300 | MHz | unquantified | `placeholder` | Hf-specific orbital/radial estimate transferred without a Th calibration; signed Petrov-to-B&C normalization remains unresolved. No calibrated ±100 MHz error bar. |
| `c_I_Th` | 0.0 | kHz | — | `held-fixed` | **No value anywhere — gap G3.** [TH] §4.6 declines to pick one and brackets it at **~1 kHz to ~1 MHz**: g_N(Th)/g_N(F) = 0.0279 pushes down, the 2700× larger electronic hyperfine factor on Th pushes up. The bracket is in the `note`, not in the value. **OPEN-19** |
| `Q_Th` | 3.11 | e·b | 0.02 | `measured` | Porsev 2021, weighted average over four Th³⁺ states ([TH] §1.3). Carried for provenance; the Hamiltonian consumes `eQq0_Th`/`eQq2_Th`, not Q |

### 5.3 ²²⁷Th¹⁹F⁺ (I_Th = ½, tentative)

| symbol | value | unit | status | source / note |
|---|---|---|---|---|
| `A_par_Th` | 7620 (magnitude) | MHz | `placeholder` | **A scaled placeholder, and a free knob.** Electronic factor from ²²⁹Th, `G_el = A∥/g_N = −4163 × 5/2 = −10 408 MHz` ([TH] §2.2), times an **ASSUMED** `g_N(²²⁷Th) = μ/I` with **μ(²²⁷Th) assumed equal to μ(²²⁹Th) = 0.366 μ_N** ⇒ g_N = 0.732 ⇒ 7.62 GHz `[derived from an arbitrary assumption]`. **No μ(²²⁷Th) exists**: Stone INDC(NDS)-0794 (2019) and the earlier Stone compilation both skip A = 227 at Z = 90 entirely, and the IAEA NDS live moments database lists only ²²⁹Th and ²³²Th (`docs/lit/lookup-227th-nuclear-moment.md`, with a ²²⁹Th positive control confirming the extraction reached the Th block). The notebook sweeps this knob rather than quoting it. **OPEN-20** |
| `eQq0_Th`, `eQq2_Th` | — | — | — | **Structurally absent, not zero-valued.** A rank-2 nuclear operator has no matrix element for I ≤ ½ — the `(I 2 I; −I 0 I)` in the denominator of B&C (9.53) has no allowed triangle ([TH] §1.4). The terms must return exactly 0.0 from the *formula*, and a gate checks that, rather than the parameter being set to zero |
| `g_N_Th` | 0.732 | — | `placeholder` | same assumption as `A_par_Th`; used only by the Th nuclear Zeeman, a kHz-scale term |
| `c_I_Th` | 0.0 | kHz | `held-fixed` | gap G3, as for ²²⁹ |

**The ENSDF (1/2⁺) is parenthesised — tentative — and the first excited state is a 5/2⁺
only 9.3(3) keV away ([TH] §1.4).** If the assignment is wrong, ²²⁷ThF⁺ acquires a
quadrupole and the whole row changes. The spec's `Spin.I` is a knob, so this costs one
argument to test. **OPEN-20.**

### 5.4 Two-photon polarisabilities

| symbol | value | unit | status | note |
|---|---|---|---|---|
| `alpha_K0_dOm0` | 1.0 | `(MHz/(V/cm))^2/MHz` | `placeholder` | **No ThF⁺ two-photon polarisability exists in any source read** ([2γ] §3.0, §5 gap 1: no published Raman/two-photon rate or line strength for a transition *within* X ³Δ₁, for either molecule, with a positive control confirming the query shape). The value 1.0 exists to make the geometry plottable |
| `alpha_K2_dOm0` | 1.0 | same | `placeholder` | as above |
| `alpha_K2_dOm2` | 1.0 | same | `placeholder` | the ΔΩ = ±2 channel, which requires an Ω = 0 intermediate ([2γ] §3.3) |
| `alpha_K1_dOm0` | **not registered until the derivation settles §3.2** | same | — | if K = 1 is identically zero in closure, the channel does not exist and no `Param` should imply it does |

---

## 6. Verification strategy — the new gates

Same three tiers as v1, same rules: every gate names the failure mode it uniquely catches,
both PASS and FAIL are reachable, no snapshot/hash/pinned-spectrum tests, constants are never
tuned to make a gate pass. New ids continue [HAM] §6's V-series.

| id | check | uniquely catches | FAIL reachable by |
|---|---|---|---|
| **V16** | **Master reduction.** With I_Th = 0: the v2 enumerator's kets, F₁ column dropped, equal `enumerate_kets(thf_spec())` element-for-element; and every v2 term's matrix equals its v1 counterpart to 1e-12 on the same block | **any wrong spectator phase or 6j column order in any of the seven re-coupled terms** — the single highest-value gate in v2 | transposing two columns of either B&C (5.174) 6j |
| V17 | Dimension closed form `(2I_Th+1)(2I_F+1)Σ_J 2(2J+1)`; blocking is a partition; every value a multiple of ½; \|Ω\| ≤ J; F₁, F triangles | enumerator bugs in the added F₁ loop (an off-by-one in the F₁ range is otherwise invisible until a spectrum looks odd) | admitting F₁ outside \|J−I_Th\|…J+I_Th |
| V18 | **Th hyperfine manifold pattern.** The diagonal B&C (9.50) coefficients over F₁ reproduce [TH] §4.1 at J = 1–4 (e.g. J = 1: −1.7500, −0.5000, +1.2500) | mis-coupling I_Th to F instead of to J — which gives a plausible-looking manifold with the wrong spacings | swapping F₁ for F in the closed form |
| V19 | **¹⁹F doublet-ordering flip.** At I_Th = 5/2, J = 1, the F = F₁+½ minus F = F₁−½ splittings in units of A∥^F are −0.4000, +0.1714, +0.5714 over F₁ = 3/2, 5/2, 7/2 ([TH] §4.5); the I_Th = 0 reference is +0.7500 for every F | a wrong phase in the B&C (5.173)/(5.174) recoupler of the outer spin — the sign *inverts* at F₁ = 3/2, so this is a sign test, not a magnitude test. Also the falsifiable ²²⁹ThF⁺ prediction the notebook displays | dropping the (−1)^{F₁′+F+I_F} phase |
| V20 | **Quadrupole ↔ Casimir.** B&C (9.53) at Ω = 0, J′ = J reproduces the textbook Casimir function with a ratio of exactly **−1** for (J, I) = (1,1), (2,1), (2,3/2), (3,5/2), all F ([TH] §3.4) | the q₀-sign convention being dropped — which flips the sign of every quadrupole splitting and nothing else notices | removing the explicit `−1` |
| V21 | **Quadrupole vanishes at I ≤ ½.** Both quadrupole terms return exactly 0.0 for I_Th = ½ and I_Th = 0, from the formula | a code path that computes a spurious quadrupole for ²²⁷ThF⁺ by dividing by a 3j that should not exist | returning the q = 0 limit instead of zero |
| V22 | **eQq₂ is ΔΩ = ±2, diagonal in J, F₁, F, m_F, and parity-even**; its J = 1 coefficients over F₁ are +0.1715, −0.1960, +0.0612 ([TH] §4.4) | confusing it with the Ω-doubling operator, which occupies the same matrix position with a different (J, F₁) dependence | giving it the ω_ef J(J+1)/4 J-dependence |
| V23 | **J-truncation is reported and non-vacuous.** For ²²⁹ThF⁺ the largest J = 1 level shift between J_max = 2 and J_max = 4 exceeds 1 MHz; for ²³²ThF⁺ the same quantity is below 1 kHz | a silently dropped Th ΔJ = ±1 hyperfine, or a `J_max` knob that is not wired through — [TH] §4.2 puts the Th ΔJ = ±1 second-order shift at 60–135 MHz against 2.6 kHz for ¹⁹F | pinning `J_max` inside the enumerator |
| V24 | **Two-photon parity.** `[T_eff, P] = 0` at zero field for every (K, ΔΩ, P) channel and every polarisation pair | a sign error in the ΔΩ = ±2 geometry, which would fake an e↔f two-photon line that does not exist | flipping the relative sign of the two q = ±2 components |
| V25 | **Channel reach.** ΔΩ = ±2 is non-zero only for K = 2; Δm_F ⊆ {0, ±1, ±2} and never ±3; under σ±-only weights Δm_F ⊆ {0, ±2} | a rank-K contraction that leaks a molecule-frame \|q\| > K component | allowing q = ±2 at K = 0 or 1 |
| V26 | **Rank-2 sum rule and reciprocity** — the B7/B7b analogues: total two-photon strength out of a state summed over final states and P is independent of m_F; and `⟨F′‖α^K‖F⟩` obeys spherical-tensor reciprocity against `⟨F‖α^K‖F′⟩` | a missing (−1)^{F′−m′_F} in the lab-frame Wigner–Eckart 3j, which leaves every magnitude right and only the m_F dependence wrong (this is exactly what V14 catches for the one-photon Stark term) | dropping the phase |
| **V27** | **Closure test (§3.3).** A fake intermediate manifold, resolved amplitudes summed with a common detuning, equals the rank-K form to machine precision | **a wrong 6j or phase in B&C (5.142)** — invisible to V24, V25 and V26, all of which survive it | transposing two columns of `{1 1 K; j′ j j″}` |
| V28 | **Sum-then-square, two-photon.** Two coherent (K, P) pathways with opposite-sign amplitudes cancel; squaring first does not | the incoherent-sum mistake, which Cossel *measured* the consequences of — his two σ pathways cancel "because of the signs of the Wigner 3j coefficients" ([2γ] §3.0) | squaring inside the channel loop |

**Tier D (opt-in, `HEFF_RUN_LITERATURE=1`)** gains nothing new, deliberately: **no Th
hyperfine constant of any kind has been measured for any ThF⁺ isotopologue** ([TH] §2.1,
gap G5), and no two-photon rate within X exists ([2γ] gap 1). So there is nothing to compare
against, and inventing a comparison against an ab initio number rescaled by this same code
would be circular.

---

## 7. Recorded physical uncertainties and parameter choices

Numbering continues [HAM] §7. Each identifier marks a place where the available
sources do not support a unique physical value or model choice.

**OPEN-16 — the sign of A∥(Th).** Skripnikov & Titov print −4163 (μ/μ_N) MHz, Denis print
+1833 MHz from the same defining equation, agreeing to 2.2 % in magnitude ([TH] §2.2, gap G4).
The sign is what **orders the F₁ manifold**, so every ²²⁹ThF⁺ level diagram depends on it.
[The sign-convention audit](../../lit/lookup-apar-th-sign-convention.md) shows
the two groups' axis conventions are opposite but that A∥ is invariant under a consistent
reversal, and that they agree on the analogous HfF⁺ constant — so **the disagreement is a real
disagreement between two calculations, not a convention mismatch**, and the audit recommends
A∥ < 0. §5 therefore records the negative branch as the default while retaining
the positive branch for explicit comparison. The convention audit does not resolve
the disagreement between the two electronic-structure calculations.

**OPEN-17 — the eQq₂ normalisation bridge.** B&C's (9.52) at q = ±2 versus Petrov 2018
Eq. (23), with its √6 and Y₂₂. [TH] §3.2 flags this **UNVERIFIED** and says every
eQq₂-derived number inherits the caveat. Until an explicit bridge is derived from
the printed definitions or confirmed by the authors, the HfF⁺-anchored value is an
order-of-magnitude estimate.

**OPEN-18 — eQq₀ and eQq₂ defaults.** There is **no published ThF⁺ or ThO quadrupole
coupling constant and no EFG at Th, for any isotope or state** ([TH] gap G2, with the
queries recorded). §5.2 defaults them to the HfF⁺-anchored estimates because a zero default
would hide a term that [TH] §4.4 argues **dominates the Ω-doublet structure of ²²⁹ThF⁺**.
The estimates therefore carry `status='estimate'` and must remain distinguishable
from measured molecular constants.

**OPEN-19 — c_I(Th).** Unconstrained over three decades, ~1 kHz to ~1 MHz; [TH] §4.6
declines to pick, and at the top of the range it would exceed the entire ¹⁹F hyperfine.
The reference default is zero; any bracket sweep is a labeled sensitivity study.

**OPEN-20 — ²²⁷Th: the spin and the moment.** I = (1/2⁺) is a *tentative* ENSDF assignment
with a 5/2⁺ level 9.3 keV above it, and **no magnetic moment exists in any of the three
compilations checked** ([lookup](../../lit/lookup-227th-nuclear-moment.md)). §5.3's
A∥ = 7.62 GHz rests on the arbitrary assumption μ(²²⁷Th) = μ(²²⁹Th). Alternative
spin or moment choices define separate sensitivity models, not updated nuclear data.

**OPEN-21 — does the K = 1 two-photon channel exist?** §3.2's argument says the
antisymmetric part vanishes in exact closure and survives only through detuning asymmetry.
That is an inference. A K = 1 channel requires a derivation beyond exact closure;
the closure model exposes only K ∈ {0, 2}.

**OPEN-22 — J_max for ²²⁹ThF⁺.** The Th ΔJ = ±1 hyperfine is a ~2 GHz off-diagonal element
whose second-order shift is 60–135 MHz and, unlike the ²³²Th case, **is not absorbable into
B₀ because it depends on F₁** ([TH] §4.2). [TH] infers that J = 1–5 or 1–6 may be needed for
kHz-level J = 4 energies. The package makes J_max a `StateSpec` knob; reported
results must state the cutoff and include a convergence comparison appropriate to
the claimed precision.

**OPEN-23 — the two-photon validity condition, in the notebook's own voice.** The rank-K
form is the Δ ≫ 2B ≈ 7 GHz limit; JILA runs at 0.16–1.5 GHz ([2γ] §3.4).
Any resolved-sum calculation must identify which contested intermediate ladder it uses.

---

## 8. Backward compatibility — what v1 keeps, exactly

| v1 surface | v2 status |
|---|---|
| `KET_C`, `enumerate_kets(thf_spec())` → 96 kets in the same order | **unchanged, byte-identical** |
| `thf_spec()` (no arguments) | **unchanged**; `thf_spec('232')` must compare equal to it |
| `heff.terms.REGISTRY` and its eleven terms | **unchanged**; v2 terms live in `REGISTRY_C2` |
| `elements_c.py` | **unchanged**; v2 re-uses `dipole_geometry` and the two hyperfine closed forms through adapters |
| `thf_v1()` | **unchanged**; `thf_v2('232')` must agree with it on every v1 symbol |
| `Conventions` | four fields added; `stamp()` gains four keys |
| `Rules` | one optional field `dF1=None`; `allows()` skips it when `None` or when the dtype has no `F1` |
| `Ctx` | one field `spins=()` with a default |
| `conventions.parity_operator` | keys on every dtype field instead of four literals; identical output on `KET_C` |
| `spectra.dipole_matrix` | gains `geometry=` keyword defaulting to today's behaviour |
| `assemble`, `engine`, `track`, `observe`, `wigner`, `formalism` | **untouched** |
| the original scientific checks | remain applicable to the one-spin limit |

---

## 9. Risks and unknowns, ranked

1. **Every ²²⁹ThF⁺ number is only as good as an unmeasured constant.** A∥(Th) is ab initio
   with an unresolved sign; eQq₀ and eQq₂ are HfF⁺-anchored estimates through a bridge
   flagged UNVERIFIED; c_I(Th) spans three decades. Nothing in the package is wrong because
   of this — but no ²²⁹ThF⁺ level diagram is a prediction in the sense a ²³²ThF⁺ one is, and
   the notebook must say so on the figure, not in a footnote.
2. **The recoupling phases are the transplant risk, exactly as the convention transplant was
   in v1.** Mitigation is V16 (master reduction), V19 (a *sign*-sensitive prediction, not a
   magnitude one), and the fact that Th and ¹⁹F reach the same physics by two independent
   routes — closed form with F → F₁ for the inner spin, genuine recoupler for the outer —
   so a phase error common to both is the only residual.
3. **eQq₂ versus ω_ef versus the Stark energy.** At the estimated 200–400 MHz, eQq₂'s
   Ω-mixing element is 35–80 MHz at J = 1, against ω_ef J(J+1)/4 = 2.65 MHz **and** against
   the 50.9 MHz Stark polarisation energy at 60 V/cm ([TH] §4.4). Two of the three are the
   same size, so the ²²⁹ThF⁺ doublet structure is a genuine competition, not a hierarchy —
   which is why the notebook sweeps it rather than quoting a level diagram.
4. **The two-photon operator is in the wrong limit for the real experiment** (§3.1). This is
   a stated validity condition, not a bug, but it is the thing most likely to be
   over-read from a pretty spectrum panel.
5. **Basis growth.** 576 states (1152 at J_max = 6) is 6–12× v1. The per-block build is the
   only cost that grows, and it is paid once. If it becomes the bottleneck the sparsity mask
   and the block-aware cache key are already in place; nothing needs redesigning.
6. **Two dtypes is a seam.** Anything that reads ket fields by position, or that assumes four
   fields, breaks silently on `KET_C2`. Mitigation: read by name everywhere (v1 already
   does), and the `parity_operator` generalisation is the one place v1 assumed four literals.
7. **`status='placeholder'` could be mistaken for a value.** It is added precisely so the
   report layer can distinguish it, and nothing gates on status (v1 rule, unchanged) — but a
   plot that renders a placeholder α without a label would mislead. The notebook labels it.
