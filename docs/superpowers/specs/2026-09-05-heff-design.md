# Reference architecture — molecular effective-Hamiltonian package

This document records the initial `heff` architecture adopted in September 2026.
Its scientific basis is the [thesis convention digest](../../digest-thesis-effective-hamiltonian.md)
and [ThF⁺ literature digest](../../digest-literature-thf-plus.md). The maintained
[architecture guide](../../architecture.md) describes the implemented surface.

---

## 1. Problem restatement, goals, non-goals

**Problem.** Two codebases hold complementary halves of an engine and neither can compute ThF⁺. Molecule-Structure has the angular-momentum algebra behind a nine-key dispatch table that pins S = 1/2 and rebuilds every matrix element on every parameter change (0.88 s at dim 64, digest §7). C2V-Molecules has a Hamiltonian-agnostic scan/derivative engine behind an asymmetric-top basis it cannot leave (four coupling lines, digest Q1). Neither has case (c) or a typed parameter record.

**Goal — v1 slice.** ²³²Th¹⁹F⁺ X ³Δ₁, v = 0, J = 1–4, ¹⁹F hyperfine (I = ½), Ω-doubling, Stark (E_z), Zeeman (B_z), per-m_F block, over a user-supplied collection of parameter sets. Outputs: eigenvalues and eigenvectors; g-factors and molecule-frame dipoles by exact Hellmann–Feynman derivative; the Ω-doublet differential Δg = g^u − g^ℓ; PT-odd matrix elements (W_d, W_T,P); E1 TDMs within X. A parameter change is a weighted resum of cached term matrices, never a matrix-element rebuild.

**Goal — beyond v1.** Case (b), the linear-polyatomic ket |v, ℓ, K, P⟩ and the asymmetric top are new modules against an unchanged spine. The thesis convention contract (digest §4, 21 items) is stated and tested once, not re-derived per element.

**Non-goals for v1.** No GPU (the slice is ~96 kets over 10 m_F blocks, §3.3). No fitting (the analytic Jacobian stays available, §3.6). No disk cache, dynamics, OBE, plotting or CLI. No asymmetric top, no Renner–Teller. No non-collinear field *wired into a production sweep* — but E_x, B_x, ω_rot F_x and AC intensity are declared knobs from day one, because the JILA scheme is a rotating bias field (lit §3.1) and retrofitting the blocking layer is expensive (§3.1).

---

## 2. Three candidate architectures

### 2.0 Literature basis for the initial case-(c) slice

The literature digest landed mid-draft and moves the v1 target basis. The JILA/Ng treatment of ThF⁺ X ³Δ₁ is **Hund's case (c)**: the basis is the product |J, Ω = ±1, F, m_F⟩ (Ng thesis App. C.3.1, p. 323 — lit digest §3.1), with rotation `B J(J+1)` and **no −Ω² term** (lit §5.12), hyperfine collapsed to the single projection-theorem constant A∥ (Eq. C.2), and Ω-doubling as the closed-form phenomenological operator `(−1)^J (ω_ef/2) · [J(J+1)/2] (|+1⟩⟨−1| + h.c.)` (Eq. C.3, lit §3.1 item 2).

The literature review established three consequences:

1. **The "missing ³Δ₁ Ω-doubling operator" is not a v1 blocker.** Synthesis §1 and Molecule-Structure digest §HEADLINE both rank it as the one genuinely missing piece of physics. In case (c) it exists in closed form with a stated J-scaling and a stated parity alternation. It only reappears as an open item if Arian wants the microscopic case (a) form `½(o_Δ + 3p_Δ + 6q_Δ)(S₊²J₊² + S₋²J₋²)` (Leanhardt 2011 Eq. 16, lit §3.2). That demotes risk #1 in the synthesis to a *later-milestone* physics question (§4 Q2, §5 R2).
2. **The v1 basis is one neither repo has.** Molecule-Structure supports exactly `aBJ`, `bBJ`, `bBS` (`molecule_library_class.collect_all_cases:61–82`, digest §HEADLINE(c)). Case (c) is a *new* subengine — but a cheap one: no Σ, no Λ, no 9j, ~96 kets for J = 1–4. This is good news for milestone ordering (§3.10): M1 does not depend on the lifted case (a) library at all.
3. **A single-electronic-state effective Hamiltonian provably cannot reproduce the headline observable.** Ng's own 32-level model gives δg/g = −0.00223 against a measured −0.00255(6), a 15 % gap he attributes to unmodelled X ³Δ₁ – ³Δ₂ coupling (Ng thesis p. 85); Petrov & Skripnikov's model, which carries ¹Σ⁺, 1³Δ₂, ³Π₀± in the basis, lands on the measurement (lit §7.8). **Whether the state spec admits several electronic states with off-diagonal couplings is therefore an architectural question, not a physics detail** — it is the one thing on this list that is a rewrite if it is not designed in. All three candidates below are judged partly on that.

### 2.1 Gaps in the earlier four-layer sketch

The four-layer sketch (state spec → operator library → assembly → engine) is right in outline and wrong or silent in seven places, each of which is a place a real bug lives in one of the two repos.

1. **No block object.** §2 conflates "the basis" with "the matrix you diagonalize". C2V learned otherwise: `generate_states_and_blocks` returns `{m_F: array}` and every stage downstream treats one block as an independent problem (digest §2b). Molecule-Structure has no blocking at all and pays ∏(2J+1) dimension for it (digest §2d). But blocking cannot be a hidden loop, because two v1 requirements cross blocks: a transverse or rotating field couples Δm_F = ±1, and an E1 TDM connects different m_F. A `Blocking` must be a named, invertible partition carrying index maps back into the full basis, with a `merge()` that reproduces C2V's concatenated-with-offsets basis (`transverse.py:83 build_full_basis`, digest §2d). Missing this is the largest structural gap.
2. **Fields and parameters are typed as different things.** `H(p, f) = Σ p_k M_k + Σ f_j V_j` is a false dichotomy: a field knob is a parameter whose matrix happens to be a field operator. Collapsing them into one flat catalogue buys three things at once — `multi_curvature` over any subset of knobs (already true in C2V, `qgt.py:180`), analytic ∂E/∂p_k for *Hamiltonian* parameters (i.e. the fit Jacobian, free, §3.6), and one cache instead of two. The distinction that actually matters — swept on a grid vs held fixed — is per-call, not per-type.
3. **Observables have no home.** §2 lists them as engine features. They are matrices in the same basis, built by the same machinery: Σ = S·n̂, ⟨M_S⟩, parity, W_d. Put them in the catalogue as a second named family and the engine takes `{name: matrix}`; leave them in the engine and you reproduce C2V's exact wart — `obs_ops = [stark, zeeman, EDM_mat]` hard-wired identically at `scan_2d.py:293`, `scan_windowed.py:120` and inside `repair._stream_diag_core` (verified this session), so a fourth observable means editing three drivers.
4. **The convention contract has no module.** The thesis's 21 items (digest §4) are properties of the basis phases and the operator library jointly. Molecule-Structure re-implements them per matrix element and carries nine `#check` / `#Check derivation` comments as a result (digest §3). One `conventions` module that owns every phase, provides the parity and time-reversal operators, and is tested once, is the whole difference.
5. **Disk caching is premature and mis-specified.** "Cached on disk keyed by the state spec" invents an invalidation surface before anything is slow. Worse, it collides with Arian's own rule: on a key mismatch a record-and-warn cache must then either use the stale entry or silently rebuild, and neither is stated. In-process `TermMatrices` covers the stated workflow; the manifest and key are designed now, the disk backend is a one-function swap later (§3.3).
6. **Ordering, tracking and labelling are one bullet and are three things.** The thesis rule is adiabatic correlation *to the zero-field state* (p.267, digest §5); C2V production runs energy order with no matcher and certifies labels separately; C2V's hardest-won lesson is that the label convention must travel *inside* the artifact (wart #9). These are an ordering policy, a tracking algorithm, and an artifact field. Separate them or they fuse (§3.5).
7. **The parameter record is a footnote inside assembly.** It caused both repos' worst live bugs — units living in a script (C2V `d_0` reduced-vs-physical, a silent 6.4 % error, digest wart #3) and sign conventions living in prose comments (`molecule_parameters.py:64, 96, 161`). The literature digest adds five more convention forks for ThF⁺ alone (n̂ direction, dipole origin, δg vs Δg factor 2, eEDM factor 2, Zeeman sign — lit §5.2, 5.4, 5.10, 5.11, 5.8). It is its own layer, upstream of assembly (§3.4).

### 2.2 The three candidates

They differ in **what the unit of composition is**, which determines what "add a new coupling case" costs.

**A — Term-matrix catalogue (dense-matrix-first).** Unit of composition: a parameter-free matrix `M_k`. A basis is a flat, ordered list of kets; an operator is a pure function of a (bra, ket) pair; a term registry maps name → (function, parameter symbol, applicable cases, declared selection rules, citation). Assembly is `H = Σ_k c_k M_k`. Adding a case = one enumerator + one set of element functions registered against the same spine. This is what QuantumStates.jl converged on (`Hamiltonian.jl`: `matrix = Σ parameters[op.param] .* op.matrix`, prior-art §2) and what both of Arian's repos already are underneath (`hamiltonian_builders.py:99` `H_func = lambda E,B: H0_num + V_E_num*E + V_B_num*B`; C2V `physics.py:617`).
*Pros:* parameter change is one BLAS call; term toggling is free, which is exactly the thesis's own deperturbation workflow ("set H_K = 0 in the model", digest §6); ∂H/∂c_k = M_k, so every derivative is exact and every fit Jacobian is analytic; the whole C2V engine plugs in unchanged because it already takes matrices. *Cons:* memory is O(n_terms · dim²) — irrelevant at dim ~10² but a real ceiling at dim ~10³–10⁴ for polyatomic full-M bases (mitigation: sparse per-term storage, dense assembled H, §3.3); and adding a coupling case still means writing that case's matrix elements by hand.

**B — Operator-algebra / spherical-tensor DSL.** Unit of composition: a symbolic tensor expression (`T¹(J)·T¹(S)`, `e^{2iφ}T²_{2q}(J,J)`, …) reduced to numbers at evaluation time by a generic Wigner–Eckart/spectator engine that knows the coupling tree. Adding a coupling case = declaring a coupling tree; matrix elements fall out.
*Pros:* the only candidate where a new case costs no new matrix elements; the coupling-order bookkeeping the thesis flags as dangerous (item 14: the spectator phases differ between the A₁-on-J₁ and A₂-on-J₂ forms, Eqs. A.39 vs A.40) gets written once instead of per element; it is what PGOPHER is, and PGOPHER is the only existing tool with genuine multi-case generality (prior-art §3).
*Cons, and they are decisive for v1:* you must get the reduction engine right before you can compute a single ThF⁺ number, and the canonical reference for the hardest part of it is known-wrong — the thesis explicitly warns against B&C Eq. 5.177 as "missing the extra factors from the spectator theorem" (digest §4 item 20b) and flags a typo in B&C's case (a)↔(b) transform (item 15) and in Hirota Tab. 2.4(2) (item 20a). Building a compiler against three known-buggy sources, with no numbers to check it against until it is finished, inverts Arian's verification style: closed-form limits validate outputs, not algebra layers. Rejected for v1.

**C — Basis-object subengines.** Unit of composition: a molecule-class module owning its ket type, enumerator, term list *and* assembler. Shared layer is only the kernels — Wigner, eigh/scan, tracking, observables, gates.
*Pros:* literally Arian's Rule 1 (distinct subengines, shared tools); matches the physical reality that case (c) and the asymmetric top share almost no matrix elements; each subengine can be verified against its own literature independently.
*Cons:* there is no shared *term* layer, so rotation, Stark, Zeeman and hyperfine get written four times, and — the real cost — the phase conventions get written four times with them. That is the status quo across the two repos, and it is why the same physical Λ-doubling term is correct in QuantumStates.jl and identically zero in Molecule-Structure (§3.2).

### 2.3 Selected architecture

**Take A as the spine; draw C's subengine boundary at exactly the triple (ket dtype, enumerator, element set); adopt B's declarative metadata and reject B's reduction engine.**

The reconciliation with Rule 1 is not a compromise, it is the literal reading: the *subengines* are the coupling cases (case (c), case (a), case (b), linear-polyatomic vibronic, asymmetric top); the *shared tools* are the Wigner kernel, the term registry, the blocking layer, the assembler, the scan engine, the observables, the artifact format and the gate suite. Nothing is collapsed into one implementation; everything below the ket dtype is shared.

The one piece of B worth taking now is cheap and pays immediately: **every registered term declares its selection rules as data, separately from its formula.** That single field turns an entire bug class into an automatable, Hamiltonian-agnostic gate — "the formula must be non-zero somewhere inside its declared-allowed set and exactly zero outside it". `LambdaDoubling_q_even_aBJ` (`matrix_elements.py:644`, read this session) fails that gate today: it gates on ΔΛ = +2q while its 3j `wigner_3j(J0,2,J1,-P0,2*q,P1)` forces ΔP = −2q, so with ΔΣ = 0 it is identically zero everywhere, and BaF A0, CaOH A000 and YbOH A000 all carry a non-zero `q` doing nothing (digest §3, verified by sweep with a Λ = 1 control giving 20 non-zero). QuantumStates.jl's version of the same operator uses −2q and works (prior-art, synthesis §3). **Both PASS and FAIL are reachable with real inputs, and we have one of each** — which is the burden-of-proof standard, satisfied in advance.

---

## 3. Reference design

### 3.1 Basis and state spec

**Decision.** A frozen `StateSpec` dataclass, per coupling case, plus a registered enumerator. The spec is data; the enumerator is the only code that knows the case.

```
StateSpec(case='c',                       # selects enumerator + element set
          electronic=[ElecState(label='X3Delta1', Omega=1, S=1, Lam=2, T0=0.0)],
          v=0, J=(1,4), spins=[Spin('19F', I=0.5, couple_to='J')],
          M='blocks', frame='rotating',   # see below
          parity='signed')                # signed-Ω primitive basis
```

**Kets are a structured numpy array with named fields**, `np.dtype([('J','f8'),('Om','f8'),('F','f8'),('mF','f8')])`, not a positional tuple and not a list of objects. This is a deliberate middle between the two repos: it keeps C2V's float64 discipline (their `(n,7)` array *must* be float64 or numpy ≥ 2 NEP-50 leaks ~1e-7 into every matrix element, `physics.py:161–169`) while removing C2V's positional-7-tuple opacity, and it makes selection-rule masks vectorised and the cache key trivially canonical. Invariant, gated: every stored value is an exact multiple of 0.5 (§3.7 A7). *Trade-off:* named-field access is ~2× slower than positional in the innermost loop; irrelevant, because the loop runs once per basis, not once per parameter set.

**³Δ₁ full manifold vs Ω = 1-only.** Three expressible options, and the choice is Arian's (§4 Q1):
- case (c), `electronic=[³Δ₁]`, `Omega=±1` — the Ng/JILA basis; 12 kets at J = 1, 32 at J ≤ 2, **96 at J ≤ 4**; Ω-doubling is the phenomenological `ω_ef` term.
- case (c), several `ElecState`s (³Δ₁, ¹Σ⁺, 1³Δ₂, ³Π₀±) with off-diagonal electronic terms — Petrov's basis; needed to reproduce δg (lit §7.8).
- case (a), `S=1, Lam=2, P_values=[1]` or `P_values=None` — the full ³Δ manifold with S-uncoupling. Molecule-Structure's `q_numbers_even_aBJ` (`quantum_numbers.py:320`, verified) already builds this: the digest ran it and got 16 kets for J = 1–4, Ω = 1.

The `electronic` field being a **list**, with a term-registry category for inter-electronic-state operators, is the one design decision here that is expensive to retrofit and nearly free to include: it costs one extra field in the ket dtype and one extra selection-rule axis. Given lit §7.8 says the headline observable needs it, include it, enumerate one state in v1.

**Later cases plug in by adding fields, not by changing the spine.** Case (b) adds `N, K`; linear polyatomic adds `v2, l` with `K = Λ + ℓ` and `P = Λ + Σ + ℓ` as declared constraints (thesis digest §7); asymmetric top adds `N, K_a/K_c` or C2V's `N, K, J, F_N, I_T`. Copy QuantumStates.jl's `constraints` idea as a list of declarative relations evaluated by the enumerator (prior-art §5.2) — `K = Λ+ℓ`, `Σ ∈ −S:S`, `P ∈ max(−J, K+Σ):min(J, K+Σ)`, `F ∈ |J−I|:|J+I|`, `m ∈ −F:F`. Verified: those are free `HalfInt` fields there, which is why S = 1 / Λ = 2 needs zero new struct code.

**Parity: signed-Ω (or signed-K) primitive basis plus an optional parity projector.** Both repos enumerate both signs and treat parity as an operator (`Energy_Levels.__init__:167`); the thesis works in parity-symmetrised kets (Table 4.2). *Trade-off:* the symmetrised basis halves the dimension and block-diagonalises, but every matrix element then needs a phase-composed variant, and the thesis's own disambiguation rule for v₂ = 3 ("write ℓ₁ states with Λ = +1 as the first ket…", digest §7) becomes load-bearing bookkeeping in the element layer. *Decision:* keep the signed primitive basis so the element functions stay in the form both repos already have, and ship the parity projector as an explicit orthogonal transform that can be applied either as a blocking or ignored. You pay 2× dimension when you don't project; at dim 96 that is nothing, and at dim 10⁴ the projector is available.

**M handling is the blocking layer, and it is a first-class object.** `M ∈ {'none', 'blocks', 'all'}` produce respectively a single m-free basis (field-free spectra), a `Blocking` of per-m_F sub-bases, and one dense basis. `Blocking` carries `{label: index_array}` into the full basis plus `merge(labels)`. Per-m_F is the default fast path (the thesis's own trick — "separately diagonalize each M_F block to avoid degeneracies at B_Z = 0", p.267); full-M is the general path, entered automatically when any active term violates Δm_F = 0.

**Frames.** The spec declares which frame m_F is quantised in. This is not pedantry: Ng defines m_F in the frame rotating with E_rot, Petrov uses M_F = M_I + M on the lab ẑ and explicitly warns "M_F is not equal to m_F" (lit §5.1). A code that mixes them mislabels every Zeeman sublevel silently.

*What could go wrong.* (i) Two specs differing only in ket ordering produce different matrices and a false cache hit — so ket order must be a deterministic, tested function of the canonicalised spec. (ii) The polyatomic |K| ambiguity at v₂ = 3 (thesis digest §7) means a spec listing only |K| is ambiguous; the ket must carry (v₂, ℓ, Λ). (iii) Half-integer float equality — exact for these magnitudes, but only if nothing ever divides; gated. (iv) A user asks for `M='blocks'` with a rotating-frame term active and gets a silently wrong answer; the assembler must detect the Δm_F violation from the declared selection rules and refuse to block (this is a structural check, so raising is correct).

### 3.2 Operator library and term catalogue

**Decision: no default values for physical quantities in operator signatures.** Molecule-Structure's `S = 1/2, I = 1/2` keyword defaults are the entire reason S is never threaded: `S_electron` is accepted at `Energy_Levels.py:33`, used only to compute default `P_values` at `:82`, and `grep -rn "e_spin"` returns five lines repo-wide, none of them a matrix-element call (digest §HEADLINE(c)). Every element therefore silently evaluates at S = ½. A signature with no defaults turns that class of bug into a `TypeError` at first call. One rule, whole bug class gone.

Operators stay pure functions `f(bra_row, ket_row, ctx) -> float | complex`, where `ctx` is a frozen record carrying S, Λ, I list, nuclear g-factors, and the convention block. This preserves the property that makes Molecule-Structure's `matrix_elements.py` the highest-value asset in either repo — uniform signature, no hidden state, individually testable (digest §6 as-is #1).

**Term registry entry:**

```
@term(name='hyperfine_A_par', param='A_par', cases=('c',),
      rules=Rules(dJ=0, dOmega=0, dF=0, dmF=0),
      hermitian=True, real=True,
      cite='Ng thesis Eq. C.2 p.319')
def hyperfine_A_par(bra, ket, ctx): ...
```

`rules` earns its place three ways: it gives the assembler a sparsity mask so only allowed (i, j) are evaluated (the dominant build cost); it answers "which terms are non-zero in this basis?" without building anything; and it is one half of gate A5 (§3.7), the formula-vs-declared-rules consistency check that catches dead operators. *Trade-off:* the rules are hand-written and can themselves be wrong. Mitigation: a wrong rule set that is *narrower* than the formula makes the gate fail loudly (non-zero outside the allowed set); a rule set that is *wider* degrades the gate to "non-zero somewhere", which is still the check that catches the ThF⁺-relevant dead-operator case. The failure mode is graceful in the direction that matters.

**Adding a new term is one decorated function plus one parameter entry.** No dispatch dicts. This is the single largest simplification versus Molecule-Structure, where a state outside the pre-existing envelope requires new keys in nine `collect_*` dicts plus relaxing the `iso_state` string construction (digest §4).

**Conventions live in exactly one module** which owns: the parity operator `E* = σ_xz R_y(π)` and the composite phase `(−1)^{J−S−ℓ+s}` (thesis items 7, 9); the two opposite-sign azimuthal phases, `⟨Λ=±1|e^{±2iφ_e}|∓1⟩ = −1` versus `⟨ℓ=±1|e^{±2iφ_n}|∓1⟩ = +1` (item 11 — the thesis calls this the most consequential sign in the code); Condon–Shortley (item 3); anomalous commutation for J and N only (item 5); `δ_ℓ = δ_v = 0`, `η_J = 0` (items 10, 13); the Wigner–Eckart form Eq. 2.6 (item 2); the case (a)↔(b) transform in the thesis's typo-corrected form, *not* B&C's (item 15); the R² convention with `Λ` never inferred (item 16); and `e/f` as `P = ±(−1)^{J−S−ℓ}` (item 21). It also owns the three "the textbook is wrong here" flags (item 20). Nothing else in the package writes a phase.

*What could go wrong.* The load-bearing one is physics, and it is an open question, not a pick: the thesis derives the composite parity phase Eq. A.15 using `S = |Σ|` (thesis digest §10 Q6). A ³Δ manifold has Σ ∈ {−1, 0, +1}, so the Σ = 0 components sit exactly where that simplification was made. Whether `(−1)^{J−S−ℓ+s}` carries over unchanged for Σ ≠ ±S must be answered before any case (a) ³Δ parity label is trusted (§4 Q3). This does not touch the case (c) v1 slice, where parity comes from the Ω-doublet superposition directly (lit §5.6) — another reason to order the milestones case (c) first.

### 3.3 Term-matrix cache and Hamiltonian assembly

**Decision.** `TermMatrices` is an immutable mapping `name → sparse matrix` for one block, plus a manifest (spec hash, term name and citation, conventions version, Wigner backend and version, build time). Assembly is `H = Σ_k c_k M_k`.

**Why the cache is the architecture, with the arithmetic.** Molecule-Structure measured 0.88 s to build dim = 64 with ~14 terms — roughly 57 k sympy-Wigner element evaluations, ~15 µs each (digest §7). The ThF⁺ v1 slice at dim 96 with ~8 terms is ~74 k evaluations, ~1.1 s. Per parameter set. Ten thousand trial parameter sets is three hours of arithmetic that should be a dot product. With the catalogue it is one build plus 10⁴ resums.

**The sweep is a single tensor contraction, and that is the design's payoff.** For a collection of parameter sets, stack coefficients as `c[n_sets, n_terms]` and term matrices as `M[n_terms, d, d]`; then `H_batch = np.tensordot(c, M, axes=1)` gives `(n_sets, d, d)` in one BLAS call, feeding straight into a batched `eigh`. Chunk over `n_sets` on a memory budget (10⁴ × 96² × 8 B ≈ 737 MB, so chunk); the chunk size is a knob, not a constant.

**Cache key and invalidation.** Key = canonical JSON of (canonicalised spec, sorted term names, conventions version). A key miss rebuilds. A key *hit* with a differing code fingerprint (term-function source hash, backend version) **prints a warning and rebuilds by default**, with `allow_stale=True` to use it anyway. Nothing raises on a hash, ever; the fingerprint is a record (Arian's standing rule; C2V's own `_identity_hash` is a cache-validity bool at `scan_io.py:287` and is the correct precedent). Reserve raising for structure: a term matrix whose shape disagrees with the block, a missing term, a non-Hermitian matrix from a term declared Hermitian.

**Storage: sparse per term, dense per assembled H.** Most terms are sparse by selection rule; the assembled H must be dense for `eigh` anyway. This is the whole mitigation for the O(n_terms · dim²) memory ceiling and it costs one `.toarray()` after the sum. Name the ceiling in the code (`# ponytail: dense assembled H; if dim > ~5e3 the sum itself needs chunking`).

**Complex vs real is decided at assembly, not at build.** Each term matrix carries its own dtype; the assembled block is promoted only if a complex term is active. So you never pay a complex `eigh` for a real Hamiltonian, and the NSD-PV operator (Molecule-Structure's only complex one, `build_PTV_NSDPV:321`) and C2V's imaginary-Hermitian `He_y/Hz_y` (`transverse.py:126–133`) both work without a global dtype decision.

**Blocking interacts with the cache**: one `TermMatrices` per block, keyed on (spec, block label). `merge()` produces the concatenated basis with offsets and re-assembles rather than trying to stitch cached blocks — cheap, and it avoids a whole class of index bug.

*What could go wrong.* (i) Memory at large dim (named above). (ii) A term silently absent from a parameter set is assembled with coefficient zero — indistinguishable from "deliberately off". Mitigation: assembly returns the active term list alongside H, and the artifact records it; the thesis's own `H_K = 0` deperturbation study (digest §6) is a legitimate use, so this must be visible, not forbidden. (iii) A block whose active terms violate its own blocking symmetry (see §3.1).

### 3.4 Parameter model

**Decision.** `Param(value, unit, uncertainty=None, source=None, status='unspecified', isotopologue=None, convention=None, note=None)` with a one-argument constructor `Param.of(7274.1)` so a scratch calculation is not a form-filling exercise. A `ParamSet` is a mapping symbol → Param plus a **conventions block** and a formalism tag.

**Units live in the data model.** Canonical internal unit MHz (both repos already); conversions at the boundary from the constants both repos already carry — cm⁻¹ × 29979.2458, D × 0.503412 → MHz/(V/cm), μ_B = 1.399624494 MHz/G, μ_N = 7.62259323e-4 MHz/G (`molecule_parameters.py:24–31`). The C2V `d_0` incident — a reduced-vs-physical dipole default living in `scan2k.make_params` while the physical override lived in a different script, silently scaling every answer (digest wart #3) — is exactly what a unit field on the value prevents.

**The conventions block is not bureaucracy; the ThF⁺ literature has five live forks.** Each is a documented factor-of-2 or sign disagreement between sources Arian will read side by side:
| Fork | Options | Source |
|---|---|---|
| `n_hat` | F→Th (JILA) vs Th→F (Petrov) | lit §5.2 — flips Ω, signed d, and E_eff's sign |
| `dipole_origin` | centre of mass vs Th nucleus | lit §5.4 — a 0.72 D / ~20 % Stark error |
| `dg_def` | δg = (g^u−g^ℓ)/2 (Ng thesis) vs Δg = g^u−g^ℓ (Petrov) | lit §5.10, confirmed numerically consistent at 0.5 % once applied, §7.3 |
| `edm_factor` | −d_e E_eff Ω/(2\|Ω\|) (Leanhardt) vs Ω/\|Ω\| (Ng, modern) | lit §5.11 |
| `zeeman_sign` | E = −g μ_B B m_F (Ng, Petrov) vs Leanhardt's absorbed sign | lit §5.8 |
Plus the two already in Molecule-Structure: `formalism ∈ {R2, N2}` (with `Lambda` mandatory — `formalism.py:9–13`, "Λ is never inferred") and `dtype_convention ∈ {BC, PGopher}` for the d-type sign flip. **The PGopher flip is applied by a converter and recorded, never hand-entered as a negative number with a comment** — the current practice (`molecule_parameters.py:64, 96, 161`) is a documented trap in the repo's own README ("Do not 'fix' a minus sign … without checking the inline source comment").

**Status** ∈ {measured, ab-initio, derived, estimate, held-fixed, stale, unspecified}. It exists because Arian already needs it: `RaF.boson.A0 'ASO': 1350` is documented-stale (`CLAUDE.md`, digest §7), and for ThF⁺ the sign of g_F is unmeasured, G∥ spans 0.034–0.048 across four sources, and E_eff spans 35.2–37.3 GV/cm (lit §7.5, §7.7). The report layer surfaces statuses; nothing gates on them (per `feedback_exploratory_estimates_no_gates` — estimator parameters are ordinary molecule parameters).

**Sweeps** are `(base_ParamSet, {symbol: array})`, not a list of full sets, so the coefficient matrix `c[n_sets, n_terms]` for §3.3 is built by broadcasting rather than by re-reading records.

*What could go wrong.* (i) Double conversion — so conversion is a `ParamSet` method that stamps the tag and is idempotent by tag check; keep `formalism.py` verbatim (118 lines, exactly invertible, 22 kernel tests, digest §6 as-is #2). (ii) A record with seven optional fields becomes friction and people bypass it with raw floats — hence `Param.of`, and the assembler accepting bare floats with a one-time warning. (iii) A `convention` mismatch between two Params in the same set is a real error and should raise (structure, not hash).

### 3.5 Engine

**Lift from C2V, rename, and add the one thing it lacks.** Verified this session: `qgt.py:180 multi_curvature(H0, verts, *, eigh=np.linalg.eigh)` — already the general n-knob form; `repair.py:827 _stream_diag_core` — the batched streaming eigh kernel; `scan_1d.py:73 batch_diagonalize_scan`; `scan_1d.py:11 track_state_ordering`; `matching/_utils.py:453 _PIN_TIE = 1e-6` — the tie tolerance that makes CPU and CUDA agree. The C2V digest verified the whole engine runs on random 8×8 matrices with `qnum=None`.

**Renames and surgery (from C2V digest §6, "reusable with surgery" 1–3):** delete the dead `qnum` parameter (11 lines, all signature/docstring); replace `stark`/`zeeman` positional slots with a `{knob: matrix}` dict; replace the hard-wired `obs_ops = [stark, zeeman, EDM_mat]` (verified at `scan_2d.py:293` and `scan_windowed.py:120`) with a caller-supplied `{name: matrix}`; replace `hamiltonians_from_params` (`repair.py:653`, which hard-codes `generate_states_and_blocks([1],[1],…)`) with the injected builder `build(spec, params) -> (H0, {knob: M}, {obs: M})`.

**Add the missing 1D entry point that returns eigenvectors.** C2V has none (digest Q3, explicit); `batch_diagonalize_scan` returns `(evals, expects)` and discards vectors. v1 needs vectors for TDMs, parity labels and Σ expectation values, so this is not optional.

**Backends: numpy/scipy default, torch optional and lazily imported.** At dim 96 over 10 m_F blocks, GPU is pure overhead. The torch path is lifted verbatim from C2V's `config.py:96 _get_eigh` when a polyatomic basis needs it. Keep `'syevj'` out — 30 % slower there and it once shipped transposed eigenvectors (C2V wart #6).

**Three orthogonal axes, per Rule 2** — nothing here is fundamentally incompatible, so nothing is forbidden:
- `order ∈ {energy, adiabatic_step, adiabatic_zero_field}`. Energy order is exact by construction and point-local; it is C2V's production default and the right v1 default. `adiabatic_zero_field` is the thesis's stated rule (p.267: order by "adiabatically correlated free field state") and **neither repo implements it** — Molecule-Structure chains step-to-step overlap (`Energy_Levels.py:1344`), C2V repairs an energy-order seed. It is a different algorithm, and it is the one the thesis's own labels assume.
- `gauge ∈ {none, pinned, matched}` — lift C2V's `apply_sign` family; `pinned` is the reproducibility-preserving one.
- `assignment ∈ {argmax, hungarian, adaptive}` — lift `_assign`; note Molecule-Structure's greedy `np.argmax(overlap, axis=1)` is not guaranteed to be a permutation and can silently duplicate a trace at a near-degeneracy (digest §2f), and that `scipy.optimize.linear_sum_assignment` is already a dependency there and unused.

**The label convention travels in the artifact.** C2V's hardest lesson (wart #9: energy-order and character labels disagree at ~3.7 % of magic nodes, and repairing a dataset invalidates a catalog built on the other convention). Every result object carries `order=`, `gauge=`, and the active term list.

*What could go wrong.* (i) Adiabatic tracking cost — overlap is O(d²) per step and Hungarian is O(d³) worst case; at d ~ 10³ over a 10⁴-point grid that dominates everything else. It is why energy order is the default and tracking is opt-in. (ii) `np.sign(0)` zeroing an eigenvector in a naive port of the phase fix (`Energy_Levels.py:312`); C2V's `_PIN_TIE` rule is the fix, take it rather than re-deriving. (iii) The thesis rule fails through a genuine avoided crossing where the adiabatic state changes character — thesis open question 2, and a real one for a K-resonance-like manifold.

### 3.6 Observables

- **Eigenvalues, eigenvectors, expectation values.** One einsum with a caller-supplied subscript (C2V `observables.py`, 29 lines, take verbatim).
- **Derivatives, all exact, no finite differences.** `multi_curvature` gives first derivatives (Hellmann–Feynman) and exact sum-over-states Hessians for any named knob subset. So a g-factor is `−(1/μ_B)∂E/∂B_z`, a dipole moment is `−∂E/∂E_z`, a polarizability is `∂²E/∂E_z²`. Molecule-Structure computes both by finite difference today (`g_eff_evecs:408`, `D_eff_evecs:422`, digest §2g); replacing that with the analytic kernel removes a step-size parameter and a re-matching call per point.
- **The new capability the architecture hands over for free: ∂E/∂p_k for any Hamiltonian parameter**, because `∂H/∂p_k = M_k` is already in the catalogue. That is the analytic Jacobian for a Levenberg–Marquardt fit — the thesis's own fitting procedure (Nelder–Mead → LM, digest §5). No optimiser ships in v1, but the design leaves the door open at zero cost, and closing it would require the fields/parameters split rejected in §2.1(2).
- **Pair differentials are first-class.** The two headline ThF⁺ observables are differences: `g_eff = g_S μ_B(⟨M_S⟩_{M=+1} − ⟨M_S⟩_{M=−1})` (thesis Fig. 5.28) and `Δg = g^u − g^ℓ` between Ω-doublets (lit §5.10). C2V has the machinery (`leg_track` keeps leg order so pair differences stay signed, and the `f ≥ 0` fold once turned a zero-g Newton objective into a V-minimum with no sign change — digest Q4). A pair-observable API that names the two partners and returns a signed difference is required, not a convenience.
- **TDMs between two separately diagonalised manifolds**, following the thesis procedure exactly (digest §5): diagonalise ground and excited separately → convert both eigenvector sets to case (a) → build the TDM matrix in the common basis → `G_evecs @ TDM @ E_evecs.T` → line positions from eigenvalue differences. **Amplitudes are summed and then squared** (p.161, so intensity-borrowing paths interfere); this is a one-line ordering choice with a factor-of-anything consequence, so it is encoded once and gated with a two-path interference toy (§3.7). Port `RaX/gen_spectra.py` — molecule-agnostic, two import rewires (digest §6 as-is #5).
- **PT-odd.** Extra term matrices (W_d d_e, W_T,P k_T,P) plus expectation values. One trap to carry forward: an imaginary-Hermitian operator has identically zero diagonal expectation value on real eigenvectors — Molecule-Structure's `PTV_shift` documents this at `Energy_Levels.py:481–483` and directs the reader to the off-diagonal element. So the observable API needs a named-pair off-diagonal mode alongside the diagonal one.

*What could go wrong.* Convention leakage into observables: whether `Δg` or `δg` comes out depends on the ParamSet's `dg_def` field, and whether the eEDM shift carries Leanhardt's ½ depends on `edm_factor`. Both must be read from the conventions block and stamped on the result, or two people comparing numbers will differ by exactly 2 and neither will know why (lit §7.3, §7.4).

### 3.7 Verification strategy

Three tiers, no snapshot tier. **A** = Hamiltonian-agnostic kernel; **B** = physics identity or limit; **D** = published-data comparison, labelled and opt-in via env flag. Every gate names the failure mode it uniquely catches, and both PASS and FAIL must be reachable — demonstrated by a fabricated failing input, not asserted.

**Tier A (first tests, in build order).**
| Gate | Check | Uniquely catches |
|---|---|---|
| A1 | `H(c) == Σ c_k M_k` on random `M_k`, random c | assembler indexing/ordering drift (generalises Molecule-Structure GATE A, `export_crossing_subspace.py:323`) |
| A2 | batched `tensordot` resum == per-point loop, random matrices | the sweep fast path diverging from the reference path |
| A3 | `multi_curvature` vs central finite difference on random Hermitian; role-swap; common-mode-vertex zero | derivative-kernel sign/self-term errors (port `test_multi_axis.py`, `test_trap_cross.py` — the template) |
| A4 | assignment permutation-ness; sign-gauge idempotence; carry is ints not vectors | tracking corruption (port `test_carry.py`, `test_gauge*.py`) |
| **A5** | **every term: non-zero somewhere inside its declared rules, exactly zero outside** | **dead operators. Positive control: `LambdaDoubling_p2q` → 20 non-zero. Negative control: `LambdaDoubling_q_even_aBJ` → 0 non-zero over 136 elements (digest §3). Both outcomes proven reachable on real code.** |
| A6 | fast Wigner backend vs sympy exact on a random argument sample including half-integers | backend precision and half-integer argument-convention mismatch |
| A7 | basis invariants: every value a multiple of 0.5; \|Ω\| ≤ J; triangle conditions; blocking is a partition (union = basis, pairwise disjoint) | enumerator bugs, which are otherwise invisible until a spectrum looks odd |
| A8 | `import <pkg>` under 1 s and imports no matplotlib/pandas/sklearn | the exact regression C2V has (6.5 s import from dead deps, `config.py:8–20`) |

**Tier B (physics identities and limits).** B1 Hermiticity of every assembled block, and realness where declared. B2 `[H₀, P] = 0` for a parity-conserving term set; `P² = 1`; parity eigenvalues ±1. B3 Kramers degeneracy at B = 0 for half-integer F. B4 **case (a)↔(b) cross-check**: the thesis's typo-corrected transform (Eq. 2.30) is unitary, and H assembled in (a) then transformed equals H assembled in (b) for a shared term set — non-tautological, because the two element sets are independently written; this is the strongest single test of the convention module. B5 S-uncoupling limit: as A_SO/B → ∞ the Ω = 1 admixture of Ω = 2, 3 scales as B/A_SO. B6 e/f ordering matches `P = ±(−1)^{J−S−ℓ}` (thesis item 21). B7 TDM sum rule (port Molecule-Structure `test_sio_stark_tdm.py:130`, m- and F-independent, closed form). B8 sum-then-square: a two-path interference toy where squaring first gives a different answer.

**Tier B, ThF⁺-specific closed forms** — all from the literature digest, all with both outcomes reachable:
- hyperfine eigenvalue `A∥[F(F+1) − I(I+1) − J(J+1)]/(2J(J+1))` ⇒ J = 1 splitting = ¾|A∥|, J = 2 = 5/12|A∥|, **ratio 9/5** (lit §3.1, §7.1);
- Ω-doublet splitting `ω_ef·J(J+1)/2` with parity alternating as `(−1)^J` — negative parity upper at J = 1, positive upper at J = 2 (lit §5.5, §5.6);
- fully polarised Stark slope `E = −m_F Ω γ_F d_mf E` with `γ_F = [J(J+1)+F(F+1)−I(I+1)]/(2F(F+1)J(J+1))`, γ_{3/2} = 1/3 (Leanhardt Eqs. 20–21, lit §3.2);
- `g_F = −G∥[F(F+1)+J(J+1)−3/4]/(2F(F+1)J(J+1)) + g_N(μ_N/μ_B)[F(F+1)−J(J+1)+3/4]/(2F(F+1))`, → (1/3)(−G∥ + g_Nμ_N/μ_B) at F = 3/2, J = 1 (Ng Eq. C.6; identical in Petrov Eqs. 2–3 — an independent restatement, so this is a genuine cross-source check);
- `δg_F/g_F = 9 d_mf E_rot/(40 B_e)` (Leanhardt Eq. 67, lit §3.2);
- avoided crossing `Δ ≈ 170 ω_ef (ω_rot/d_mf E_rot)³` for m_F = ±3/2 (Leanhardt Eqs. 42–43). The literature digest §8.6 flags reproducing this as a strong check; classify it B with a stated wide tolerance and carry Molecule-Structure's written rule verbatim — **constants are never tuned to make a gate pass**.

**Tier D, opt-in behind an env flag, labelled in the docstring.** Ng 2022 Table I constants; the 4B = 29.09733(4) GHz J = 1→2 interval; measured |g_{3/2}| = 0.0149(3); Petrov's Δg₀, Δg₁ table. Follow C2V's pattern (`skipif ATM_RUN_SLOW_PHYSICS`, `ATM_DGM_XCHECK`) — never default-on, because datasets and parameter sets churn.

**Ports as-is.** From C2V: `test_trap_cross.py` (the template, copy its docstring posture), `test_multi_axis.py`, `test_hessian.py`, `test_energy_curvature.py`, `test_carry.py`, `test_gauge*.py`, `test_sign_modes.py`, `test_scan_io.py`. From Molecule-Structure: `test_formalism.py`'s 22 kernel checks **minus** `test_raf_a0_regression_oracle` (class C); `test_sio_conventions.py` C1 (the exact 4×4 operator identity `b_F I·S + (c/3)√6 T²₀ ≡ b I·S + c I_zS_z`) and C3b (first-principles decoupled-basis reconstruction to 2e-16); `export_crossing_subspace.py` GATE A / GATE C / GATE E.

**Explicitly not gated.** Any ThF⁺ spectrum for a given parameter set; eigenvalues at pinned parameters; byte-identical output; anything hash-based. A cache-manifest fingerprint mismatch warns and rebuilds. Left behind from Molecule-Structure: `test_formalism_physics.py`, `test_raf_a0_regression_oracle`, `test_sio_crossing` G4, `test_sio_stark_tdm` G2, `sio_codegeneracy_map` G2 (all class C, digest §5). Left behind from C2V: `test_table3_coeffs.py` (66 KB of manuscript-structure assertions), `test_fig_env_overrides.py`, `test_trap_compare_loader.py` (unconditional read of an 800 MB production tree).

*What could go wrong.* The tier-B ThF⁺ closed forms come from sources that disagree with each other by exactly the factors in §3.4's conventions table. A gate written against Leanhardt and a code written against Ng will differ by 2 and the gate will look broken. Every tier-B test must state which convention block it assumes.

### 3.8 Package layout, dependencies, environment, entry points

Twelve modules, flat, no subpackages.

| Module | One line |
|---|---|
| `spec.py` | `StateSpec`, ket dtype, enumerator registry, `Blocking` and `merge` |
| `conventions.py` | the 21-item contract; parity and time-reversal operators; case (a)↔(b) transform |
| `wigner.py` | 3j/6j/9j behind one interface: cached sympy default, optional fast backend, sympy as the test oracle |
| `terms.py` | `@term` registry, `Rules`, citation metadata, per-case lookup |
| `elements_c.py` | case (c) matrix elements (v1) |
| `elements_a.py` | case (a) matrix elements, lifted from Molecule-Structure |
| `elements_b.py` | case (b) βJ / βS matrix elements, lifted |
| `params.py` | `Param`, `ParamSet`, conventions block, units, the R²↔N² converter (`formalism.py` verbatim) |
| `assemble.py` | `TermMatrices`, cache key and manifest, `H(params, knobs)`, batched sweep resum |
| `engine.py` | eigh backends, 1D/2D sweeps returning eigenvectors, batching |
| `track.py` | assignment, sign gauges, ordering policies, label recording |
| `observe.py` | expectation values, `multi_curvature` port, pair differentials, PT-odd off-diagonal mode |
| `spectra.py` | TDMs between manifolds, line strengths, branching ratios (`gen_spectra` port) |

That is thirteen; `spectra.py` arrives at M5, so v1 ships twelve. **Deferred with a named trigger:** `io.py` (persistence, ported from C2V `scan_io.py`) lands when one sweep exceeds memory — at dim 96 it never will.

**Dependencies.** Required: numpy, scipy, sympy, pytest. Optional: torch (lazy, GPU path), a fast Wigner backend. **Verified this session:** `wigners` (Rust) ships `win_amd64` and `macosx_11_0_arm64` wheels for every release **but exposes only 3j, Clebsch–Gordan and Wigner-D arrays — no 6j, no 9j** (read from its `python/wigners/__init__.py` ctypes bindings); `py3nj` has 3j/6j/9j but publishes **sdist only, no wheels at all** (seven files on PyPI, all `.tar.gz`), so it needs a Fortran toolchain on both machines. Neither is a drop-in. **Decision: sympy-behind-`lru_cache` is the default backend** — it is what both repos already use (`physics.py:12–27` records a 99.9 % hit rate and ~2.2× faster builds; `matrix_elements.py:18–23`), it is exact, and the matrix-first architecture pays the Wigner cost once per basis rather than once per parameter set, which removes the pressure that would justify a fragile dependency. `wigner.py` is an interface so a fast 3j path can be swapped in behind gate A6 if profiling ever demands it.

**Environment.** New conda env, Windows 11 + macOS, Python 3.12 (matching Molecule-Structure's pin). Keep the import light — gate A8 above.

**Entry points.** A functional API is the primary surface: `spec → blocks → term_matrices → assemble → sweep → observables`, each stage a plain function on plain data. One convenience object (`Manifold`) wraps the common path. No CLI, no notebook boilerplate, no `sys.path` injection (`config_path.py` in Molecule-Structure exists only because `Source Code/` is not a package — install this one with `pip install -e .`).

*What could go wrong.* Twelve modules is a lot for one physicist to hold; the mitigation is that seven of them (`wigner`, `elements_*`, `engine`, `track`, `observe`) are near-verbatim lifts with known-good tests, so the genuinely new surface is `spec`, `terms`, `assemble`, `params`, `conventions`.

### 3.9 Migration plan from the two repos

**Lift verbatim (rename only).** From C2V: `qgt.py` (343 lines) — the exact-derivative family, self-term-excluded-by-index and the zero-coupling-at-zero-gap guard are both hard-won; `observables.py` (29 lines); `matching/_utils.py` — `_assign`, cost kernels, `signed_boundary_overlap`, `Carry`, the sign-gauge family and `_PIN_TIE = 1e-6`; `repair.py:120 _batched_diag`, `:158 _diag_rows`, `:827 _stream_diag_core`; `config.py`'s `_get_gpu` / `_get_eigh` **minus lines 8–20** (dead pandas/sklearn/matplotlib/sympy imports); `holonomy.py` and `ci_finder.py` if degeneracy exploration is wanted; `bundle.py`'s `__main__` self-check pattern; `test_trap_cross.py` as the gate template. From Molecule-Structure: `formalism.py` (118 lines, verbatim); the `matrix_elements.py` case (a) and case (b) families; the `quantum_numbers.py` enumerators; the `record(name, passed, detail)` + exit-code gate style and the "constants are never tuned" rule; the inline citation discipline in `molecule_parameters.py` (the data and its provenance comments are hand-curated and irreplaceable).

**Lift with surgery, in this order.**
1. `quantum_numbers.py` enumerators → thread S and Λ from the spec; delete `K_mag` as a dispatch constant. Unblocks the case (a) ³Δ₁ basis immediately (already verified working when called directly).
2. `matrix_elements.py` → strip every default from every signature (§3.2); audit each body for a literal `1/2` or `5/2` — known: `H_odd_X:224` and `build_PTV_bBS:295–301` bake I = 5/2, `build_PTV_NSDPV:339` bakes I = 1/2.
3. The four `H_*` builders → **delete**, replaced by the table-driven assembler. They are four near-copies of one double loop; this also removes the unguarded `params[...]` reads at `hamiltonian_builders.py:50, 154–156, 223, 246` by construction.
4. Field layer → `{knob: matrix}`. `StarkX_bBJ` (`matrix_elements.py:258`) and `ZeemanX_bBJ` (`:196`) already exist and are gated but have no builder hook; the reason they were never retrofitted (it would change `H_function`'s signature and break the tutorial regression, per `test_sio_stark_tdm.py`'s docstring) does not apply in a new package.
5. `state_ordering` → C2V's `_assign` (`Energy_Levels.py:1354`'s greedy argmax is not guaranteed to be a permutation).
6. `molecule_parameters.py` → keep the numbers and citations, drop the `[species][boson|fermion][label]` nesting.
7. PT-odd operators → thread I; keep the physics.
8. C2V `orchestrate.py` → keep the flag validation and the compatibility rules; replace the three C2v lines with the injected builder.

**Leave behind.** `molecule_library_class.py` in its entirety (647 lines, nine `iso_state` keys × fourteen dispatch dicts, ~120 lines of commented-out alternates — every generality problem traces here); the `iso_state` string masquerade (`Energy_Levels.py:110–130`, where boson RaF X is literally `'174X000'`); `fermion_or_boson` as an architectural axis (it fuses "how many nuclear spins" with "how are they coupled"); `H_odd_A` (no field terms, and `symbolic=False` raises `NameError` at `:274`); `matrix_elements_sym.py` (binds sympy then numpy to the same name at lines 1–2); `case_bBJ_MEs.py`; `update_params` (full O(dim²) rebuild — the exact anti-pattern this package exists to remove); `CPT_Sims.py` (qutip is commented out of the env). From C2V: `physics.py` lines 55–700 (basis + all eleven MEs, all C₂ᵥ), `analysis.py`, `magic_geometry.py`, `polarizability.py`, all 20 `fig_*.py`, `v1_reference.py`, the manuscript test machinery.

**Bug list carried forward, not fixed in place** (both repos stay read-only): the eight items in synthesis §1, of which `LambdaDoubling_q_even_aBJ` is the one with a downstream consequence — three molecules carry a non-zero `q` that does nothing. Gate A5 in the new package catches it on arrival.

*What could go wrong.* Lifting `matrix_elements.py` also lifts nine `#check` / `#Check derivation` comments (`:276, 289, 337, 397, 405, 424, 447, 930, 962, 974`), several attributing a form to "Yuiki" with a noted index swap. Those are unaudited physics travelling under the appearance of working code. Gate A5 catches dead ones; it does not catch a wrong-but-non-zero phase. §4 Q5 asks Arian which he has since audited.

### 3.10 Milestones

Each ends in a runnable check. The ThF⁺ slice is M4, not the end.

- **M1 — spine, one to two days.** `wigner.py` + `spec.py` (case (c) enumerator) + `terms.py` + three terms (rotation `B J(J+1)`, hyperfine A∥, Ω-doubling ω_ef) + `assemble.py` + numpy `eigh`. **Check:** A1, A2, A6, A7, B1, plus the closed-form hyperfine 9/5 ratio and the ω_ef·J(J+1)/2 scaling with `(−1)^J` parity alternation. Reachable in a day because case (c) needs no Σ, no Λ and no 9j, and the whole J = 1–4 basis is 96 kets.
- **M2 — parameters and fields.** `params.py` with the conventions block; Stark and Zeeman as knobs; per-m_F blocking. **Check:** γ_F = 1/3 full-polarisation Stark slope; the g_F-from-G∥ formula against Petrov's independent restatement; B2 parity; B3 Kramers; A5 across the whole registry.
- **M3 — engine and derivatives.** Batched sweeps returning eigenvectors; `multi_curvature` port; pair differentials. **Check:** A3 against finite difference; `δg_F/g_F = 9 d_mf E_rot/(40 B_e)`; A4 tracking invariants.
- **M4 — the ThF⁺ vertical slice.** 2D (E_z, B_z) map per m_F block over a parameter-set collection; Δg between Ω-doublets; PT-odd matrix elements. **Check:** the tier-D opt-in comparison to Ng 2022 Table I and to Petrov's Δg₀/Δg₁, both labelled. Note in advance (lit §7.8): a single-electronic-state model is expected to miss δg/g by ~15 %; that is a *finding*, not a failure, and adding ³Δ₂ to `electronic` is the follow-up.
- **M5 — spectra.** `gen_spectra.py` port; TDMs within X. **Check:** B7 sum rule, B8 sum-then-square.
- **M6 — case (a) subengine.** Lift `matrix_elements.py` case (a) and `q_numbers_even_aBJ`; build the ³Δ manifold. **Check:** B4 (case (a) vs case (b) cross-basis agreement) and B5 (S-uncoupling limit). This is the milestone that actually validates the convention module, and it is also where the case (a) ³Δ₁ Ω-doubling operator (§4 Q2) and the Σ ≠ ±S parity phase (§4 Q3) must be settled.
- **M7 — linear polyatomic ket.** Extend the ket dtype with (v₂, ℓ, K); reuse everything else. **Check:** the ℓ-doubling and Renner–Teller term forms against the thesis Table 4.2 structure.

---

## 4. Recorded scientific and design uncertainties

**Physics.**

**Q1. Which basis for the v1 slice?** (a) Case (c) |J, Ω=±1, F, m_F⟩ with phenomenological ω_ef — the Ng/JILA formulation, 96 kets for J = 1–4, matches every published ThF⁺ number, buildable in M1. (b) Case (c) with several electronic states (³Δ₁, ¹Σ⁺, 1³Δ₂, ³Π₀±) and off-diagonal couplings — the Petrov formulation; required to reproduce the measured δg/g (lit §7.8), costs an inter-state term category and ab initio off-diagonal matrix elements you do not control. (c) Case (a) with S = 1, Λ = 2, full ³Δ manifold — connects to the lifted matrix-element library and to the polyatomic roadmap, but needs the ΔΛ = ±4 Ω-doubling operator (Q2) and the Σ ≠ ±S parity phase (Q3) settled first. *Cost of choosing (a) now:* case (a) is deferred to M6, which is where it belongs anyway. *Cost of choosing (c) now:* M1 slips and the first number is blocked behind two open physics items.

**Q2. Ω-doubling operator for a ³Δ.** (a) Phenomenological, case (c): `(−1)^J (ω_ef/2)[J(J+1)/2](|+1⟩⟨−1| + h.c.)`, ω_ef = 5.29(5) MHz measured (Ng Eq. C.3). (b) Microscopic, case (a): `½(o_Δ + 3p_Δ + 6q_Δ)(S₊²J₊² + S₋²J₋²)` with `ω_ef = 4õ_Δ` (Leanhardt Eq. 16 + p.14). (c) Generated, not fitted: Petrov's route, where the doubling comes out of ⟨³Δ₁|L̂₊ − g_SŜ₊|Ω=0⟩ matrix elements tuned to reproduce the measured 5.29 MHz. *Costs:* (a) is free and exact for v1 but carries no J-dependence beyond J(J+1) and no centrifugal correction (none is published — lit §8.4); (b) is the only form that generalises to other Δ states and needs a derivation and a source; (c) needs electronic structure inputs. Note the three constants ω_ef, Gresh's k, and õ_Δ differ by factors 2 and 4 (lit §5.5) — whichever you pick, the conversion goes in the parameter record, not in a comment.

**Q3. Does the composite parity phase `(−1)^{J−S−ℓ+s}` survive Σ ≠ ±S?** The thesis derives Eq. A.15 using `S = |Σ|` (thesis digest §10 Q6), and a ³Δ manifold has Σ = 0 components exactly where that step was taken. Options: (a) you have already checked it and it carries over; (b) it needs re-derivation before any case (a) ³Δ parity label is trusted; (c) sidestep by working in case (c) for v1, where parity comes from the Ω-doublet superposition. *Cost of not settling it:* every case (a) triplet parity label is unverified, and B6 (e/f ordering) has no reference.

**Q4. n̂ direction, and therefore the sign of everything.** JILA points n̂ from F to Th (Leanhardt p.15, Ng p.318); Petrov/Skripnikov point it from Th to F (arXiv:2503.02840 p.3). Ω, the signed dipole, and the sign convention for E_eff all flip. *Cost:* pick one and assert it in `conventions.py`, with a converter for the other — or every cross-check against the other lineage silently sign-flips.

**Q5. Which of the flagged matrix elements do you trust?** Nine functions in `matrix_elements.py` carry `#check`, `#check phase factor`, `#Check derivation` or `#Check sign`, several attributing the form to "Yuiki" with a noted index swap. Options: (a) name the audited ones and lift only those; (b) lift all and re-derive on use; (c) lift all and rely on B4 (case (a)↔(b) agreement) to catch the ones that matter. *Cost of (c):* B4 only catches terms that exist in both cases.

**Q6. Field geometry in v1.** (a) Collinear E_z, B_z only, with E_x/B_x/ω_rot declared but unwired. (b) Rotating-frame `ħω_rot F_x` wired from M2 — it is one more knob, but it breaks m_F conservation, so the default blocking becomes full-M and the v1 dimension goes from ~10 per block to 96. At this size that is free; at polyatomic sizes it is not. *Cost of deferring:* the Berry-phase and avoided-crossing physics (Leanhardt Eqs. 36–43) is the JILA experiment's actual operating regime.

**Q7. Hyperfine scope.** The ThF⁺ literature gives only A∥; no dipolar, contact or e_Δ constants exist for ThF⁺ anywhere (lit §8.2), and no rotational g-factor (lit §8.3, estimated at ~6 % of the total in ThO). Options: (a) A∥ only, matching every published model; (b) add the full Frosch–Foley set with estimated values marked `status='estimate'`; (c) add them as registered terms with zero default so they can be switched on later. *Recommendation implied by the term registry:* (c) is nearly free — but it is your call whether estimated values ship with them.

**Numerical and API choices recorded at the design stage.**

**Q8. Ordering policy default.** Energy order (exact, point-local, C2V production default) vs the thesis's adiabatic correlation to the zero-field state (p.267 — which neither repo implements as stated). *Cost:* the thesis rule is the one your published labels assume; energy order is the one that cannot silently mislabel.

**Q9. Precision floor.** Molecule-Structure rounds eigenvalues and eigenvectors (`Energy_Levels.py:1375 diagonalize(..., round=10)`, with the state object's own default at 6 → a 1 µHz floor in MHz units). For a sub-Hz PT-odd splitting, keep no rounding at all?

## 5. Risks and unknowns, ranked

1. **Convention transplant.** The lifted operator library's phase conventions are documented in code comments, carry nine unaudited `#check` flags, and contain one proven sign bug that made an operator identically zero for years. Every number the new package produces inherits them. *Mitigation:* one `conventions` module; gate A5 on arrival; gate B4 (case (a) built and transformed vs case (b) built directly) as the strongest available cross-check. *Residual:* B4 cannot see a wrong phase that is common to both element sets.
2. **The ³Δ₁ Ω-doubling operator for Λ = 2 in case (a).** Demoted from the synthesis's #1 because the case (c) form exists in closed form (§2.0), but unresolved for M6 and for any other Δ state. Molecule-Structure's ΔΛ = ±2 operators are identically zero on a Δ state (verified by exhaustive sweep with a Λ = 1 control); QuantumStates.jl's is also ΔΛ = ±2; the thesis points at Brown, Cheung & Merer 1987 but never writes it. *Mitigation:* registry entry, phenomenological alternative always available.
3. **Wigner backend.** Verified this session: `wigners` has wheels everywhere but no 6j/9j; `py3nj` has 6j/9j but no wheels at all, so it needs a Fortran toolchain on Windows and macOS both. The sympy-plus-cache default removes this from the critical path — but only because the matrix-first architecture pays the Wigner cost once per basis. If that assumption breaks (a polyatomic basis where even one build is minutes), the fast-backend problem returns with no good option, and the fallback is writing 6j/9j on top of a fast 3j.
4. **Multi-electronic-state requirement.** A single-state effective Hamiltonian misses the measured δg/g by 15 % (Ng thesis p.85 vs lit §7.8). Reproducing that observable requires the `electronic` list and ab initio off-diagonal electronic matrix elements (G∥ alone spans 0.034–0.048 across four sources, lit §7.7).
5. **Adiabatic tracking cost at larger bases.** O(d²) overlap per step plus up to O(d³) assignment, times the grid. At v1 dimensions invisible; at polyatomic full-M dimensions it dominates. C2V needed 967 lines of certified repair machinery to make it trustworthy on a 2D grid and still runs production on energy order. *Mitigation:* energy order default; tracking opt-in; do not port `repair.py` until something needs it.
6. **Over-engineering the parameter record.** Seven fields plus a conventions block is real friction if every scratch value must be wrapped. The failure mode is people bypassing it with raw floats, which is worse than not having it. *Mitigation:* `Param.of(x)`; bare floats accepted with a one-time warning. Watch this one — it is the place this design is most likely to be wrong in the direction of too much.
7. **Memory ceiling of the term catalogue** at polyatomic dimensions — O(n_terms · dim²). Sparse per-term storage buys roughly an order of magnitude; beyond that the architecture needs block-sparse or on-the-fly assembly for the largest terms, which would be a real change.
8. **The literature's own open items propagate into the package's defaults**: the sign of g_F is unmeasured (theory forces g_F < 0), no centrifugal correction to ω_ef exists, no rotational g-factor for ThF⁺ exists, and Petrov's printed body-fixed dipole appears to be a typo by a factor 10 (lit §7.2). These are `status` fields and switches, not bugs — but they mean no result is better than its inputs.
9. **Polyatomic scope.** The ket dtype and constraint machinery can represent ℓ/K/v before all corresponding operators exist. Each backend therefore states its implemented term set explicitly.
