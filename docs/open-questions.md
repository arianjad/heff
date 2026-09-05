# Open questions raised during implementation

## OQ-A — the Ω-doubling sign and the parity of the upper doublet component

Raised: Task 5 of `docs/superpowers/plans/2026-09-05-heff-v1-thf-tutorial.md`.
**Status: resolved in code, pending Arian's confirmation of one statement.**

### The physical invariant (convention-free)

The **upper** Ω-doublet component has parity `(−1)^J` at every J — equivalently,
**e lies `ω_ef J(J+1)/2` above f, uniformly in J** under Brown's 1975 rule.

Two primary sources, read independently:

- Ng 2022 Fig. 2, rendered and read as an image ([HAM] §2.3). Caption: positive
  (negative) parity levels are black (grey). At **J = 1** the grey line is above
  the black in both F = 3/2 and F = 1/2; at **J = 2** the black is above the grey
  in both F = 5/2 and F = 3/2.
- Gresh 2016's `k″ < 0` for the X ³Δ₁ lower state, which is the same ordering
  ([HAM] §2.3, third overturn note).

### Why heff's element is `−ω_ef J(J+1)/4` and not Ng's `(−1)^J ω_ef J(J+1)/4`

heff's parity operator (`heff.conventions.parity_phase`, thesis Eq. A.15 = B&C
Eq. 6.234) acts as

```
E* |J, Ω⟩ = (−1)^{J−S−ℓ+s} |J, −Ω⟩ = (−1)^{J−1} |J, −Ω⟩     (³Δ: S = 1, ℓ = s = 0)
```

so the combination `(|+1⟩ + σ|−1⟩)/√2` has parity `σ·(−1)^{J−1}`. The eigenvalues
of an off-diagonal `c` in this two-state space are `±|c|`; the **upper** one
(`+|c|`) is the symmetric combination (σ = +1) when `c > 0` and the antisymmetric
one (σ = −1) when `c < 0`.

- With the **J-independent negative** element `c = −ω_ef J(J+1)/4 < 0`: σ = −1 at
  every J, so the upper parity is `−(−1)^{J−1} = (−1)^J` at every J. ✔ invariant.
  (J = 1: phase +1, antisymmetric has parity −1, and it is the upper one, +ω_ef/2.
  J = 2: phase −1, antisymmetric has parity +1, still upper.)
- With **Ng Eq. C.3 transcribed literally**, `c = (−1)^J ω_ef J(J+1)/4`: σ = (−1)^J,
  so the upper parity is `(−1)^J(−1)^{J−1} = −1` at **every** J. ✘ The two
  alternations cancel instead of composing; the figure's alternation is lost.

Ng's `(−1)^J` is the same physics written in a ket-phase convention where
`E*|J, Ω⟩ = |J, −Ω⟩` (no `(−1)^{J−S}`), i.e. his `|J, Ω = −1⟩` carries a
J-dependent phase relative to the ket heff's parity operator acts on. Nothing in
any source read for [HAM] writes that phase down, which is why this was flagged
rather than assumed.

The **splitting magnitude** `ω_ef J(J+1)/2` is identical under either form and is
unaffected by this question.

### Status in code

- `heff/elements_c.py:omega_doubling` returns `−J(J+1)/4` (coefficient `ω_ef`),
  with the transposition recorded in its `cite`.
- Both halves of [HAM] V4 are **hard gates** in
  `tests/test_elements_c_fieldfree.py`: the splitting law
  (`test_V4_omega_doubling_splitting_law`) and the parity ordering
  (`test_V4_upper_doublet_parity_alternates_as_minus_one_to_the_J`). The parity
  gate fails at even J on Ng's literal form and at odd J on a globally flipped
  sign, so both outcomes are reachable.
- No constant, tolerance or parity phase was adjusted.

### The question for Arian

> Confirm that the upper Ω-doublet component of ThF⁺ X ³Δ₁ has parity `(−1)^J` at
> every J (e above f uniformly, Brown 1975). If yes, the code is correct as
> written and this item closes. If the ordering is instead parity `−1` at every J,
> Ng Eq. C.3 is literal in heff's convention and the sign in `omega_doubling`
> flips back — one character, caught either way by the V4 parity gate.

A one-line question to K. B. Ng about the `|J, Ω = −1⟩` ket phase would settle the
convention side of it, and is the same message that settles OPEN-3.
