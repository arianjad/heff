# Checkpoint — two-photon notation vs the Pipi papers (2026-09-17)

Task: check `heff`'s effective two-photon (2 × E1) notation against Anastasia
Pipi's two arXiv papers and their appendices. Deliverable:
`docs/lit/2026-09-17-twophoton-notation-vs-pipi.md`. Read-only on `heff/`,
`tests/`, `scripts/`, `notebooks/`, `results/`.

## Milestone 1 — PDFs downloaded, two-photon sections located (DONE)

Downloaded to the session tmp dir (not committed):
`/Users/arianjadbabaie/.claude/jobs/260d4de9/tmp/pipi2024.pdf` (41 pp.,
11.4 MB) and `pipi2026.pdf` (19 pp., 2.3 MB). Layout text extracts left
**untracked** at `docs/lit/pipi2024-arXiv2410.11839-molecular-control.txt` and
`docs/lit/pipi2026-arXiv2608.03702-fno-control.txt`: `.gitignore` line 17 is
`docs/lit/**/*.txt` and the existing Petrov extracts are untracked for the same
reason, so the extracts were not force-added past that rule.

`pdftotext -layout` interleaves the two columns, so every quoted equation was
re-read with `pdftotext -f N -l N <pdf> -` (true reading order) before being
cited.

Locations found (all verified at source, not via subagent relay):

| paper | what | PDF page |
|---|---|---|
| 2410.11839v5 | Eq. (S4)/(S5), `H_int` with the `Omega/2` prefactor, Sec. SB | 11 |
| 2410.11839v5 | **Eq. (S10)**, the Raman-Rabi frequency, Sec. SC | 12 |
| 2410.11839v5 | no-RWA justification, `|w_iM - w1|` vs `|w_iM + w2|` | 12 |
| 2410.11839v5 | CaH+ reference Rabi rate `2pi x 2.087 kHz` | 12 |
| 2410.11839v5 | Table S2, per-pulse Rabi rates and `D = pi/(lambda_LD Omega)` | 18 |
| 2410.11839v5 | **Table S4**, H3O+ Rabi rates + selection rules in the caption | 28 (to 34) |
| 2410.11839v5 | **Fig. S11 caption**, `pi`(abs.) + `sigma-`(emit.) -> `dm = +1` | 39 |
| 2608.03702v2 | Fig. 1 caption, sideband detuning `D = nu_f + w_J'J - w` | 2 |
| 2608.03702v2 | Eq. (2) `H_int`; Eqs. (3), (4) resonance and detuning | 4 |
| 2608.03702v2 | **Appendix A 1**, Eqs. (A1)-(A10) — the defining block | 16 |
| 2608.03702v2 | Appendix A 2, Eqs. (A11)-(A13), block-diagonal + RWA | 17 |

Appendix/supplement structure: 2410.11839v5 main text pp. 1-9, supplement from
p. 10 (Sec. SA p. 10, SB p. 11, SC p. 12, SD p. 13, SE figures/tables p. 17+).
2608.03702v2 main text pp. 1-13, Appendix A p. 16, Appendix B p. 17,
Appendix C p. 18.

## Milestone 2 — zeros confirmed (DONE)

Neither paper contains any polarizability tensor, any Wigner-Eckart /
Clebsch-Gordan / 3j / Condon-Shortley statement, any rank-0/1/2 decomposition,
any polarization-vector component definition, any explicit field amplitude
(`E_0 cos wt`, `E_rms`, intensity), or any beam-propagation / helicity
statement. Verified with my own case-insensitive `grep -c -iF` over both
extracts against the control term `Rabi` (27 hits in 2024, 7 in 2026), so the
query shape reaches populated content.

## Milestone 3 — numerical checks run (DONE)

Scripts (tmp, not committed): `chk.py`, `chk2.py`; output `chk.out`,
`chk2.out`. Verbatim output is pasted into §4 of the deliverable.

- CHECK 1, `sigma+-` <-> `Delta m_F` for `(pi` absorbed, `sigma+-` emitted):
  heff's `reading='raman'` matches Pipi 2024 Fig. S11 and contradicts Pipi 2026
  Eq. (A10). Both PASS and FAIL reachable (the reading flip inverts the
  assignment).
- CHECK 2, `dJ` reach inside one parity eigenspace, 229ThF+ `J_max=3` (360
  kets), channel (K=2, dOmega=0): `|dJ| = 1` elements are 0.144-0.189 with
  `[T, Par] = 0` exactly. So heff's operator permits `dJ = +-1`; Pipi 2026
  Eq. (A4) forbids it. Traced to the molecule, not to heff.

## Milestone 4 — deliverable written and committed (DONE)

`docs/lit/2026-09-17-twophoton-notation-vs-pipi.md`.

## STOP items raised for Arian (not fixed here)

1. The two Pipi papers state **opposite** `sigma+- <-> Delta m_F` maps for the
   same `(pi` absorbed, `sigma+-` emitted) pair. Neither defines `sigma+-` in
   components and neither states the beam geometry, so the conflict cannot be
   resolved from the papers.
2. Pipi's Eq. (A5)/(S10) keeps both time orderings with unequal denominators,
   which retains a rank-1 (antisymmetric) two-photon part that heff excludes by
   construction (OPEN-21). Same physics as the existing OPEN-21 escalation, now
   with an external example.
