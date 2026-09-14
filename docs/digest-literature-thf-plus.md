# Digest: ThF⁺ X ³Δ₁ literature — source coverage and absence record

This file is the source-coverage and absence record for the ThF⁺ X ³Δ₁ effective-Hamiltonian literature review conducted 2026-09-04. The physics — Hamiltonian terms, conventions, and parameter values — lives in [`thf-plus-x3delta1-effective-hamiltonian.md`](thf-plus-x3delta1-effective-hamiltonian.md), which re-read every sign, definition, and number at source; this file records only which sources were searched, which were read, and what was searched for and not found anywhere.

---

## Search scope and coverage

**Library search.** Title and attachment-name queries run on 2026-09-04:
  `*ThF*`, `*HfF*`, `*Gresh*`, `*Cairncross*`, `*Leanhardt*`, `*Cornell*`, `*Skripnikov*`, `*Petrov*`, `*Roussy*`, `*Loh*`, `*Zhou*`, `*Denis*`, `*Fleig*`, `*eEDM*`, `*electron electric dipole*`, `*thorium*`, `*Ng -*`, `*velocity modulation*`, `*trapped molecular*`, `*Omega-doubl*`, `*Grau*`, `*Meyer*`, `*Vutha*`, `*g factor*`, `*g-factor*`, `*spin-rotational*`, `*Titov*`, `*Mosyagin*`, `*thesis*`.
The broader `*thesis*` query found author-title attachments missed by author-name
queries. Exact-title searches covered ThF spectroscopy, velocity-modulation
spectroscopy, coherence, electric-field-dependent g factors, and the HfF⁺
spin-rotational Hamiltonian.

**Online search.** arXiv queries covered the exact ThF⁺ spectroscopy title,
ThF⁺ hyperfine and Ω-doubling, ThF g factors, Skripnikov/Petrov ThF⁺ work,
and Petrov's Zeeman papers. Hit counts and exclusions below preserve the search
boundary.

**Sources obtained as full text:**

| Source | Route | Note |
|---|---|---|
| Ng et al., PRA 105, 022823 (2022) | arXiv:2202.01346 (the JILA-hosted `PhysRevA.105.022823.pdf` URL in the brief returns an HTML error page, 34 kB, not a PDF) | read in full |
| Ng PhD thesis (JILA), 325 pp. | Zotero `QX2C5ZA4/KiaBoonNgThesis.pdf` | App. B (p. 316–317), App. C (p. 318–325), Ch. 2.4 (p. 35–40), Ch. 4.1 (p. 76–96) read in full |
| Gresh et al., J. Mol. Spectrosc. 319, 1 (2016) | arXiv:1509.03682 | read in full incl. Tables 1–3 (rendered) |
| Gresh PhD thesis (JILA), 145 pp. | Zotero `YM4P3CV6/dgresh_thesis.pdf` | TOC only |
| Leanhardt et al., J. Mol. Spectrosc. 270, 1 (2011) | arXiv:1008.2997 | Sec. III–IV read in full (pp. 12–25 of the arXiv PDF) |
| Petrov & Skripnikov, arXiv:2503.02840 (2025) | arXiv | read in full, pp. 3–4 rendered |
| Petrov, Skripnikov & Titov, arXiv:2302.02856 (PRA 107, 062814) | arXiv | read in full |
| Petrov, Skripnikov & Titov, PRA 96, 022508 (2017) | arXiv:1704.06631 | Theory section read in full |
| Skripnikov & Titov, PRA 91, 042504 (2015) | arXiv:1503.01001 | Tables I–II + PT-odd definitions read |
| Denis et al., New J. Phys. 17, 043005 (2015) | Zotero `KC5JC4TW` | abstract, Tables 5/10/11, conclusions |
| Cairncross et al., PRL 119, 153001 (2017) + Supplement | Zotero `E7RPWZE2` (18 pp., supplement included) | scanned; supplement is data collection/systematics, **not** the Hamiltonian |
| Cairncross PhD thesis (JILA), 260 pp. | Zotero `JKS9NAEE` | Ch. 2 pp. 36–50 extracted |
| Roussy et al., Science 381, 46 (2023) | arXiv:2212.11841 + Zotero `9VBXBDEZ` | skimmed only |
| Zhou et al., PRL 124, 053201 (2020) | Zotero `R553MJBH` | skimmed only |
| Loh, Grau, Stutz theses (JILA) | Zotero | copied, not read |

**Not obtained:** Loh et al., Science 342, 1220 (2013) and its supplement (paywalled; no arXiv version located — searched `au:Petrov AND ti:"Zeeman interaction"` and the ThF⁺ title queries above; the Loh thesis in Zotero is the substitute and was not read). Caldwell et al., PRA 108, 012804 (2023) (systematics; cited by Petrov 2025 as [51]). Petrov, PRA 108, 062804 (2023) (cited as [23] for the S₂/S₃ formulas). Skripnikov JCP 147, 021101 (2017) — this is HfF⁺ E_eff, not ThF⁺, so it was deliberately not pursued. Fleig & Nayak: the Zotero copy (`TCMGSTQ9`) is **ThO**, not ThF⁺; no Fleig ThF⁺ E_eff paper was located separately from Denis et al. 2015 (on which Fleig is a co-author). The Tier-2 item "Rotational splittings in diatomic molecules of interest to searches for new physics (2025)" was not pursued.

---

## Absence records

Quantities searched for and not found in any source read, with the `OPEN-*` item in [`open-questions.md`](open-questions.md) each one backs.

| Quantity | Sources searched | Consequence |
|---|---|---|
| Sign of g_{F=3/2} | Ng 2022 §II E p. 5, Ng thesis §4.1.4 p. 87, arXiv:2503.02840 p. 4 — all state the experiment is insensitive to it | `OPEN-4`: the signed g_F defaults to the theory-supported negative branch, not a measurement |
| ¹⁹F dipolar / contact hyperfine constants (a, b_F, c, e_Δ separately) | Ng 2022, Ng thesis App. C, Gresh 2016, arXiv:2503.02840 | `OPEN-5`: only the axial combination A∥ = 2a − b_F − (2/3)c is known; `OPEN-8`: e_Δ is omitted, so no measured hyperfine-dependent Ω-doubling exists to check against |
| Second-order Zeeman / rotational g-factor g_r for ThF⁺ | Ng 2022 §II D p. 5 says "It might be of interest to compute the rotational g-factor of ThF⁺" and notes it is ~6% of the total in ThO, citing Petrov et al. PRA 89, 062505 | `OPEN-7`: g_r is absorbed into G⊥ = (G_xx + G_yy)/2 of the Zeeman tensor, taken from Petrov's second-order sums rather than from a measurement (`docs/lit/lookup-effective-zeeman-tensor.md` §4.1) |

**`OPEN-9` (NSD-PV) is not backed by this digest.** The search scope above never queried nuclear-spin-dependent parity violation, W_a, or W_P — no digest source was read for that topic, so `OPEN-9`'s "not found in the sources reviewed" statement in `open-questions.md` rests on the reference document's own search, not on this digest.
