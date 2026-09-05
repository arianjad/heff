# Lookup: ground-state magnetic dipole moment of ²²⁷Th

Retrieved 2026-09-05. Question: does any compilation or paper give a measured/estimated
μ(²²⁷Th) and confirm the ENSDF-tentative (1/2⁺) ground state? Answer: **no** — three
independent nuclear-moment sources were queried and all three omit ²²⁷Th entirely.

## Sources reached

1. **N.J. Stone, "Table of Recommended Nuclear Magnetic Dipole Moments: Part I,
   Long-Lived States," INDC(NDS)-0794, IAEA Nuclear Data Section, November 2019.**
   PDF obtained from `https://nds.iaea.org/records/x2773-kj397/files/indc-nds-0794.pdf`
   (the plain `indc-nds-0794.pdf` URL 302-redirects to an InvenioRDM landing page,
   not the file itself — had to follow the record page to the real `/files/...` link).
   56 pages, extracted with PyMuPDF. Retrieved 2026-09-05.

2. **N.J. Stone, "Table of Nuclear Magnetic Dipole and Electric Quadrupole
   Moments"** (Oxford Physics, Clarendon Laboratory), literature search stated
   complete to **early 1998** — an earlier draft/edition than the well-known
   ADNDT 90, 75 (2005) print version, mirrored by PSI at
   `https://www.psi.ch/sites/default/files/import/low-energy-muons/DocumentsEN/nuclear-moments.pdf`.
   108 pages, extracted with PyMuPDF. Retrieved 2026-09-05. (Caveat: this is not
   verified to be byte-identical to the 2005 ADNDT paper or to the intermediate
   INDC(NDS)-0658 (2014) edition named in the task brief — I could not locate a
   directly downloadable copy of INDC(NDS)-0658 in this session — but it is a
   genuine Stone compilation predating INDC(NDS)-0794, useful as a second,
   independent cross-check.)

3. **IAEA NDS "Nuclear Electromagnetic Moments" live database** (Mertzimekis,
   distinct tool from the IAEA Live Chart `ground_states` endpoint the user
   already queried), `https://www-nds.iaea.org/nuclearmoments/isotope_measurement_results.php?A=227&Z=90`.
   Retrieved 2026-09-05 via curl (a WebFetch attempt on this URL returned
   "HTTP 402 Payment Required" — spurious; direct curl returned HTTP 200 with a
   normal page, so the WebFetch tool's proxy is the thing that choked, not the
   IAEA server).

4. Secondary literature: WebSearch queries "227Th magnetic moment",
   "thorium-227 nuclear magnetic dipole moment mu spin 1/2+", "odd-A thorium
   isotopes magnetic moments hyperfine anomaly laser spectroscopy", "Kälber
   thorium isotopes laser spectroscopy stored ions 227Th magnetic moment
   nuclear radii Zeitschrift Physik A", "227Th hyperfine structure magnetic
   moment ... nuclear". No PDF access (pdf-mcp and zotero MCP both down per
   this session's connection-failure notice; Springer paywalled the one
   directly relevant paper, redirecting to an IdP login).

## Verbatim Z = 90 rows

### INDC(NDS)-0794 (2019), page 42

Full Z = 89–91 context (verbatim, line-by-line as extracted):

```
89 Ac 227
0
21.77 y
3/2-
+1.1(1)/1.22(17)
O/R
1955Fr26/2017Gr18
PR 98 1514 (55)/PR C96 054331 (2017)

90 Th 229
0
7880 y
5/2+
+0.46(4)
O
1974Ge06
JPPa 35 483 (74)

91 Pa 228
0
22 h
3+
3.5(5)
NO/S
1989He07
NP A493 83 (89)
```

**The Z = 90 block contains exactly one row: ²²⁹Th. There is no ²²⁷Th row.**
The table goes directly from ²²⁷Ac to ²²⁹Th — A = 227 is skipped entirely at
Z = 90.

### PSI-hosted Stone table (literature complete to early 1998), page 149

```
89 Ac 227
0
21.77 y
3/2-
+1.1(1)
O
1955Fr26
PR 98 1514 (55)/PR 111 1747 (58)
+1.7(2)
O
1955Fr26
PR 98 1514 (55)/PR 111 1747 (58)

90 Th 229
0
7340 y
5/2+
+0.46(4)
[239Pu]
O
1974Ge06
JPPa 35 483 (74)
+4.3(9)
O
1974Ge06
JPPa 35 483 (74)

90 Th 232 gsband
g(18-24)>g(10-16)
TF
1992Ha03
PRL 48 383 (82)
g(av)=0.28(2)

91 Pa 228
0
22 h
(3+)
3.5(5)
NO/S
1989He07
NP A493 83 (89)
```

Same conclusion: Z = 90 has entries only for ²²⁹Th (ground state) and ²³²Th
(g-band), both bracketed directly by ²²⁷Ac below and ²²⁸Pa above. **No ²²⁷Th
row.** (Half-life for ²²⁹Th differs, 7340 y vs 7880 y in the 2019 table — an
old-vs-updated adopted half-life, doesn't affect the conclusion. The
magnetic-moment value +0.46(4) traces to the same 1974Ge06 reference in both
editions.)

### IAEA NDS live "Nuclear Electromagnetic Moments" database, query A=227&Z=90

The query resolves to the Thorium (Z=90) index page rather than an isotope
page, and that index lists exactly two isotope links:

```html
<td class="nuc2"><strong>Thorium (Z=90)</strong></td>
...
<td class="nuc"><a href="isotope_measurement_results.php?A=229&Z=90"><sup>229</sup>Th</a></td>
<td class="nuc"><a href="isotope_measurement_results.php?A=232&Z=90"><sup>232</sup>Th</a></td>
```

No ²²⁷Th link is present — this database, like the two Stone-table PDFs, has
no electromagnetic-moment record at all for ²²⁷Th.

## ²²⁹Th positive control (both PDFs)

- INDC(NDS)-0794: `90 Th 229, ground state, 7880 y, 5/2+, μ = +0.46(4) μ_N,
  method O, ref 1974Ge06 (J. Phys. (Paris) 35, 483 (1974))`.
- PSI-hosted table: `90 Th 229, ground state, 7340 y, 5/2+, μ = +0.46(4) μ_N
  [239Pu reference], method O, ref 1974Ge06`, plus a second line `+4.3(9), O,
  1974Ge06` (evidently the accompanying quadrupole-moment-type entry from the
  same 1974 measurement, not a magnetic moment — this table interleaves μ and
  Q entries under one isotope header).

Both values (0.46(4) μ_N) fall inside the 0.36–0.46 band named as the sanity
check, and the 1974Ge06 citation is identical in both editions — confirming
the extraction pipeline (grep on the raw PyMuPDF text, `90 Th ###` pattern)
correctly reaches and reads the Th block in both PDFs. The grep is not
returning an empty result because it is broken; it is empty for A=227 because
the source has nothing there.

(Note per the existing repo digest `docs/digest-literature-th-hyperfine.md`
§1.2: the 1974 value 0.46(4) μ_N for ²²⁹Th is itself superseded by two later
laser-spectroscopy/atomic-theory determinations — Safronova et al. 2013,
0.360(7) μ_N, and Porsev/Safronova/Kozlov 2021 (arXiv:2107.14723), 0.366(6)
μ_N — but that revision applies to ²²⁹Th, not ²²⁷Th, and is irrelevant to the
absence question this lookup answers.)

## Secondary literature: WebSearch summary

- No WebSearch query returned a numeric μ(²²⁷Th) value, measured or estimated,
  from any source.
- Kälber, Rink, Bekk et al., "Nuclear radii of thorium isotopes from laser
  spectroscopy of stored ions," Z. Phys. A 334, 103–108 (1989) — collinear
  laser spectroscopy on stored Th⁺ ions covering ²²⁷Th–²³⁰Th and ²³²Th.
  Search-engine abstract summaries state that isotope shifts (→ charge radii)
  were extracted for the full ²²⁷–²³²Th set, but that hyperfine splittings are
  reported "for 229Th for 3 electronic levels" specifically — not for ²²⁷Th.
  I could **not** verify this from the primary text: Springer paywalls the
  article (redirects to an IdP login), and pdf-mcp/zotero are both down this
  session. Flagging this as an open item rather than a confirmed absence — if
  ²²⁷Th HFS was measured in this paper (its odd-A hyperfine structure would in
  principle exist even though the search-engine synopsis omits it), it is not
  reflected in either Stone compilation above, which would be odd for a
  180°-different result but is not impossible if Kälber et al. did not attempt
  a magnetic-moment extraction (e.g., an unresolved/weak splitting, or no
  suitable atomic-structure calculation available in 1989 to convert A to μ).
- No other laser-spectroscopy, hyperfine-anomaly, or nuclear-structure paper
  surfaced in four targeted WebSearch queries reporting a ²²⁷Th moment.
- ENSDF's (1/2⁺) ground-state spin/parity assignment for ²²⁷Th (as already
  stated in the task brief and in the repo's existing digest) is not
  corroborated or contradicted by any moment measurement here, because none
  exists — that assignment must rest on decay-scheme/alpha-transfer
  systematics, not on a magnetic-moment or hyperfine determination.

## Conclusion

**No measured or estimated magnetic dipole moment for ²²⁷Th exists in the
literature as of the 2019 (INDC(NDS)-0794) and the earlier (~1998-cutoff)
Stone compilations, nor in the IAEA NDS live nuclear-moments database.** All
three sources place a Z=90 electromagnetic-moment entry at ²²⁹Th (μ =
+0.46(4) μ_N, 5/2⁺, ref. Gerstenkorn et al. 1974) and at ²³²Th (g-factor,
0⁺ ground state), with nothing at all for A=227 — the Z=90 row set jumps
directly over it. Four targeted WebSearch queries for a ²²⁷Th magnetic
moment or hyperfine-derived estimate returned no numeric hit; the one
plausibly relevant paper (Kälber et al. 1989, collinear laser spectroscopy of
stored Th⁺ ions across ²²⁷–²³²Th) could not be read in full this session
(Springer paywall, pdf-mcp/zotero both down) and is flagged, not resolved —
its abstract synopsis suggests hyperfine analysis was reported for ²²⁹Th
only. This corroborates and formally sources the existing repo statement
(`docs/digest-literature-th-hyperfine.md` §1.1, "no value") that ²²⁷ThF⁺
hyperfine work has no Th magnetic hyperfine constant available to predict, for
lack of any nuclear μ input.
