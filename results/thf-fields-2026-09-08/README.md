# ThF+ isotope level and field plots

Work in progress, 2026-09-08. Reproducer: `scripts/plot_thf_isotopes.py`.
Separate exploratory parameter sets preserve the package defaults.

- 232Th: audited native values.
- 229Th: Zitzer 2025 nuclear moment, Skripnikov 2015 negative hyperfine factor;
  quadrupole omitted in the baseline because its signed coefficients are unknown.
  A separate legacy-quadrupole sensitivity case retains the unvalidated defaults.
- 227Th: Minkov 2024 deformed-nucleus moment, with the existing molecular factor.

All parent levels J=1–3 and all magnetic sublevels are included. Stark scans have
B=0; Zeeman scans have E=0. Calculations use J<=8 and compare J<=7. Basis
convergence is numerical convergence within this effective model, not a physical
uncertainty estimate. Shared odd-isotope B, D, dipole, and fluorine inputs remain
unscaled transfers. No quadrupole sign or normalization is invented.
