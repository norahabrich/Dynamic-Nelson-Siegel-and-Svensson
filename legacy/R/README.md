# Original R code (Werley Cordeiro)

These scripts are the original two-step Dynamic Nelson-Siegel and Svensson code
written by [Werley Cordeiro](https://github.com/werleycordeiro) in 2019-2020.
They are kept here for reference and credit. The Python package in `dnsmodel/`
is a rewrite and extension of this work.

Changes made here, so that the scripts run from this folder:

- data are read from `../../data/bonds.csv` instead of a Dropbox link;
- `DNSS-TS.R` now passes the loading matrix to `Yhat.betas()` and the sample
  size to `VARcoeff()` (both arguments were missing and the script stopped);
- the unused `highfrequency` package is no longer required.

Run with `Rscript DNS-TS.R` from this folder (needs `xts`, `vars`, `YieldCurve`).
