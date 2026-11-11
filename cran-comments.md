## Submission

This is a minor update (2.1.2) of tseLCA, which is currently on CRAN as 2.1.1.
It changes the printed output and the help index only (per one reviewer's request from JSS). No estimates, standard
errors, or function signatures change.

* `print()` of a fitted structural model is now brief (header, fit statistics,
  and estimates with standard errors). It previously duplicated `summary()`,
  which remains the full report.
* The deprecated `three_step()` is now marked `\keyword{internal}` so it no
  longer appears in the package help index. It is still exported and works as
  before.
* `tse_twostep()` has a clearer help title.

No functions were added or removed. See NEWS.md.

## Test environments

* local: Windows 11, R 4.5.1
* win-builder: R-devel (2026-11-10 r90874 ucrt), Windows Server 2022

## R CMD check results

0 errors | 0 warnings | 0 notes on win-builder (R-devel) and locally
(R 4.5.1, `--as-cran`).

## Downstream dependencies

There are no reverse dependencies on CRAN.
