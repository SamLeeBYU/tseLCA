## Submission

This is a patch release (2.1.1) of tseLCA, which is currently on CRAN as 2.1.0
(published 2026-10-24). It fixes one estimation error in Step 3. We submit it
soon after 2.1.0 because the error affects reported standard errors.

* The variance matrix of the maximum likelihood estimator of multinomial distal
  outcomes was computed as H^-1 M H^-1. Because H is not symmetric for this
  estimator, the right-hand factor must be the transpose of H^-1. The variance
  matrix was not symmetric and the standard errors were wrong. Point estimates
  and the other distal families are unchanged.

No functions were added or removed. See NEWS.md.

## Test environments

* local: Windows 11, R 4.5.1
* win-builder: R-devel (2026-10-30 r90805 ucrt), Windows Server 2022

## R CMD check results

0 errors | 0 warnings | 0 notes on win-builder (R-devel) and locally
(R 4.5.1, `--as-cran`).

## Downstream dependencies

There are no reverse dependencies on CRAN.
