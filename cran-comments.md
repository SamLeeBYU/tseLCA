## Submission

This is a major release (2.0.0) of tseLCA, which is currently on CRAN as
1.1.1. The package now provides a step-wise interface with a formula interface
and S3 classes and methods for each step (`tse_lca()`, `tse_classify()`,
`tse_covariate()`, `tse_distal()`, `tse_twostep()`, and the one-call
`tseLCA()`). The 1.x function `three_step()` is kept as a deprecated wrapper.
Several estimation bugs are fixed; see NEWS.md.

This update follows 1.1.1 (published 2026-09-23) closely because 1.1.1
returns incorrect results in several cases:

* bias-corrected standard errors were wrong (Step-2 Jacobian column order and
  modal-assignment derivative) or `NA` when an item probability was on the
  boundary;
* classes and distal parameters were mislabelled when the reference class
  was not the first class;
* maximum likelihood distal-outcome estimates were biased under proportional
  assignment, and the Gaussian variance was fixed at 1;
* AIC/BIC used an incorrect number of parameters.

The package is also the subject of a manuscript being revised for the
Journal of Statistical Software, and the revision depends on this version.

## Test environments

* local: Windows 11, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no reverse dependencies on CRAN.
