## Submission

This is a minor release (2.1.0) of tseLCA, which is currently on CRAN as 2.0.0
(published 2026-10-01). It fixes four estimation errors in Step 3 and one flaw
in `as_tse_lca()`, which was new in 2.0.0.

* The standard errors of the BCH estimator of Poisson and Binomial distal
  outcomes used the inverse Hessian of the Gaussian case. They now use the
  variance function of the family. Point estimates are unchanged.
* The safeguard of the BCH Newton-Raphson loop for covariate effects checked
  the wrong sign of the Hessian and stopped with an error whenever the loop ran
  past a fifth of `step3.maxit` iterations. Fits that do not converge now stop
  with "BCH NR reached maximum iterations".
* In the combined covariate and distal model, the respondents were reassigned
  to classes with the covariate-adjusted posteriors, which biased the class
  means of the distal outcome. The model now keeps the Step-2 assignments.
* `tse_classify()` reports the entropy R^2 relative to the entropy of the class
  proportions, as the measurement model does.
* For a measurement model given by its parameters, the Step-1 variance used for
  the corrected standard errors was recomputed on the Step-3 data. That is only
  valid when the parameters are the maximum likelihood estimates for those
  data. `as_tse_lca()` now has a `vcov` argument for the variance matrix of the
  Step-1 parameters, which the corrected standard errors use. Without `vcov`,
  Step 3 warns and reports robust standard errors.

No functions were added or removed. See NEWS.md.

## Test environments

* local: Windows 11, R 4.5.1
* win-builder: R-devel (2026-10-23 r90760 ucrt), Windows Server 2022

## R CMD check results

0 errors | 0 warnings | 0 notes on win-builder (R-devel).

Locally (R 4.5.1, `--as-cran`) there are 2 notes, both because pandoc and 'V8'
are not installed on the check machine (README.md and NEWS.md could not be
checked; math rendering in the HTML manual was skipped).

## Downstream dependencies

There are no reverse dependencies on CRAN.
