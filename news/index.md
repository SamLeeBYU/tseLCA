# Changelog

## tseLCA 1.1.0 (development version)

### Externally supplied Step-1 starting values

- Added
  [`lca_step1_startval()`](https://samleebyu.github.io/tseLCA/reference/lca_step1_startval.md),
  a wrapper around
  [`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
  that injects a user-supplied Step-1 starting value and sets
  `kmea = FALSE`, bypassing `multilevLCA`’s default deterministic
  k-means-on-principal- components initialization. This addresses
  feedback that the default initialization can consistently converge to
  a local optimum of the Step-1 log-likelihood on some datasets.
  `startval` accepts either:
  - an integer classification vector (`1..n_classes`, one entry per row
    of `data`), e.g. the modal class from an external solver run with
    many random starts (‘StepMix’, ‘poLCA’); or
  - a numeric matrix of conditional item-response probabilities
    `P(Y_h = k | X = t)`, from which a classification is derived
    internally (naive-Bayes argmax under a flat class prior). This is
    the natural format for an externally estimated Step-1 solution that
    isn’t tied to the current sample, e.g. `poLCA`’s `probs` output or a
    published item-response table. Either way, users can pass their
    external Step-1 solution straight through instead of hand-assembling
    a `tseLCA` measurement object.
- Added a matching `startval` argument to
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  and
  [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md),
  so every place a measurement model is fit from raw data can use an
  external starting value instead of `multilevLCA`’s k-means
  initialization. Supplying `startval` skips the
  `iter.measurement`/`R2.threshold` random-restart logic, since
  restarting from fresh k-means seeds would defeat the purpose of a
  user-vetted start.
- Added a `startval` argument to
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md),
  forwarded to
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  (and to
  [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md)
  when `get.twostep.vcov = TRUE`) so the measurement model can be fit
  from an external starting value in a single call.

### Multiple random-start Step-1 initialization (`n_init`)

- Added an `n_init` argument to
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md),
  [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md),
  and
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md):
  the unconditional multi-random-start analog of `n_init` in `StepMix`
  or `nrep` in `poLCA`. When supplied, the measurement model is fit
  `n_init` times from independent uniform-random classifications
  (`kmea = FALSE`, not `multilevLCA`’s deterministic k-means-on-PCA
  path), and the highest-log-likelihood fit is kept. This is a separate
  mechanism from the existing `iter.measurement`/`R2.threshold` restart
  logic, which reruns `multilevLCA`’s own k-means initialization and
  only when entropy R^2 is low; `n_init` restarts always run and never
  use k-means.
- `step1`, `startval`, and `n_init` are mutually exclusive ways of
  controlling Step 1 in
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
  (and `startval`/`n_init` in
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  and
  [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md));
  supplying more than one errors.

### Bug fixes

- Fixed the orientation of the BCH weight matrix used throughout
  `use.bch = TRUE` estimation (both covariate and distal-outcome
  models). `pwx[s, t] = P(W = s | X = t)` is column-stochastic; the
  correct BCH weight matrix is `w.is %*% t(pwx)^-1` (Mplus Web Note 21:
  the row of `H^-1` for each case’s most likely class, where
  `H = t(pwx)` is row-stochastic), not `w.is %*% pwx^-1`, which the code
  had been computing. The two orientations only agree when `pwx` is
  symmetric, so the bug was largely invisible in well-separated,
  balanced test cases; with genuine classification-error asymmetry it
  biased BCH point estimates and could produce the “negative column
  sums” error the package warns about (users were advised to fall back
  to `use.bch = FALSE`). With the corrected orientation, every row of
  the weight matrix sums to 1 and the class totals it implies exactly
  match the posterior class sizes, The weight-matrix computation is now
  consolidated into a single internal helper (`bch_weight_matrix()`)
  used by all four call sites that previously duplicated it, with a
  regression test asserting `rowSums(w.it) == 1` under a deliberately
  asymmetric classification-error matrix.

## tseLCA 1.0.0

- Initial submission to CRAN.

### Core Estimation Framework

- Implemented BCH and ML bias-adjusted three-step estimators for latent
  class analysis (LCA).
- Added support for structural models containing covariates ($`Z_p`$),
  distal outcomes ($`Z_o`$), and combined models (estimating the
  relationship between $`Z_p`$ and the latent class first, followed by
  the distal outcome adjusting for covariate-adjusted posteriors).
- Implemented analytic sandwich variance estimation to correctly
  propagate measurement uncertainty from the first-step LCA through
  classification-error correction in the final step.
- Added a robust standard error option (`use.simple.cov = TRUE`) that
  bypasses the measurement-uncertainty correction for faster computation
  in large, well-separated samples.

### Measurement Model (Step 1) Integration

- Integrated with the ‘multilevLCA’ package for efficient Step-1
  measurement model estimation.
- Added support for polytomous indicator items (0-based integer coding).
- Implemented Full Information Maximum Likelihood (FIML) to handle
  missing data in the measurement model via the `incomplete = TRUE`
  argument (using a two-pass row-filtering strategy).
- Added the ability to pass a pre-fitted measurement model (via the
  `step1` argument) to reuse across multiple structural models or apply
  to different sample subsets.
- Implemented automated random restarts for the measurement model
  triggered when entropy $`R^2`$ falls below a user-specified threshold.

### Algorithmic Flexibility & Structural Models

- Added support for both modal and proportional (soft) posterior class
  assignment (`use.modal.assignment`).
- Integrated Gaussian, Poisson, and binomial families for distal outcome
  estimation.
- Added the `rebase` argument to allow users to easily change the
  reference latent class for the multinomial logit parameterization
  while maintaining invariant log-likelihoods.
- Implemented two-step EM estimation
  ([`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md))
  to generate stable starting values for the three-step structural
  model.

### Utilities and Methods

- Included standard S3 methods for `tseLCA` objects:
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (which
  delegates to ‘multilevLCA’ for item-profile visualization).
- Built a data-generating process
  ([`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md))
  that replicates the Bakk & Kuha (2018) simulation study design for
  both covariates and distal outcomes under varying separation
  conditions.

## tseLCA 1.0.1

### CRAN resubmission

- Removed single quotes around acronyms in DESCRIPTION; added
  explanations of BCH, ML, and LCA.
- Replaced `T`/`F` with `TRUE`/`FALSE` throughout internal codebase
- Added `\value` tags to all exported functions missing them, including
  `bk2018_params`.
- `inst/examples`: examples now write to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) instead of the
  home filespace.
- `inst/examples`: commented out `rm(list = ls())` calls.
- `inst/examples`: commented out
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  calls.

## tseLCA 1.0.2

CRAN release: 2026-07-11

### CRAN resubmission

- Title: corrected hyphenated compound capitalization to “Three-Step”
  per reviewer request.
