# Estimate two-step covariate model with multilevLCA (optional reference path)

Calls
[`multilevLCA::multiLCA`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
with `fixedpars = 1` and `Z = Zp.names` to fit the two-step covariate
model. This is the original multilevLCA approach and is used when
`get.twostep.vcov = TRUE` in
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
to obtain multilevLCA's corrected standard errors for the two-step gamma
estimates.

## Usage

``` r
fitZ_from_multiLCA(
  data,
  Y.names,
  n_classes,
  Zp.names,
  maxIter.measurement,
  measurement.tol,
  covariate.tol,
  iter.measurement,
  R2.threshold,
  incomplete = FALSE,
  rebase = "C1",
  startval = NULL,
  n_init = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- n_classes:

  Integer. Number of latent classes.

- Zp.names:

  Character vector of covariate column names.

- maxIter.measurement:

  Maximum EM iterations.

- measurement.tol:

  Convergence tolerance.

- covariate.tol:

  NR tolerance for the covariate model.

- iter.measurement:

  Number of random restarts.

- R2.threshold:

  Entropy R\\^2\\ restart threshold.

- incomplete:

  Logical. FIML for partially missing indicators. See the `Missing Data`
  section of `vignette("tseLCA", package = "tseLCA")`. Default `FALSE`.

- rebase:

  Character or integer. Reference class for column naming of `$mGamma`.
  Must match the `rebase` used in
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
  so coefficient labels are consistent. Default `"C1"`.

- startval:

  Optional starting classification for the measurement portion of this
  `multiLCA(fixedpars = 1)` fit – an integer vector or a conditional
  item-response probability matrix, as described in
  [`lca_step1_startval()`](https://samleebyu.github.io/tseLCA/reference/lca_step1_startval.md)
  – e.g. the same value passed to
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  for the primary Step-1 fit. When supplied, `kmea = FALSE` is used and
  `iter.measurement`/`R2.threshold` restarts are skipped, for the same
  reasons as in
  [`lca_step1_startval()`](https://samleebyu.github.io/tseLCA/reference/lca_step1_startval.md).
  Mutually exclusive with `n_init`. Default `NULL`.

- n_init:

  Optional positive integer. If supplied, fits this
  `multiLCA(fixedpars = 1)` model `n_init` times from independent
  uniform-random classifications (`kmea = FALSE`) and keeps the fit with
  the highest log-likelihood, as in
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)'s
  `n_init` argument. `iter.measurement`/`R2.threshold` restarts are
  skipped. Mutually exclusive with `startval`. Default `NULL`.

- verbose:

  Logical.

## Value

A list with the following elements:

- `mGamma`:

  (Q+1) x (T-1) numeric matrix of multinomial logit coefficients. Rows
  are named by covariate (including `"Intercept"`), columns by
  non-reference class (e.g. `"C2"`, `"C3"`).

- `mPhi`:

  Item parameter matrix (items x classes) from the fixed-parameter
  multilevLCA fit.

- `vOmega`:

  Length-T vector of marginal class proportions, computed as the average
  of the fitted class probability matrix (`vPi_avg` in multilevLCA
  output).

- `LLKSeries`:

  Matrix of observed-data log-likelihoods across EM iterations, passed
  through directly from the multilevLCA fit.

- `raw_fit`:

  The full
  [`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
  output object, including `$Varmat_cor` (corrected variance matrix) and
  `$SEs_cor_gamma` (corrected standard errors for `mGamma`) if
  available.

## Deprecated

Deprecated as of tseLCA 2.0.0: use
[`tse_twostep()`](https://samleebyu.github.io/tseLCA/reference/tse_twostep.md).
It keeps working and warns once per session when called directly.

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Two-step estimation with multiLCA (fixedpars = 1)
fZ_ml <- fitZ_from_multiLCA(
  data                = d,
  Y.names             = paste0("Y", 1:6),
  n_classes           = 3,
  Zp.names            = "Zp",
  maxIter.measurement = 5000L,
  measurement.tol     = 1e-8,
  covariate.tol       = 1e-6,
  iter.measurement    = 10L,
  R2.threshold        = 0.70
)
#> Warning: fitZ_from_multiLCA() is deprecated as of tseLCA 2.0.0; use tse_twostep(se = TRUE). (Shown once per session; see NEWS for the new interface.)
fZ_ml$mGamma           # two-step estimates
#>                    C2         C3
#> (Intercept)  1.990672 -3.1319909
#> Zp          -1.018352  0.9190157
fZ_ml$raw_fit$Varmat_cor   # multilevLCA corrected vcov
#>             [,1]         [,2]        [,3]         [,4]
#> [1,]  0.30784133 -0.126954327  0.06746616 -0.020557798
#> [2,] -0.12695433  0.072299967 -0.01991673  0.008336796
#> [3,]  0.06746616 -0.019916726  0.69429491 -0.177175541
#> [4,] -0.02055780  0.008336796 -0.17717554  0.048705842
# }
```
