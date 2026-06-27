# Estimate two-step covariate model via multilevLCA (optional reference path)

Calls
[`multilevLCA::multiLCA`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
with `fixedpars = 1` and `Z = Zp.names` to fit the two-step covariate
model. This is the original multilevLCA approach and is used when
`get.twostep.vcov = TRUE` in `teLCA::three_step()` to obtain
multilevLCA's corrected standard errors for the two-step gamma
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

  Logical.

- rebase:

  Character or integer. Reference class for column naming of `$mGamma`.
  Must match the `rebase` used in `teLCA::three_step()` so coefficient
  labels are consistent. Default `"C1"`.

- verbose:

  Logical.

## Value

A list with `$mGamma`, `$mPhi`, `$vOmega`, `$LLKSeries`, and `$raw_fit`
(the full multilevLCA output, including `$Varmat_cor` and
`$SEs_cor_gamma` if available).

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Two-step estimation via multiLCA (fixedpars = 1)
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
fZ_ml$mGamma           # two-step estimates
#>                   C2        C3
#> Intercept  1.8862448 -5.343414
#> Zp        -0.7520118  1.411197
fZ_ml$raw_fit$Varmat_cor   # multilevLCA corrected vcov
#>              [,1]         [,2]         [,3]         [,4]
#> [1,]  0.226606015 -0.080420786  0.022531107 -0.008956779
#> [2,] -0.080420786  0.036221462  0.004946794  0.001391765
#> [3,]  0.022531107  0.004946794  1.501589310 -0.345447542
#> [4,] -0.008956779  0.001391765 -0.345447542  0.082865708
# }
```
