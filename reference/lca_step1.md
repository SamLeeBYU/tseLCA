# Fit the LCA measurement model (Step 1)

Estimates the latent class measurement model via multilevLCA and,
optionally, fixes `mPhi` and estimates covariate effects (two-step
initialization) via
[`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).

## Usage

``` r
lca_step1(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-08,
  covariate.tol = 1e-06,
  iter.measurement = 10L,
  R2.threshold = 0.7,
  use.two.step = TRUE,
  estimate.one.step = TRUE,
  incomplete = FALSE,
  maxIter.fitZ = 200L,
  include.intercept = TRUE,
  rebase = "C1",
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame containing at minimum the indicator columns.

- Y.names:

  Character vector of item column names.

- n_classes:

  Integer. Number of latent classes.

- Zp.names:

  Character vector of covariate column names, or `NULL`.

- maxIter.measurement:

  Maximum EM iterations. Default `5000L`.

- measurement.tol:

  Convergence tolerance. Default `1e-8`.

- covariate.tol:

  Convergence tolerance for the `fitZ` BFGS M-step.

- iter.measurement:

  Number of random restarts when entropy R\\^2\\ is low.

- R2.threshold:

  Entropy R\\^2\\ below which restarts are triggered.

- use.two.step:

  Logical. If `TRUE`, also estimate `fitZ` via
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).

- estimate.one.step:

  Logical. If `FALSE`, skip the unconditional EM and only compute
  `fitZ`.

- incomplete:

  Logical. FIML for partially missing indicators.

- maxIter.fitZ:

  Maximum BFGS-EM iterations for
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).

- include.intercept:

  Logical. Prepend intercept to covariate design matrix.

- rebase:

  Character or integer specifying the reference latent class. Use
  `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
  measurement model is permuted so this class becomes column 1, making
  it the reference for all downstream multinomial logit
  parameterizations.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A list with `$fit0` (multilevLCA measurement model) and `$fitZ`
(two-step covariate model from `fitZ_from_fit0`, or `NULL`).

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Measurement model only
s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
s1$fit0$vPi    # estimated class prevalences
#>                
#> P(C1) 0.3202142
#> P(C2) 0.3663893
#> P(C3) 0.3133965
s1$fit0$mPhi   # item-response probabilities
#>                C1         C2        C3
#> P(Y1|C) 0.8671917 0.84615644 0.0870791
#> P(Y2|C) 0.8662547 0.91259281 0.1539543
#> P(Y3|C) 0.9669370 0.95701501 0.1746446
#> P(Y4|C) 0.9591301 0.09549431 0.1368390
#> P(Y5|C) 0.8865452 0.12842286 0.1725065
#> P(Y6|C) 0.8980182 0.13578510 0.1043140

# With two-step covariate initialization
s1z <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3,
                 Zp.names = "Zp", use.two.step = TRUE, verbose = TRUE)
#> fitZ EM converged in 8 iterations.
s1z$fitZ$mGamma   # two-step gamma estimates
#>                   C2        C3
#> Intercept  1.8853450 -5.337580
#> Zp        -0.7516079  1.409796
# }
```
