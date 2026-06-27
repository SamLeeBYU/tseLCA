# Estimate covariate effects with measurement parameters fixed (two-step EM)

Fixes `mPhi` at `fit0$mPhi` and estimates multinomial logit coefficients
`mGamma` (Q x (T-1)) via an EM algorithm with a BFGS M-step.

## Usage

``` r
fitZ_from_fit0(
  fit0,
  data,
  Y.names,
  Zp.names,
  tol = 1e-06,
  maxIter = 200L,
  incomplete = FALSE,
  include.intercept = TRUE,
  rebase = "C1",
  starting_val = NULL,
  verbose = FALSE
)
```

## Arguments

- fit0:

  Output of `lca_step1()$fit0`.

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- Zp.names:

  Character vector of covariate column names.

- tol:

  Convergence tolerance.

- maxIter:

  Maximum EM iterations.

- incomplete:

  Logical.

- include.intercept:

  Logical.

- rebase:

  Character or integer. Reference class for the multinomial logit
  parameterization (e.g. `"C1"`, `"C2"`, or an integer). Default `"C1"`.
  Must match the `rebase` used in
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  so class column ordering is consistent.

- starting_val:

  Optional Q x (T-1) starting value matrix for `mGamma`.

- verbose:

  Logical. Print convergence messages. Default `FALSE`.

## Value

A list with `$mGamma` (Q x (T-1)), `$mPhi`, `$vOmega`, `$LLKSeries`,
`$converged`, `$n_obs`.

## Examples

``` r
# \donttest{
d  <- generate_data(200, "high", "covariate", seed = 1)
s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)

# Estimate two-step gamma with mPhi fixed at Step-1 values
fZ <- fitZ_from_fit0(
  fit0     = s1$fit0,
  data     = d,
  Y.names  = paste0("Y", 1:6),
  Zp.names = "Zp",
  verbose  = TRUE
)
#> fitZ EM converged in 8 iterations.
fZ$mGamma   # Q x (T-1) coefficient matrix
#>                   C2        C3
#> Intercept  1.8853450 -5.337580
#> Zp        -0.7516079  1.409796
fZ$converged
#> [1] TRUE
# }
```
