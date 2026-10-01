# Two-step estimates of covariate effects

Estimates the multinomial logit of class membership on covariates with
the measurement-model parameters held fixed at their Step-1 values (Bakk
and Kuha 2018). Unlike the three-step estimators, the indicators enter
the Step-2 likelihood directly, so no classification step is needed.

## Usage

``` r
tse_twostep(object, formula, ref = 1, se = FALSE, control = NULL)
```

## Arguments

- object:

  A measurement model from
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
  (it must keep its data).

- formula:

  One-sided covariate formula.

- ref:

  Reference class of the multinomial logit.

- se:

  Logical. If `TRUE`, the estimates and their standard errors (corrected
  for the Step-1 uncertainty) are obtained with the two-step estimator
  of multilevLCA, initialized at this model's classes; its measurement
  model is checked against `object`. If `FALSE` (default), only the
  estimates are computed, and the variance is `NA`.

- control:

  Estimation settings; default: those of `object`.

## Value

A `tseLCA_twostep` object (also a `tseLCA_covariate`).

## References

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
coef(tse_twostep(m, ~ Zp))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.1934130     -0.9411383     -3.4524271      0.8971774 
```
