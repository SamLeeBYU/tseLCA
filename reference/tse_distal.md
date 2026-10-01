# Relate latent classes to a distal outcome (Step 3)

Estimates the class-specific distribution of a distal outcome,
correcting for the classification error of the Step-2 class assignments
(Bakk, Tekle, and Vermunt 2013). Given a covariate model from
[`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md),
the class prior depends on the covariates and the covariate-model
uncertainty is propagated to the distal estimates.

## Usage

``` r
tse_distal(
  object,
  formula,
  family = "gaussian",
  method = NULL,
  se = NULL,
  control = NULL,
  data = NULL
)
```

## Arguments

- object:

  A classification from
  [`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md),
  or a covariate model from
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md)
  (combined model).

- formula:

  `Zo ~ 1`, with the distal outcome on the left-hand side.

- family:

  Distribution of the outcome within classes: `"gaussian"` (default),
  `"poisson"`, `"binomial"`, `"multinomial"` (nominal outcome), or the
  corresponding family object
  ([`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html); canonical links
  only).

- method, se:

  As for
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).
  For a combined model they default to those of the covariate model.

- control:

  Estimation settings; default: those of `object`.

- data:

  Optional data frame with the distal outcome, as in
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).

## Value

A `tseLCA_distal` object, or a `tseLCA_both` object when `object` is a
covariate model.

## Details

The class-specific parameters are means (`gaussian`, with a common
within-class variance, reported as `$sigma2`), log means (`poisson`),
logits (`binomial`), or category probabilities (`"multinomial"`). The
estimators and standard errors are as for
[`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).
Use
[`omnibus_test()`](https://samleebyu.github.io/tseLCA/reference/omnibus_test.md)
to test whether the outcome differs across classes.

## References

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

## Examples

``` r
d <- generate_data(500, "high", "distal", seed = 2)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
fd <- tse_distal(tse_classify(m, assignment = "proportional"), Zo ~ 1)
summary(fd)
#> Three-step latent class model: distal outcome
#>   Classes: 3   Estimator: ML   Family: gaussian   N: 500
#>   Log-lik: -2220.3460 (df = 24)   AIC: 4488.69   BIC: 4589.84
#> 
#> Distal outcome means by class:
#>       Estimate Std. Error z value Pr(>|z|)    
#> mu_C1 -0.96436    0.09113 -10.582   <2e-16 ***
#> mu_C2  1.02629    0.07548  13.596   <2e-16 ***
#> mu_C3  0.10413    0.09979   1.044    0.297    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
omnibus_test(fd)
#> 
#>  Wald test of equal distal outcome distributions across latent classes
#> 
#> data:  gaussian distal outcome, 3 classes
#> W = 278.59, df = 2, p-value < 2.2e-16
#> 
```
