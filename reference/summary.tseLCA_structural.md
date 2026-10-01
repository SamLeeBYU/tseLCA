# Summarize a fitted tseLCA model

[`summary()`](https://rdrr.io/r/base/summary.html) collects fit
statistics and coefficient tables; printing the result formats the
tables with
[`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html).
The coefficient table (columns `Estimate`, `Std. Error`, `z value`,
`Pr(>|z|)`) can be extracted with `coef(summary(fit))`.

## Usage

``` r
# S3 method for class 'tseLCA_structural'
summary(object, ...)

# S3 method for class 'summary.tseLCA_structural'
coef(object, ...)

# S3 method for class 'summary.tseLCA_structural'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'tseLCA_structural'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'tseLCA_measurement'
summary(object, ...)

# S3 method for class 'summary.tseLCA_measurement'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'tseLCA_measurement'
print(x, ...)
```

## Arguments

- object:

  A fitted `tseLCA` object.

- ...:

  Further arguments passed to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html).

- x:

  A `summary.tseLCA_structural` or `summary.tseLCA_measurement` object,
  or a fitted `tseLCA` object (for `print`).

- digits:

  Number of significant digits to print.

- signif.stars:

  Logical; print significance stars?

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.tseLCA_structural"` or `"summary.tseLCA_measurement"`.
Print methods return their argument invisibly.

## Examples

``` r
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
summary(fit)
#> Three-step latent class model: covariates
#>   Classes: 3   Estimator: ML   N: 200
#>   Log-lik: -548.6403 (df = 22)   AIC: 1141.28   BIC: 1213.84
#>   Entropy R² (covariate-adjusted): 0.8589
#> 
#> Covariate effects on class membership (multinomial logit):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):C2   2.2334     0.6258   3.569 0.000359 ***
#> Zp:C2           -1.1570     0.3002  -3.854 0.000116 ***
#> (Intercept):C3  -3.2742     0.7191  -4.553 5.29e-06 ***
#> Zp:C3            0.9401     0.1896   4.959 7.10e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
printCoefmat(coef(summary(fit)))
#>                Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept):C2  2.23338    0.62582  3.5688 0.0003587 ***
#> Zp:C2          -1.15699    0.30017 -3.8545 0.0001160 ***
#> (Intercept):C3 -3.27422    0.71914 -4.5529 5.290e-06 ***
#> Zp:C3           0.94007    0.18958  4.9587 7.098e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
