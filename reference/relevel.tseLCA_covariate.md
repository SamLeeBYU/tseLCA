# Change the reference class of a covariate model

Refits the model with another reference class of the multinomial logit.

## Usage

``` r
# S3 method for class 'tseLCA_covariate'
relevel(x, ref, ...)
```

## Arguments

- x:

  A `tseLCA_covariate` object from
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).

- ref:

  The new reference class (number or label such as `"C2"`).

- ...:

  Unused.

## Value

The refitted `tseLCA_covariate` object.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
fit <- tse_covariate(tse_classify(m), ~ Zp)
coef(stats::relevel(fit, ref = "C3"))
#> (Intercept):C1          Zp:C1 (Intercept):C2          Zp:C2 
#>      3.4835609     -0.8984976      5.5246368     -1.7805775 
```
