# Class enumeration results

Methods for the `tseLCA_select` object returned by
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
with several numbers of classes. `best_model()` returns the fitted model
that minimizes an information criterion; `x[[k]]` returns the `k`-class
model.

## Usage

``` r
best_model(object, ...)

# S3 method for class 'tseLCA_select'
best_model(object, criterion = c("BIC", "AIC", "SABIC"), ...)

# S3 method for class 'tseLCA_select'
x[[i, ...]]

# S3 method for class 'tseLCA_select'
as.data.frame(x, ...)

# S3 method for class 'tseLCA_select'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'tseLCA_select'
plot(x, which = c("AIC", "BIC", "SABIC"), ...)
```

## Arguments

- ...:

  Further arguments passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html)
  (`plot`) or unused.

- criterion:

  Information criterion to minimize: `"BIC"` (default), `"AIC"`, or
  `"SABIC"`.

- x, object:

  A `tseLCA_select` object.

- i:

  Number of classes of the model to extract.

- digits:

  Number of significant digits to print.

- which:

  Criteria to plot.

## Value

`best_model()` and `[[`: a `tseLCA_measurement` object.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): the
enumeration table. [`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html): `x`,
invisibly.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
sel <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 1:4)
as.data.frame(sel)
#>   nclass    logLik npar      AIC      BIC    SABIC entropy.R2 min.class nobs
#> 1      1 -1965.822    6 3943.645 3968.932 3949.888         NA 1.0000000  500
#> 2      2 -1615.803   13 3257.606 3312.396 3271.133  0.8251766 0.4326791  500
#> 3      3 -1455.505   20 2951.010 3035.303 2971.821  0.8780302 0.3122150  500
#> 4      4 -1452.154   27 2958.308 3072.103 2986.403  0.7889633 0.1569247  500
best_model(sel)
#> Latent class measurement model
#>   Classes: 3   N: 500
#>   Log-lik: -1455.5052 (df = 20)   AIC: 2951.01   BIC: 3035.30
#>   Entropy R²: 0.8780
sel[[2]]
#> Latent class measurement model
#>   Classes: 2   N: 500
#>   Log-lik: -1615.8031 (df = 13)   AIC: 3257.61   BIC: 3312.40
#>   Entropy R²: 0.8252
```
