# Class-membership probabilities from a covariate model

The fitted class prior \\P(X = t \mid Z)\\ for the rows of `newdata`, or
of the estimation data.

## Usage

``` r
# S3 method for class 'tseLCA_covariate'
predict(object, newdata = NULL, type = c("prob", "class"), ...)
```

## Arguments

- object:

  A `tseLCA_covariate` object from
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).

- newdata:

  Optional data frame with the covariates.

- type:

  `"prob"` (default) for the n x T probability matrix, or `"class"` for
  the most likely class.

- ...:

  Unused.

## Value

A matrix (rows with missing covariates are `NA`) or integer vector.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
fit <- tse_covariate(tse_classify(m), ~ Zp)
predict(fit, newdata = data.frame(Zp = 1:5))
#>          C1         C2         C3
#> 1 0.2346248 0.74768656 0.01768866
#> 2 0.3993275 0.52673529 0.07393719
#> 3 0.4998232 0.27289595 0.22728081
#> 4 0.4268484 0.09646544 0.47668621
#> 5 0.2606745 0.02438452 0.71494097
```
