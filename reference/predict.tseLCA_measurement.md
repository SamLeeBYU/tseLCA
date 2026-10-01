# Class membership predictions from a measurement model

Posterior class-membership probabilities P(X = t \| Y) for the rows of
`newdata` (or the estimation sample), or their modal class. Indicators
in `newdata` are coded with the categories stored in the model; missing
indicator values are skipped (the posterior uses the observed ones), and
rows with no observed indicator get `NA`.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
predict(object, newdata = NULL, type = c("posterior", "class"), ...)

# S3 method for class 'tseLCA_measurement'
fitted(object, ...)
```

## Arguments

- object:

  A `tseLCA_measurement` object.

- newdata:

  Optional data frame with the indicator columns. Omitted: the
  estimation sample.

- type:

  `"posterior"` (default) for an n x T matrix of probabilities, or
  `"class"` for the modal class of each row.

- ...:

  Unused.

## Value

A matrix (`type = "posterior"`) or integer vector (`"class"`).

## Examples

``` r
d <- generate_data(300, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
predict(m, newdata = d[1:5, ])
#>             C1           C2           C3
#> 1 0.0633839383 0.9350917888 1.524273e-03
#> 2 0.0009948482 0.9461235966 5.288156e-02
#> 3 0.9990289946 0.0008413575 1.296479e-04
#> 4 0.9520799155 0.0478234960 9.658856e-05
#> 5 0.0117728922 0.1877192401 8.005079e-01
predict(m, newdata = d[1:5, ], type = "class")
#> [1] 2 2 1 1 3
```
