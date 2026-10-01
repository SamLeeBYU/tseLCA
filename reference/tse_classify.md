# Assign observations to latent classes (Step 2)

Computes posterior class-membership probabilities from a fitted
measurement model, assigns observations to classes, and estimates the
classification-error probabilities \\D\_{ts} = P(W = s \mid X = t)\\
between the true class \\X\\ and the assigned class \\W\\. These error
probabilities are what the bias-adjusted Step-3 estimators (BCH and ML)
correct for.

## Usage

``` r
tse_classify(
  object,
  newdata = NULL,
  assignment = c("modal", "proportional"),
  control = NULL
)

# S3 method for class 'tseLCA_classify'
print(x, digits = max(3L, getOption("digits") - 4L), ...)
```

## Arguments

- object:

  A measurement model from
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
  (two or more classes).

- newdata:

  Optional data frame to classify. Omitted: the data the measurement
  model was estimated on.

- assignment:

  `"modal"` (default): each observation is assigned to its most likely
  class. `"proportional"`: each observation is assigned to every class
  with its posterior probability as weight.

- control:

  Estimation settings; default: those of `object`. See
  [`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md).

- x:

  A `tseLCA_classify` object.

- digits:

  Number of significant digits to print.

- ...:

  Unused.

## Value

A `tseLCA_classify` object with components `posteriors` (n x T),
`classifications` (modal classes), `weights` (the assignment weights
\\P(W = s \mid Y_i)\\), `D` (the T x T classification-error matrix),
`entropy.R2` (computed from these posteriors), and `data`. Pass it to
the Step-3 functions.

## Details

The measurement model is held fixed. With `newdata`, observations from
another sample are classified with it, e.g. to relate the classes to
covariates observed only in a subsample; the uncertainty of the
measurement model is then still that of the sample it was estimated on.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
cl <- tse_classify(m)
cl
#> Latent class assignment (Step 2)
#>   Classes: 3   Assignment: modal   N: 500   Entropy R²: 0.8782
#> 
#> Classification error probabilities P(W = s | X = t)
#> (rows: true class X; columns: assigned class W)
#>       W=C1  W=C2  W=C3
#> X=C1 0.963 0.032 0.005
#> X=C2 0.040 0.942 0.018
#> X=C3 0.004 0.029 0.968
#> 
#> Class proportions: estimated (measurement model) and assigned
#>              C1    C2    C3
#> estimated 0.357 0.331 0.312
#> assigned  0.358 0.332 0.310
tse_classify(m, assignment = "proportional")
#> Latent class assignment (Step 2)
#>   Classes: 3   Assignment: proportional   N: 500   Entropy R²: 0.8782
#> 
#> Classification error probabilities P(W = s | X = t)
#> (rows: true class X; columns: assigned class W)
#>       W=C1  W=C2  W=C3
#> X=C1 0.935 0.058 0.007
#> X=C2 0.063 0.899 0.038
#> X=C3 0.007 0.040 0.952
#> 
#> Class proportions: estimated (measurement model) and assigned
#>              C1    C2    C3
#> estimated 0.357 0.331 0.312
#> assigned  0.357 0.331 0.312

# classify another sample with the same measurement model
tse_classify(m, newdata = d[1:200, ])
#> Latent class assignment (Step 2)
#>   Classes: 3   Assignment: modal   N: 200   Entropy R²: 0.8703
#> 
#> Classification error probabilities P(W = s | X = t)
#> (rows: true class X; columns: assigned class W)
#>       W=C1  W=C2  W=C3
#> X=C1 0.958 0.036 0.006
#> X=C2 0.048 0.945 0.007
#> X=C3 0.000 0.041 0.959
#> 
#> Class proportions: estimated (measurement model) and assigned
#>              C1    C2    C3
#> estimated 0.357 0.331 0.312
#> assigned  0.360 0.345 0.295
```
