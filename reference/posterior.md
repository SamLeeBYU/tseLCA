# Posterior class-membership probabilities and modal class assignments

`posterior()` returns the n x T matrix of posterior class-membership
probabilities used by a fitted model; `classes()` returns the modal
(most likely) class of each observation.

## Usage

``` r
posterior(object, ...)

# S3 method for class 'tseLCA'
posterior(object, ...)

classes(object, ...)

# S3 method for class 'tseLCA'
classes(object, ...)
```

## Arguments

- object:

  A fitted `tseLCA` object.

- ...:

  Further arguments (currently unused).

## Value

`posterior()`: a numeric n x T matrix. `classes()`: an integer vector of
length n with values in `1..T`.

## Examples

``` r
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
head(posterior(fit))
#>                C1          C2           C3
#> [1,] 9.985230e-01 0.001476673 3.472258e-07
#> [2,] 9.985230e-01 0.001476673 3.472258e-07
#> [3,] 5.264000e-03 0.992771673 1.964327e-03
#> [4,] 5.264000e-03 0.992771673 1.964327e-03
#> [5,] 7.220078e-06 0.002770271 9.972225e-01
#> [6,] 9.985082e-01 0.001458223 3.353288e-05
table(classes(fit))
#> 
#>  1  2  3 
#> 68 59 73 
```
