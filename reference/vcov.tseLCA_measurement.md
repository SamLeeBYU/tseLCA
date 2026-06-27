# Extract the variance-covariance matrix from a tseLCA model object

Extract the variance-covariance matrix from a tseLCA model object

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
vcov(object, ...)

# S3 method for class 'tseLCA_covariate'
vcov(object, which = c("three_step", "two_step"), ...)

# S3 method for class 'tseLCA_distal'
vcov(object, ...)

# S3 method for class 'tseLCA_both'
vcov(
  object,
  which = c("covariate", "distal", "both"),
  step = c("three_step", "two_step"),
  ...
)
```

## Arguments

- object:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- ...:

  Further arguments (currently unused).

- which:

  Character. `"three_step"` (default) or `"two_step"` for covariate
  models; `"covariate"`, `"distal"`, or `"both"` for both models.

- step:

  Character. For `tseLCA_both`: `"three_step"` (default) or
  `"two_step"`.

## Value

A named square matrix. The two-step vcov is only available when
`get.twostep.vcov = TRUE` was set in
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
vcov(fit_m)   # returns NULL with a message
#> No variance-covariance matrix available for measurement-only models.
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
vcov(fit)   # Q*(T-1) x Q*(T-1) vcov matrix with named rows/cols
#>              Intercept:C2         Zp:C2  Intercept:C3        Zp:C3
#> Intercept:C2  0.240113667 -0.0840137884  0.0197404469 -0.009578627
#> Zp:C2        -0.084013788  0.0357512662 -0.0001681731  0.003017456
#> Intercept:C3  0.019740447 -0.0001681731  1.6077453606 -0.400263916
#> Zp:C3        -0.009578627  0.0030174563 -0.4002639161  0.103368072
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
vcov(fit)   # T x T vcov matrix with mu_C1..mu_CT row/col names
#>               mu_C1         mu_C2         mu_C3
#> mu_C1  0.0140394227  0.0002115275 -0.0001400022
#> mu_C2  0.0002115275  0.0141216292 -0.0006978379
#> mu_C3 -0.0001400022 -0.0006978379  0.0228157389
# }
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 0.5)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", Zo.name = "Zo",
                  use.simple.cov = TRUE)
vcov(fit, which = "covariate")
#>              Intercept:C2         Zp:C2  Intercept:C3        Zp:C3
#> Intercept:C2  0.240113667 -0.0840137884  0.0197404469 -0.009578627
#> Zp:C2        -0.084013788  0.0357512662 -0.0001681731  0.003017456
#> Intercept:C3  0.019740447 -0.0001681731  1.6077453606 -0.400263916
#> Zp:C3        -0.009578627  0.0030174563 -0.4002639161  0.103368072
vcov(fit, which = "distal")
#>              mu_C1        mu_C2        mu_C3
#> mu_C1 3.863995e-03 4.056130e-05 1.833309e-05
#> mu_C2 4.056130e-05 4.412576e-03 2.247326e-05
#> mu_C3 1.833309e-05 2.247326e-05 5.646743e-03
# }
```
