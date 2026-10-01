# Variance-covariance matrix of a fitted tseLCA model

Row and column names match
[`coef()`](https://rdrr.io/r/stats/coef.html).

## Usage

``` r
# S3 method for class 'tseLCA_structural'
vcov(
  object,
  component = c("all", "covariate", "distal"),
  step = c("three_step", "two_step"),
  ...
)

# S3 method for class 'tseLCA_measurement'
vcov(object, boundary.tol = 0.01, ...)
```

## Arguments

- object:

  A fitted `tseLCA` object.

- component:

  For `tseLCA_both` objects: `"all"` (default; covariate then distal
  coefficients), `"covariate"`, or `"distal"`.

- step:

  `"three_step"` (default) or `"two_step"` (the two-step estimates used
  to initialize Step 3; covariate models only).

- ...:

  Further arguments (currently unused).

- boundary.tol:

  Measurement models only: parameters within this tolerance of 0 or 1
  are treated as fixed. Default `1e-2`.

## Value

A named square matrix.

## Details

- Measurement models: the BHHH variance matrix of the log-ratio
  parameters (attribute `"parameterization"` records the scale).

- Covariate and distal-outcome models: the Step-3 variance matrix, which
  includes the correction for Step-1 uncertainty unless the model was
  fitted with `use.simple.cov = TRUE`. For `family = "multinomial"` it
  is on the probability scale and rank-deficient (each class's
  probabilities sum to one).

- `tseLCA_both` with `component = "all"`: the covariate and distal
  blocks on the diagonal; the cross-covariances between the two sets of
  parameters are not computed and are `NA`.

## Examples

``` r
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
vcov(fit)
#>                (Intercept):C2        Zp:C2 (Intercept):C3        Zp:C3
#> (Intercept):C2    0.391644881 -0.173583653    0.001643327 -0.002599746
#> Zp:C2            -0.173583653  0.090099886    0.016315352 -0.002300251
#> (Intercept):C3    0.001643327  0.016315352    0.517169347 -0.130664301
#> Zp:C3            -0.002599746 -0.002300251   -0.130664301  0.035941355
```
