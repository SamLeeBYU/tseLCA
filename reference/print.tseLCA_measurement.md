# Print a tseLCA model object

Compact one-line or table summary printed to the console.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
print(x, ...)

# S3 method for class 'tseLCA_covariate'
print(x, digits = 4, ...)

# S3 method for class 'tseLCA_distal'
print(x, digits = 4, ...)

# S3 method for class 'tseLCA_both'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- ...:

  Further arguments.

- digits:

  Integer. Number of decimal places for coefficient tables.

## Value

Invisibly returns `x`.

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
print(fit_m)
#> tseLCA -- measurement model
#>   Classes: 3   Log-lik: -299.8180   AIC: 639.64   BIC: 691.74
#>   Entropy R²: 0.8631
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
print(fit)
#> tseLCA -- three-step covariate model
#>   Classes: 3   Estimator: ML   Log-lik: -542.3379   AIC: 1164.68   BIC: 1296.61
#> 
#> Covariate coefficients (three-step):
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.0295    0.4900  4.1416 < 0.001 ***
#> Zp:C2         -0.8192    0.1891 -4.3323 < 0.001 ***
#> Intercept:C3  -5.4875    1.2680 -4.3278 < 0.001 ***
#> Zp:C3          1.4540    0.3215  4.5223 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
print(fit)
#> tseLCA -- three-step distal outcome model
#>   Classes: 3   Estimator: ML   Family: gaussian
#> 
#> Distal outcome means by class:
#>              Estimate Std.Error z.value     p.value
#> mu_C1 (mean)  -0.8837    0.1185 -7.4579 < 0.001 ***
#> mu_C2 (mean)   0.9948    0.1188  8.3712 < 0.001 ***
#> mu_C3 (mean)   0.1488    0.1510  0.9852 0.3245     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
