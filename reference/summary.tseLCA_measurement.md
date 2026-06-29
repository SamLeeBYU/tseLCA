# Summarize a tseLCA model object

Verbose summary including model fit, class prevalences, item-response
probabilities, and coefficient tables with standard errors and p-values.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
summary(object, ...)

# S3 method for class 'tseLCA_covariate'
summary(object, digits = 4, ...)

# S3 method for class 'tseLCA_distal'
summary(object, digits = 4, ...)

# S3 method for class 'tseLCA_both'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- ...:

  Further arguments (currently unused).

- digits:

  Integer. Number of decimal places for coefficient tables.

## Value

Invisibly returns `object`.

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
summary(fit_m)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -299.8180
#> AIC            : 639.6360
#> BIC            : 691.7394
#> Entropy R²     : 0.8631
#> 
#> Class prevalences:
#>             
#> P(C1) 0.2614
#> P(C2) 0.3641
#> P(C3) 0.3745
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8503 0.9439 0.2512
#> P(Y2|C) 0.6930 0.9069 0.2633
#> P(Y3|C) 0.8137 0.9954 0.2800
#> P(Y4|C) 0.9197 0.1017 0.0335
#> P(Y5|C) 0.8086 0.1011 0.0583
#> P(Y6|C) 0.9830 0.0571 0.1128
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
summary(fit)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -542.3379
#> AIC            : 1164.6758
#> BIC            : 1296.6084
#> Entropy R²     : 0.8772  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.8853 -5.3376
#> Zp        -0.7516  1.4098
#> 
#> Three-step estimates:
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
summary(fit)
#> -- tseLCA Three-Step Distal Outcome Model -------------------
#> Latent classes : 3
#> Estimator      : ML
#> Family         : gaussian
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error z.value     p.value
#> mu_C1 (mean)  -0.8837    0.1185 -7.4579 < 0.001 ***
#> mu_C2 (mean)   0.9948    0.1188  8.3712 < 0.001 ***
#> mu_C3 (mean)   0.1488    0.1510  0.9852 0.3245     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
