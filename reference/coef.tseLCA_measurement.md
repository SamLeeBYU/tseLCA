# Extract coefficients from a tseLCA model object

Extract coefficients from a tseLCA model object

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
coef(object, ...)

# S3 method for class 'tseLCA_covariate'
coef(object, which = c("three_step", "two_step"), ...)

# S3 method for class 'tseLCA_distal'
coef(object, ...)

# S3 method for class 'tseLCA_both'
coef(
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

  Character. For covariate and both models: `"three_step"` (default) or
  `"two_step"`. For both models also accepts `"covariate"`, `"distal"`,
  or `"both"`.

- step:

  Character. For `tseLCA_both`: `"three_step"` (default) or
  `"two_step"`.

## Value

The coefficient matrix (covariate models), named numeric vector (distal
models), or a named list of both (measurement or both models).

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
coef(fit_m)   # returns list with $prevalences and $item_probs
#> $prevalences
#>                
#> P(C1) 0.2613812
#> P(C2) 0.3641288
#> P(C3) 0.3744900
#> 
#> $item_probs
#>                C1         C2         C3
#> P(Y1|C) 0.8503009 0.94386031 0.25116991
#> P(Y2|C) 0.6930271 0.90692158 0.26334344
#> P(Y3|C) 0.8136952 0.99542535 0.27998700
#> P(Y4|C) 0.9197189 0.10173196 0.03353693
#> P(Y5|C) 0.8085871 0.10111305 0.05829891
#> P(Y6|C) 0.9829910 0.05714552 0.11283694
#> 
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
coef(fit)                      # three-step estimates
#>                   C2        C3
#> Intercept  2.0294548 -5.487549
#> Zp        -0.8191555  1.453963
coef(fit, which = "two_step")  # two-step starting values
#>                   C2        C3
#> Intercept  1.8853450 -5.337580
#> Zp        -0.7516079  1.409796
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
coef(fit)   # named vector of class means
#>      mu_C1      mu_C2      mu_C3 
#> -0.8836699  0.9947896  0.1488161 
# }
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", Zo.name = "Zo",
                  use.simple.cov = TRUE)
coef(fit, which = "covariate")
#>                   C2        C3
#> Intercept  2.0294548 -5.487549
#> Zp        -0.8191555  1.453963
coef(fit, which = "distal")
#>        mu_C1        mu_C2        mu_C3 
#> -1.068160593  0.006733996  1.191464437 
coef(fit, which = "both")
#> $covariate
#>                   C2        C3
#> Intercept  2.0294548 -5.487549
#> Zp        -0.8191555  1.453963
#> 
#> $distal
#>        mu_C1        mu_C2        mu_C3 
#> -1.068160593  0.006733996  1.191464437 
#> 
# }
```
