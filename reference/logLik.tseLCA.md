# Log-likelihood, number of observations, and information criteria

[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the
log-likelihood of a fitted model with its number of free parameters
(`df`) and observations (`nobs`), so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work directly. For a
measurement model this is the Step-1 log-likelihood. For structural
models it is the log-likelihood of the joint model for the indicators
and the structural variables evaluated with the Step-1 parameters held
fixed; for a `tseLCA_both` object it is the distal-outcome component,
which conditions on the covariates.

## Usage

``` r
# S3 method for class 'tseLCA'
logLik(object, ...)

# S3 method for class 'tseLCA'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `tseLCA` object.

- ...:

  Further arguments (currently unused).

## Value

[`logLik()`](https://rdrr.io/r/stats/logLik.html): an object of class
`"logLik"`. [`nobs()`](https://rdrr.io/r/stats/nobs.html): an integer.

## Examples

``` r
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
logLik(fit)
#> 'log Lik.' -548.6403 (df=22)
AIC(fit)
#> [1] 1141.281
BIC(fit)
#> [1] 1213.844
nobs(fit)
#> [1] 200
```
