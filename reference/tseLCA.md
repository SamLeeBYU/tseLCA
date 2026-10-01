# Three-step latent class analysis in one call

Fits the measurement model, classifies the observations, and relates the
classes to covariates and/or a distal outcome, all from one formula.
This is a convenience wrapper around the step-wise functions
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md),
[`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md),
[`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md),
and
[`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md),
which remain available for inspecting each step (see
[`measurement()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)).

## Usage

``` r
tseLCA(
  formula,
  data,
  nclass,
  family = "gaussian",
  method = c("ML", "BCH", "none"),
  se = c("corrected", "robust"),
  assignment = c("modal", "proportional"),
  ref = 1,
  missing = c("listwise", "fiml"),
  start = NULL,
  control = tse_control()
)
```

## Arguments

- formula:

  `cbind(indicators) ~ covariates | distal outcome`; see Details.

- data:

  A data frame.

- nclass:

  Number of latent classes.

- family:

  Distribution of the distal outcome; see
  [`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md).

- method:

  Step-3 estimator: `"ML"`, `"BCH"`, or `"none"`; see
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).

- se:

  Standard errors: `"corrected"` or `"robust"`.

- assignment:

  Step-2 class assignment: `"modal"` or `"proportional"`.

- ref:

  Reference class of the covariate model.

- missing:

  Missing indicator values: `"listwise"` or `"fiml"`; see
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

- start:

  Optional Step-1 starting values; see
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

- control:

  Estimation settings, see
  [`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md).

## Value

A `tseLCA_measurement`, `tseLCA_covariate`, `tseLCA_distal`, or
`tseLCA_both` object, depending on the formula. Its components are
available with
[`measurement()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`classification()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`covariate()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
and
[`distal()`](https://samleebyu.github.io/tseLCA/reference/measurement.md).

## Details

The formula has up to three parts: `indicators ~ covariates | outcome`.

- Left-hand side: the indicators, `cbind(Y1, Y2, ...)`.

- First right-hand side part: the covariates of class membership, with
  the usual formula syntax (`1` for none).

- Optional second part: the name of a distal outcome.

For example, `cbind(Y1, Y2, Y3) ~ age + sex | income` relates the
classes to the covariates `age` and `sex` and to the distal outcome
`income`; `cbind(Y1, Y2, Y3) ~ 1 | income` has only the distal outcome;
and `cbind(Y1, Y2, Y3) ~ 1` fits the measurement model alone.

The number of classes is chosen beforehand from the measurement model,
for example with `tse_lca(..., nclass = 1:6)`.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)

fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
summary(fit)
#> Three-step latent class model: covariates and distal outcome
#>   Classes: 3   Estimator: ML   Family: gaussian   N: 500
#>   Log-lik: -2053.7961 (df = 26)   AIC: 4159.59   BIC: 4269.17
#>   Entropy R² (covariate-adjusted): 0.8693
#> 
#> Covariate effects on class membership (multinomial logit):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):C2   2.0411     0.3264   6.254 4.01e-10 ***
#> Zp:C2           -0.8821     0.1430  -6.168 6.91e-10 ***
#> (Intercept):C3  -3.4836     0.6211  -5.609 2.04e-08 ***
#> Zp:C3            0.8985     0.1485   6.049 1.46e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Distal outcome means by class:
#>       Estimate Std. Error z value Pr(>|z|)    
#> mu_C1 -1.03985    0.07961 -13.062   <2e-16 ***
#> mu_C2  0.94707    0.08158  11.609   <2e-16 ***
#> mu_C3  0.09331    0.08508   1.097    0.273    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
measurement(fit)
#> Latent class measurement model
#>   Classes: 3   N: 500
#>   Log-lik: -1455.5052 (df = 20)   AIC: 2951.01   BIC: 3035.30
#>   Entropy R²: 0.8780
classification(fit)
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

# the same model, step by step
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
fc <- tse_covariate(tse_classify(m), ~ Zp)
fb <- tse_distal(fc, Zo ~ 1)
```
