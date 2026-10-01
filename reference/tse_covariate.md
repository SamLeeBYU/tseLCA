# Relate latent classes to covariates (Step 3)

Estimates a multinomial logistic regression of latent class membership
on covariates, \\P(X = t \mid Z) \propto \exp(Z\gamma_t)\\, correcting
for the classification error of the Step-2 class assignments. The
measurement model is held fixed, so the covariates cannot change the
classes.

## Usage

``` r
tse_covariate(
  object,
  formula,
  method = c("ML", "BCH", "none"),
  se = c("corrected", "robust"),
  ref = 1,
  start = NULL,
  control = NULL,
  data = NULL
)
```

## Arguments

- object:

  A classification from
  [`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md).

- formula:

  One-sided formula for the covariates, e.g. `~ age + sex`. Factors,
  interactions, and transformations are allowed.

- method:

  Step-3 estimator: `"ML"`, `"BCH"`, or `"none"` (see Details).

- se:

  Standard errors: `"corrected"` or `"robust"` (see Details).

- ref:

  Reference class of the multinomial logit: a class number or label such
  as `"C2"`.

- start:

  Optional starting values: a (Q+1) x (T-1) coefficient matrix (Q
  covariates plus the intercept, T classes). By default, the two-step
  estimates (Bakk and Kuha 2018) are used.

- control:

  Estimation settings; default: those of the measurement model. See
  [`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md).

- data:

  Optional data frame with the covariates: the classified data (the same
  rows, in the same order) with any additional columns. Omitted: the
  data stored in `object`.

## Value

A `tseLCA_covariate` object; see
[`coef.tseLCA_structural()`](https://samleebyu.github.io/tseLCA/reference/coef.tseLCA_structural.md),
[`summary.tseLCA_structural()`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md),
[`predict.tseLCA_covariate()`](https://samleebyu.github.io/tseLCA/reference/predict.tseLCA_covariate.md),
and
[`anova.tseLCA_covariate()`](https://samleebyu.github.io/tseLCA/reference/anova.tseLCA_covariate.md).
Pass it to
[`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)
to also model a distal outcome.

## Details

Estimators (`method`):

- `"ML"` (default): the bias-adjusted maximum likelihood estimator of
  Vermunt (2010), which treats the assigned class as an indicator of the
  true class with known classification-error probabilities.

- `"BCH"`: the Bolck-Croon-Hagenaars estimator (Bolck, Croon, and
  Hagenaars 2004; Vermunt 2010), which reweights the assignments with
  the inverse of the classification-error matrix.

- `"none"`: the uncorrected three-step estimator (a weighted multinomial
  logit of the assigned classes), which is biased toward zero; provided
  for comparison.

Standard errors (`se`): `"corrected"` (default) adds the uncertainty of
the Step-1 measurement model (Bakk, Oberski, and Vermunt 2014) to the
robust (sandwich) Step-3 variance; `"robust"` omits it. For `"BCH"` the
robust variance is used, which accounts for the Step-1 uncertainty
through the weights (Vermunt 2010); for `"none"` the robust variance is
used.

## References

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, 22(4), 520–540.
[doi:10.1093/pan/mpu003](https://doi.org/10.1093/pan/mpu003)

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, 12(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
cl <- tse_classify(m)
fit <- tse_covariate(cl, ~ Zp)
summary(fit)
#> Three-step latent class model: covariates
#>   Classes: 3   Estimator: ML   N: 500
#>   Log-lik: -1339.0650 (df = 22)   AIC: 2722.13   BIC: 2814.85
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
confint(fit)
#>                     2.5 %     97.5 %
#> (Intercept):C2  1.4013837  2.6807691
#> Zp:C2          -1.1623650 -0.6017952
#> (Intercept):C3 -4.7008849 -2.2662383
#> Zp:C3           0.6073547  1.1896410
predict(fit, newdata = data.frame(Zp = 1:5))
#>          C1         C2         C3
#> 1 0.2346248 0.74768656 0.01768866
#> 2 0.3993275 0.52673529 0.07393719
#> 3 0.4998232 0.27289595 0.22728081
#> 4 0.4268484 0.09646544 0.47668621
#> 5 0.2606745 0.02438452 0.71494097

# BCH, and the uncorrected estimator for comparison
coef(tse_covariate(cl, ~ Zp, method = "BCH"))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.0073032     -0.8634467     -3.3018545      0.8558238 
coef(tse_covariate(cl, ~ Zp, method = "none"))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      1.6510780     -0.7017868     -3.0087390      0.7850891 
```
