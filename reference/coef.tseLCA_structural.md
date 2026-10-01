# Coefficients of a fitted tseLCA model

Returns a named coefficient vector whose names match the rows and
columns of [`vcov()`](https://rdrr.io/r/stats/vcov.html), so that
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) and other
generic tools work.

## Usage

``` r
# S3 method for class 'tseLCA_structural'
coef(
  object,
  component = c("all", "covariate", "distal"),
  step = c("three_step", "two_step"),
  matrix = FALSE,
  ...
)

# S3 method for class 'tseLCA_measurement'
coef(object, ...)
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

- matrix:

  Logical. If `TRUE`, return the coefficients in their natural matrix
  layout: (Q+1) x (T-1) for covariate models, T x C for multinomial
  distal outcomes, a list of both for `tseLCA_both`.

- ...:

  Further arguments (currently unused).

## Value

A named numeric vector (or matrix / list if `matrix = TRUE`).

## Details

- Measurement models: class-size log-ratios \\\log(\pi_t/\pi_1)\\ and
  item-response log-ratios \\\log(P(Y=k \mid t)/P(Y=0 \mid t))\\. Use
  [`class_sizes()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
  and
  [`item_probs()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
  for the probability scale.

- Covariate models: multinomial-logit coefficients named
  `"<covariate>:<class>"`, e.g. `"Zp:C2"`.

- Distal-outcome models: class-specific means (`gaussian`), log means
  (`poisson`), logits (`binomial`), or category probabilities
  (`multinomial`, named `"<class>:<category>"`).

## Examples

``` r
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
#> Warning: three_step() is deprecated as of tseLCA 2.0.0; use tseLCA() or the step-wise tse_lca(), tse_classify(), tse_covariate(), and tse_distal(). (Shown once per session; see NEWS for the new interface.)
coef(fit)
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.2333837     -1.1569878     -3.2742157      0.9400712 
coef(fit, matrix = TRUE)
#>                    C2         C3
#> (Intercept)  2.233384 -3.2742157
#> Zp          -1.156988  0.9400712
confint(fit)
#>                     2.5 %     97.5 %
#> (Intercept):C2  1.0068081  3.4599593
#> Zp:C2          -1.7453032 -0.5686724
#> (Intercept):C3 -4.6837138 -1.8647177
#> Zp:C3           0.5684972  1.3116452
```
