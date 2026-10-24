# Use a measurement model with given parameters (Step 1)

Creates a measurement model from given class sizes and item-response
probabilities, evaluated on `data`, without estimating it. This allows
Steps 2 and 3 to be based on a measurement model estimated elsewhere: in
another program, reported in a publication, or saved from an earlier
analysis.

## Usage

``` r
as_tse_lca(
  formula,
  data,
  class_sizes,
  item_probs,
  vcov = NULL,
  missing = c("listwise", "fiml"),
  control = tse_control()
)
```

## Arguments

- formula:

  `cbind(Y1, Y2, ...) ~ 1`, as in
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

- data:

  A data frame.

- class_sizes:

  Class proportions, one per class (they are normalized to sum to one).

- item_probs:

  Item-response probabilities in the layout of
  [`item_probs()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md):
  one column per class, and one row per binary item (\\P(Y = 1 \mid X =
  t)\\, where 1 is the item's second category) or per category of a
  polytomous item (\\P(Y = k \mid X = t)\\), in the order of the
  indicators.

- vcov:

  Optional variance matrix of the Step-1 parameters, for corrected
  standard errors in Step 3. It refers to the unconstrained log-ratio
  parameters, in the layout of
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) of a
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
  model: the class-size log-ratios \\\log(\pi_t / \pi_1)\\, \\t = 2,
  \dots, T\\, then, class by class and item by item, \\\log(P(Y = k \mid
  X = t) / P(Y = 0 \mid X = t))\\ for each category \\k \ge 1\\. Classes
  and items are in the order of `class_sizes` and of the indicators. Row
  and column names are ignored.

- missing, control:

  As in
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

## Value

A `tseLCA_measurement` object, usable like one from
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

## Details

Corrected standard errors in Step 3 (`se = "corrected"`) need the
variance matrix of the Step-1 parameters, which can only come from the
sample they were estimated on. A model from `as_tse_lca()` receives only
the parameter values, so it uses the variance matrix given in `vcov`.
Without `vcov`, Step 3 warns and reports robust standard errors
(`se = "robust"`), which do not include the Step-1 uncertainty.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) of the model returns
`vcov`, or an error if none was given.

If the parameters are the maximum likelihood estimates for `data` (a
model estimated on these data and saved), `vcov(m)` of the original
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
model is the matching `vcov`. For parameters estimated on another
sample, take the variance matrix from that estimation (e.g. Latent GOLD
or Mplus output), or use
[`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)
with `newdata` on the original model.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)

# the same measurement model from its parameters and their variance
m2 <- as_tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d,
                 class_sizes = class_sizes(m), item_probs = item_probs(m),
                 vcov = vcov(m))
all.equal(logLik(m2), logLik(m), tolerance = 1e-6)
#> [1] TRUE
coef(tse_covariate(tse_classify(m2), ~ Zp))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.0410764     -0.8820801     -3.4835616      0.8984978 

# without `vcov`, Step 3 warns and reports robust standard errors
m3 <- as_tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d,
                 class_sizes = class_sizes(m), item_probs = item_probs(m))
fit3 <- tse_covariate(tse_classify(m3), ~ Zp)
#> Warning: The measurement model was created by as_tse_lca() without `vcov`, so the standard errors cannot be corrected for the Step-1 uncertainty; robust standard errors are reported. Supply the variance matrix of the Step-1 parameters with as_tse_lca(vcov = ).
```
