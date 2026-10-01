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

- missing, control:

  As in
  [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

## Value

A `tseLCA_measurement` object, usable like one from
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md).

## Details

The Step-1 variance used for corrected standard errors in Step 3 is
computed on `data` at the given parameters, which is valid when they are
the maximum likelihood estimates for `data` (e.g. a model estimated on
these data and saved). For parameters estimated on another sample, use
`se = "robust"` in Step 3, or refit with
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
using the parameters as `start`.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)

# the same measurement model from its parameters
m2 <- as_tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d,
                 class_sizes = class_sizes(m), item_probs = item_probs(m))
all.equal(logLik(m2), logLik(m), tolerance = 1e-6)
#> [1] TRUE
coef(tse_covariate(tse_classify(m2), ~ Zp))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.0410764     -0.8820801     -3.4835616      0.8984978 
```
