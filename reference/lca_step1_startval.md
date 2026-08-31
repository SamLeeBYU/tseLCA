# Fit the LCA measurement model from an externally supplied classification

A thin wrapper around multilevLCA's deterministic initialization path.
[`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)'s
default Step-1 initialization (k-means on principal components) is
deterministic given the data and, on some datasets, converges to a local
rather than global optimum of the Step-1 log-likelihood. If you have
already found a better solution with an external solver run with many
random starts (e.g. StepMix, poLCA, or similar), this function lets you
inject that classification directly: it writes `startval` into a
temporary column of `data` and calls
`multiLCA(..., startval = <that column>, kmea = FALSE)`, which skips
k-means entirely and initializes the EM algorithm from the supplied
classification instead.

## Usage

``` r
lca_step1_startval(
  data,
  Y.names,
  n_classes,
  startval,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-08,
  incomplete = FALSE,
  rebase = "C1",
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame containing at minimum the indicator columns.

- Y.names:

  Character vector of item column names.

- n_classes:

  Integer. Number of latent classes.

- startval:

  Either of the following, giving a starting classification for the
  measurement model:

  An integer vector

  :   Length `nrow(data)`, a starting class assignment (`1..n_classes`)
      for every row of `data`, typically obtained from an external
      latent class solver run with many random starts (e.g. the modal
      class from many-random-start posterior probabilities, as in
      StepMix or poLCA).

  A numeric matrix

  :   A conditional item-response probability matrix \\P(Y_h = k \mid X
      = t)\\ with one row per (item, category) pair – items in `Y.names`
      order, categories `0..K_h-1` within each item, matching the column
      order of `expand_Y(data[, Y.names], ivItemcat)` – and one column
      per class. A per-row classification is derived internally by
      naive-Bayes argmax under a flat class prior (see
      `classify_from_phi()`). This is the natural format for an
      externally estimated Step-1 solution that isn't tied to this
      specific sample, e.g. poLCA's `probs` output or a published
      item-response table.

  No automatic random restarts are performed on top of this starting
  value (contrast
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)'s
  `iter.measurement`/`R2.threshold` restart logic, which applies only to
  multilevLCA's own k-means initialization, and its `n_init` argument,
  which does run unconditional random restarts but from independent
  random classifications rather than a single fixed one).

- maxIter.measurement:

  Maximum EM iterations before giving up on convergence. Default
  `5000L`.

- measurement.tol:

  Convergence tolerance. Default `1e-8`.

- incomplete:

  Logical. FIML for partially missing indicators. See the `Missing Data`
  section of `vignette("tseLCA", package = "tseLCA")`. Default `FALSE`.

- rebase:

  Character or integer specifying the reference latent class. Use
  `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
  measurement model is permuted so this class becomes column 1, making
  it the reference for all downstream multinomial logit
  parameterizations.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A list with `$fit0`
([`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
measurement model, rebase-permuted per `rebase`) and `$fitZ = NULL`.
This is the same shape as
[`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)'s
return value, and matches
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)'s
`$measurement_model` when `startval` is passed there directly.

## Details

Most users should not need to call this function directly. Pass
`startval` to
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
(for structural estimation) or
[`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
(for a measurement-only fit) instead – both implement the same mechanism
and return the fitted measurement model as `$measurement_model$fit0` /
`$fit0` respectively. This function is documented mainly to describe
what `startval` accepts and how it is used internally.

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Recommended: pass `startval` to three_step() (or lca_step1() for a
# measurement-only fit) rather than calling this function directly --
# both use this same mechanism internally.

# A starting classification from an external solver (here, the DGP's own
# true classes, standing in for e.g. a StepMix solution with many
# random starts):
fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
                  startval = d$X)
fit$measurement_model$fit0$vPi
#>                
#> P(C1) 0.3495135
#> P(C2) 0.2915225
#> P(C3) 0.3589640

# Equivalently, supply a conditional item-response probability matrix
# (one row per item-category pair, in Y.names order -- since all 6 items
# here are binary, each contributes 2 rows: P(Y=0|C), P(Y=1|C)). In
# practice this would come from an external solver (e.g. poLCA's `probs`
# or a published item-response table); here a quick first-pass fit
# stands in for that external source.
fit_ref <- three_step(d, paste0("Y", 1:6), n_classes = 3)$measurement_model$fit0
phi <- matrix(0, nrow = 12, ncol = 3)
for (h in 1:6) {
  phi[2 * h - 1, ] <- 1 - fit_ref$mPhi[h, ] # P(Y_h = 0 | C)
  phi[2 * h,     ] <- fit_ref$mPhi[h, ]     # P(Y_h = 1 | C)
}
fit_phi <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
                      startval = phi)
# }
```
