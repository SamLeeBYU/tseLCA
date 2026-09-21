## Test environments
* local machine (Windows, WSL2), R 4.6.1
* win-builder R-devel (x86_64-w64-mingw32, 2026-09-15 r90540)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are no downstream dependencies for this package.

## Changes since last CRAN version (1.0.2)

### Bug fix — BCH weight matrix orientation (breaking numerical change)

The BCH weight matrix was computed with the wrong orientation in all
`use.bch = TRUE` calls. The corrected formula follows Mplus Web Note 21:
`H^{-1}` is taken row-wise for each case's modal class, where
`H = t(pwx)` is row-stochastic. The previous code used the column-stochastic
`pwx` directly, which only agrees with the correct form when `pwx` is
symmetric. Numerical output from `use.bch = TRUE` will differ from 1.0.2
in any dataset with asymmetric classification errors; this is intentional
and correct behavior.

### New features

- `lca_step1_startval()`: wrapper around `multilevLCA::multiLCA()` that
  injects a user-supplied Step-1 starting value (integer classification
  vector or item-response probability matrix), bypassing the default
  k-means initialization. Allows users to pass an external solution from
  `poLCA`, `StepMix`, or a published item-response table directly.
- `startval` argument added to `lca_step1()`, `fitZ_from_multiLCA()`, and
  `three_step()`.
- `n_init` argument added to `lca_step1()`, `fitZ_from_multiLCA()`, and
  `three_step()`: fits the measurement model `n_init` times from independent
  uniform-random starts and keeps the highest-log-likelihood result (the
  unconditional multi-start analog of `nrep` in `poLCA`).
- `family = "multinomial"` added to `three_step()` for nominal categorical
  distal outcomes; estimates the full `T x C` class-conditional probability
  matrix with correct sandwich variance.
- `omnibus_test()`: generalized Wald test of `H0: theta_1 = ... = theta_T`
  for any `tseLCA_distal` or `tseLCA_both` object. Uses Moore-Penrose
  pseudoinverse so it remains valid for the rank-deficient multinomial case;
  degrees of freedom recover `(T-1)*(C-1)` for multinomial and `T-1` for
  scalar-parameter families.
