# Individual-level BHHH Varmat for binary and polytomous LCA

Computes the outer-product (BHHH) information matrix and
variance-covariance matrix for LCA measurement model parameters in the
unconstrained (logit/log-ratio) space, matching multilevLCA's `$Varmat`.

## Usage

``` r
lca_indiv_varmat(
  Y.exp,
  mDesign.exp,
  fit0,
  ivItemcat,
  boundary.tol = 0.01,
  use.freq = TRUE
)
```

## Arguments

- Y.exp:

  Expanded one-hot matrix (N x sum(K_h)).

- mDesign.exp:

  Expanded design matrix (same dims), or `NULL`.

- fit0:

  Step-1 fit with `$vPi` and `$mPhi`.

- ivItemcat:

  Integer vector of category counts per item.

- boundary.tol:

  Scalar tolerance for boundary detection. Default `1e-2`.

- use.freq:

  Collapse duplicate score rows before cross-product. Default `TRUE`.

## Value

A list containing `Infomat`, `Varmat`, `SEs`, and the individual score
matrix `mScore`. Individual-level BHHH Varmat for binary and polytomous
LCA

Computes the outer-product (BHHH) information matrix and
variance-covariance matrix for LCA measurement model parameters in the
unconstrained (logit/log-ratio) space.

Assumes `fit0$mPhi` has the following structure (from multilevLCA):

- Dichotomous item h (`ivItemcat[h] == 2`): 1 row = \\P(Y=1\|C)\\. The
  base level \\P(Y=0\|C)\\ is excluded.

- Polytomous item h (`ivItemcat[h] > 2`): `K_h` rows = \\P(Y=0\|C),
  \ldots, P(Y=K_h-1\|C)\\. The base level IS included.

`expand_Y` produces one-hot columns in the same order so that
`expand_Phi(fit0$mPhi, ivItemcat)` aligns column-wise with
`expand_Y(mY, ivItemcat)`.

Free (estimable) parameters per item:

- Dichotomous: the single row of `mPhi` (\\P(Y=1\|C)\\).

- Polytomous: rows 2..K_h of `mPhi` (\\P(Y=1\|C), \ldots\\). Row 1
  (\\P(Y=0\|C)\\) is the reference and is not a free parameter.

Boundary parameters (within `boundary.tol` of 0 or 1) are treated as
fixed: their score columns are zeroed so they do not contribute to the
information matrix.

List with `$Infomat`, `$Varmat`, `$SEs`, `$mScore`.

## Details

The score in unconstrained space is \\s\_{it} = u\_{it}(y_i - d_i \circ
p\_{it})\\, where \\d_i\\ is the missing-data design indicator matrix.
