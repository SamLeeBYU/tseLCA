# Extract Y.exp, mDesign, posteriors from a multilevLCA mU matrix

`fit0$mU` from multilevLCA stores data already in one-hot expanded form:
each item k occupies R_k consecutive columns (one per category),
followed by T columns of posterior class probabilities.

## Usage

``` r
extract_Y_from_mU(fit0, ivItemcat = NULL)
```

## Arguments

- fit0:

  Raw multilevLCA fit object with `$mU`, `$mPhi`, `$vPi`.

- ivItemcat:

  Integer vector of category counts per item (length K). If `NULL`,
  inferred from `fit0$mPhi` dimensions.

## Value

A list with:

- Y.exp:

  n x sum(R_k) expanded one-hot matrix (NAs replaced with 0).

- mDesign:

  n x sum(R_k) design/mask matrix. `NULL` if no missing.

- ivItemcat:

  Integer vector of category counts per item.

- u_post:

  n x T posterior class probability matrix from `mU`.

## Details

For dichotomous items (R_k=2) the two columns are stored. For polytomous
items (R_k\>2) all R_k columns are stored. This function first
compresses the expanded Y back to integer codes, then re-expands
consistently with `expand_Y` so downstream functions receive the correct
n x sum(R_k) matrix.
