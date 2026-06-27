# Individual-level BHHH Varmat matching teLCA's lca_em.cpp

Computes the standard outer-product (BHHH) sandwich variance on the N
individual observations, exactly replicating what LCA_teLCA returns in
\$Varmat and \$mScore.

## Usage

``` r
lca_indiv_varmat(mY, T, pi, phi)
```

## Arguments

- mY:

  N x H raw 0/1 matrix (original, not one-hot expanded)

- T:

  Number of classes

- pi:

  Length-T class prevalences (at converged estimates)

- phi:

  H x T item-response probabilities (at converged estimates) NOTE: pass
  the n_free x T matrix from fit0\$mPhi, i.e. one row per item
  (P(Y=1\|t) for binary items)

## Value

list with \$Infomat, \$Varmat, \$SEs, \$mScore (N x p score matrix)

## Details

Infomat = S'S / N, Varmat = psinv(Infomat) / N = \\(S'S)^{-1}\\
