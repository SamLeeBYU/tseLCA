# Compress a one-hot expanded Y matrix back to integer codes

Inverse of `expand_Y`. Takes a one-hot expanded matrix where each item
occupies `K_h` consecutive columns (one per category, 0-based) and
returns an N x H integer matrix of category codes (0, 1, ..., K_h-1).

## Usage

``` r
compress_Y(mY_exp, ivItemcat)
```

## Arguments

- mY_exp:

  N x sum(K_h) one-hot matrix (as stored in `fit0$mU`).

- ivItemcat:

  Integer vector of category counts per item (length H).

## Value

N x H integer matrix of category codes.

## Details

Rows where all K_h columns for an item are `NA` are returned as `NA` for
that item.
