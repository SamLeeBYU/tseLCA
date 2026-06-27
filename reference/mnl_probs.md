# Compute multinomial logistic class probabilities given a covariate

Evaluates P(X = t \| Zp) for each observation using a multinomial logit
with one covariate and class-specific intercepts and slopes.

## Usage

``` r
mnl_probs(Zp, params)
```

## Arguments

- Zp:

  Numeric vector of length n. Covariate values.

- params:

  List with elements `$b0` (length-T intercepts, reference = 0) and `$b`
  (length-T slopes, reference = 0). See
  `bk2018_params$covariate_params`.

## Value

An n x T matrix of class probabilities (rows sum to 1).

## Examples

``` r
# Class membership probabilities for Zp = 1..5
mnl_probs(1:5, bk2018_params$covariate_params)
#>           [,1]       [,2]       [,3]
#> [1,] 0.2037951 0.78188368 0.01432126
#> [2,] 0.3842557 0.54234331 0.07340104
#> [3,] 0.4905617 0.25471421 0.25472410
#> [4,] 0.3842487 0.07339685 0.54235449
#> [5,] 0.2037890 0.01432027 0.78189074
```
