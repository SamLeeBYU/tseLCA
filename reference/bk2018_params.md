# Default population parameters for the Bakk & Kuha (2018) simulation

A list of pre-specified parameters used by
[`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md)
and related functions. All elements correspond to the values stated in
the paper (Section 3, p. 879).

## Usage

``` r
bk2018_params
```

## Details

- `class_props`:

  Length-3 vector of equal class proportions (1/3 each).

- `separation_levels`:

  Named vector mapping `"low"`, `"mid"`, `"high"` to the probability of
  a "likely" response (0.70, 0.80, 0.90).

- `covariate_params`:

  List with `$b0` (intercepts) and `$b` (slopes) for the multinomial
  logit P(X=t \| Zp). Intercepts `b02` and `b03` are set so that
  marginal class sizes average to 1/3 when Zp ~ Uniform{1..5}.

- `distal_params`:

  List with `$mu` (class means, c(-1, 1, 0)) and `$sigma` (residual SD
  = 1) for the distal outcome model.

## Examples

``` r
# True item-response probabilities for high separation
bk2018_params$rho_high
#> NULL

# Covariate model parameters (intercepts and slopes)
bk2018_params$covariate_params
#> $b0
#> [1]  0.000000  2.344591 -3.655370
#> 
#> $b
#> [1]  0 -1  1
#> 

# Distal outcome parameters
bk2018_params$distal_params
#> $mu
#> [1] -1  1  0
#> 
#> $sigma
#> [1] 1
#> 
```
