# Omnibus Wald test of class equality for a distal outcome

Tests whether the distal outcome distribution differs across latent
classes at all – \\H_0: \theta_1 = \theta_2 = \dots = \theta_T\\ for the
class-specific distal parameters \\\theta_t\\ (class means for
`family = "gaussian"`, log-rates for `"poisson"`, logits for
`"binomial"`, or the full length-C category-probability vector for
`"multinomial"`) – using a generalized Wald test with
[`vcov()`](https://rdrr.io/r/stats/vcov.html). This answers whether the
outcome's distribution is associated with class membership at all,
before drilling into which classes differ (a natural next step – not yet
implemented in tseLCA – is pairwise post-hoc comparisons with
multiplicity correction). The degrees of freedom equal the rank of the
contrast covariance (`T - 1` for scalar outcomes; `(T - 1) * (C - 1)`
for multinomial, i.e. the textbook chi-squared test of homogeneity in a
\\T \times C\\ table), computed with a Moore-Penrose pseudo-inverse so
the test remains valid despite `multinomial`'s inherently singular
covariance (each class's category probabilities sum to 1).

## Usage

``` r
omnibus_test(object, ...)

# S3 method for class 'tseLCA_distal'
omnibus_test(object, ...)

# S3 method for class 'tseLCA_both'
omnibus_test(object, ...)

# S3 method for class 'tseLCA_omnibus'
print(x, digits = 4, ...)
```

## Arguments

- object:

  A `tseLCA_distal` object, or a `tseLCA_both` object (tests its
  `$distal` component).

- ...:

  Unused; present for S3 method consistency.

- x:

  A `tseLCA_omnibus` object.

- digits:

  Integer. Number of decimal places for the test statistic.

## Value

An object of class `"tseLCA_omnibus"`: a list with `$statistic` (the
Wald chi-squared statistic), `$df`, and `$p.value`.

## Examples

``` r
# \donttest{
d <- generate_data(300, "high", "distal", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
omnibus_test(fit)
#> Omnibus Wald test of class equality (distal outcome)
#>   Family: gaussian   Classes: 3
#>   W(2) = 196.5943, p < 0.001
# }
```
