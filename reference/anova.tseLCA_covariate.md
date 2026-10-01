# Wald tests of covariate terms

Tests, term by term, that all multinomial-logit coefficients of a
covariate term are zero (for all non-reference classes), using the
model's variance matrix.

## Usage

``` r
# S3 method for class 'tseLCA_covariate'
anova(object, ...)
```

## Arguments

- object:

  A `tseLCA_covariate` object from
  [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md).

- ...:

  Unused.

## Value

An `anova` table with the Wald statistic, degrees of freedom, and
p-value of each term.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
anova(tse_covariate(tse_classify(m), ~ Zp))
#> Wald tests of covariate terms (all class contrasts)
#> 
#>    Df  Chisq Pr(>Chisq)    
#> Zp  2 70.768  4.294e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
