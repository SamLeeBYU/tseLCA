# Class sizes and item-response probabilities of the measurement model

`class_sizes()` returns the estimated class proportions and
`item_probs()` the class-conditional item-response probabilities of the
Step-1 measurement model underlying any fitted `tseLCA` object.

## Usage

``` r
class_sizes(object, ...)

# S3 method for class 'tseLCA'
class_sizes(object, se = FALSE, ...)

item_probs(object, ...)

# S3 method for class 'tseLCA'
item_probs(object, se = FALSE, ...)
```

## Arguments

- object:

  A fitted `tseLCA` object.

- ...:

  Further arguments (currently unused).

- se:

  Logical. If `TRUE`, also return standard errors.

## Value

`class_sizes()`: a named numeric vector of length T summing to one.
`item_probs()`: a matrix with one row per item (binary items: \\P(Y = 1
\mid X = t)\\) or per item category (polytomous items: \\P(Y = k \mid X
= t)\\) and one column per class. With `se = TRUE`, a list with elements
`estimate` and `se` of that form.

## Details

With `se = TRUE`, their standard errors are returned as well. They are
obtained by the delta method from the variance of the measurement
model's log-ratio parameters
([`vcov()`](https://rdrr.io/r/stats/vcov.html) of the measurement
model): class sizes are the softmax of \\\log(\pi_t / \pi_1)\\, and the
response probabilities of an item in class \\t\\ the softmax of
\\\log(P(Y = k \mid X = t) / P(Y = 0 \mid X = t))\\. Parameters on the
boundary of the parameter space are treated as fixed and get a standard
error of zero.

## Examples

``` r
d <- generate_data(200, "high", "covariate", seed = 1)
m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
class_sizes(m)
#>        C1        C2        C3 
#> 0.3495138 0.2915216 0.3589645 
item_probs(m)
#>                C1         C2         C3
#> P(Y1|C) 0.8702096 0.79456644 0.12317767
#> P(Y2|C) 0.9016604 0.88525528 0.10247858
#> P(Y3|C) 0.8743309 0.87570434 0.06720021
#> P(Y4|C) 0.8565891 0.09127798 0.06686104
#> P(Y5|C) 0.8909744 0.09780804 0.02807791
#> P(Y6|C) 0.8206322 0.13853263 0.09135284
item_probs(m, se = TRUE)$se
#>                 C1         C2         C3
#> P(Y1|C) 0.04554151 0.06597822 0.04799340
#> P(Y2|C) 0.03714323 0.05935448 0.04764421
#> P(Y3|C) 0.04439277 0.06221734 0.04390390
#> P(Y4|C) 0.05534999 0.05443704 0.03138966
#> P(Y5|C) 0.05008769 0.05959757 0.02279529
#> P(Y6|C) 0.05707720 0.05842492 0.03747299
```
