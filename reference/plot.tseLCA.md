# Plot item-response probability profiles for a tseLCA model

Delegates to `plot.multiLCA` from multilevLCA, which draws the
class-specific item-response probability profiles of the Step-1
measurement model.

## Usage

``` r
# S3 method for class 'tseLCA'
plot(x, horiz = FALSE, clab = NULL, ...)
```

## Arguments

- x:

  A fitted `tseLCA` object.

- horiz:

  Logical. If `TRUE`, item labels are drawn horizontally.

- clab:

  Optional character vector of length T giving class labels.

- ...:

  Further arguments passed to `plot.multiLCA`.

## Value

Called for its side effect (a base-graphics plot). Invisibly returns
`NULL`.

## Examples

``` r
d     <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
plot(fit_m)

plot(fit_m, clab = c("Low", "Mixed", "High"))
```
