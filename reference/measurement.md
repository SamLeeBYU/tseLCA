# Components of a fitted tseLCA model

Extract the step-wise components of a model fitted with
[`tseLCA()`](https://samleebyu.github.io/tseLCA/reference/tseLCA.md) or
the step-wise functions: the Step-1 measurement model, the Step-2
classification, and the Step-3 covariate and distal outcome models.

## Usage

``` r
measurement(x, ...)

# S3 method for class 'tseLCA'
measurement(x, ...)

classification(x, ...)

# S3 method for class 'tseLCA'
classification(x, ...)

covariate(x, ...)

# S3 method for class 'tseLCA'
covariate(x, ...)

distal(x, ...)

# S3 method for class 'tseLCA'
distal(x, ...)
```

## Arguments

- x:

  A fitted tseLCA object.

- ...:

  Unused.

## Value

`measurement()`: a `tseLCA_measurement` object. `classification()`: a
`tseLCA_classify` object. `covariate()`: a `tseLCA_covariate` object.
`distal()`: a `tseLCA_distal` object.

## Examples

``` r
d <- generate_data(500, "high", "covariate", seed = 1)
d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)
fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
class_sizes(measurement(fit))
#>       C1       C2       C3 
#> 0.356968 0.330817 0.312215 
classification(fit)$D
#>             W=C1       W=C2        W=C3
#> X=C1 0.963004187 0.03202622 0.004969594
#> X=C2 0.039681675 0.94195543 0.018362900
#> X=C3 0.003557478 0.02867658 0.967765940
coef(covariate(fit))
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.0410764     -0.8820801     -3.4835616      0.8984978 
omnibus_test(distal(fit))
#> 
#>  Wald test of equal distal outcome distributions across latent classes
#> 
#> data:  gaussian distal outcome, 3 classes
#> W = 319, df = 2, p-value < 2.2e-16
#> 
```
