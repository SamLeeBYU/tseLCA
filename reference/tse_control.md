# Estimation settings for tseLCA models

Collects the numerical settings of the three estimation steps. Pass the
result as the `control` argument of the model-fitting functions.

## Usage

``` r
tse_control(
  step1.maxit = 5000L,
  step1.tol = 1e-08,
  step1.restarts = 10L,
  step1.restart.R2 = 0.7,
  n_init = NULL,
  step3.maxit = 200L,
  step3.tol = 1e-06,
  boundary.tol = 0.01,
  hessian = c("observed", "opg"),
  verbose = FALSE
)
```

## Arguments

- step1.maxit:

  Maximum number of EM iterations for the Step-1 measurement model. A
  fit that reaches it is retried with twice as many.

- step1.tol:

  Convergence tolerance of the Step-1 EM algorithm (change in
  log-likelihood).

- step1.restarts:

  Number of additional random starts tried when the default Step-1 fit
  has an entropy R\\^2\\ below `step1.restart.R2`; the fit with the
  highest log-likelihood is kept.

- step1.restart.R2:

  Entropy R\\^2\\ threshold that triggers `step1.restarts`.

- n_init:

  Optional number of Step-1 fits from independent random
  classifications, bypassing the default k-means initialization; the fit
  with the highest log-likelihood is kept. Unlike `step1.restarts`,
  these always run. `NULL` (default) uses the default initialization.

- step3.maxit:

  Maximum number of iterations for the Step-3 (structural model) EM or
  Newton-Raphson algorithm.

- step3.tol:

  Convergence tolerance of the Step-3 algorithm.

- boundary.tol:

  Step-1 probabilities within this distance of 0 or 1 are treated as
  fixed when computing the Step-1 variance.

- hessian:

  Step-3 information matrix used for standard errors: `"observed"`
  (default), or `"opg"` for the outer product of the case-wise scores.
  `"opg"` applies to ML covariate models.

- verbose:

  Logical. Print progress and convergence messages.

## Value

A list of class `"tse_control"`.

## Examples

``` r
tse_control()
#> tseLCA estimation settings
#>                  value   
#> step1.maxit      5000    
#> step1.tol        1e-08   
#> step1.restarts   10      
#> step1.restart.R2 0.7     
#> n_init           NULL    
#> step3.maxit      200     
#> step3.tol        1e-06   
#> boundary.tol     0.01    
#> hessian          observed
#> verbose          FALSE   
tse_control(step1.maxit = 10000, n_init = 20)
#> tseLCA estimation settings
#>                  value   
#> step1.maxit      10000   
#> step1.tol        1e-08   
#> step1.restarts   10      
#> step1.restart.R2 0.7     
#> n_init           20      
#> step3.maxit      200     
#> step3.tol        1e-06   
#> boundary.tol     0.01    
#> hessian          observed
#> verbose          FALSE   
```
