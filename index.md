# tseLCA

## Overview

**tseLCA** (*Three-Step Estimation for Latent Class Analysis*) relates
latent classes to covariates and distal outcomes by bias-adjusted
three-step estimation.

1.  **Measurement model**
    ([`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)).
    The latent classes are estimated from the indicators alone, and the
    number of classes is chosen from a class-enumeration table.
2.  **Classification**
    ([`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)).
    Observations are assigned to classes, and the classification error
    is estimated.
3.  **Structural model**
    ([`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md),
    [`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)).
    The classes are related to covariates and/or distal outcomes, with
    the ML (Vermunt 2010) or BCH (Bolck, Croon & Hagenaars 2004)
    correction for classification error.

Because the measurement model is fixed before any structural variable
enters, covariates and outcomes cannot change what the classes mean.
This is the key difference from one-step estimation (e.g. **poLCA**),
where the class solution can shift with every change to the structural
model. The standard errors of the structural estimates account for the
uncertainty of the measurement model (Bakk, Oberski & Vermunt 2014).
Measurement models are estimated with **multilevLCA**.

## Features

- Class enumeration with AIC, BIC, SABIC, entropy, and class sizes.
- ML and BCH three-step estimators, with modal or proportional
  assignment. The uncorrected estimator is also available, for
  comparison.
- Standard errors corrected for the uncertainty of the measurement
  model.
- Covariate formulas with factors, interactions, and transformations.
  Wald tests by term, predicted class probabilities, and any reference
  class.
- Gaussian, Poisson, binomial, and multinomial distal outcomes, alone or
  combined with covariates, and an omnibus test of equality across
  classes.
- Measurement models estimated on one sample and applied to another.
- Indicators as factors, logicals, characters, or numeric codes.
  Full-information maximum likelihood for missing indicator values.
- Standard R methods throughout:
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`update()`](https://rdrr.io/r/stats/update.html).

## Installation

``` r

install.packages("tseLCA")

# development version
# install.packages("pak")
pak::pak("SamLeeBYU/tseLCA")
```

## Example

``` r

library(tseLCA)

# Simulated data: six binary indicators, a covariate Zp, and a distal outcome Zo
d <- generate_data(n = 1000, separation = "high", scenario = "covariate", seed = 1)
d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)

# Step 1: choose the number of classes from the measurement model
sel <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 1:4)
sel
m <- best_model(sel, criterion = "BIC")

# Step 2: classification
cl <- tse_classify(m)

# Step 3: covariate and distal outcome models
fc <- tse_covariate(cl, ~ Zp)
summary(fc)
fb <- tse_distal(fc, Zo ~ 1)
summary(fb)

# The same model in one call
fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
```

See the [introductory
vignette](https://SamLeeBYU.github.io/tseLCA/articles/tseLCA-workflow.html)
for the full workflow.

## Upgrading from tseLCA 1.x

[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
still works, with the same estimates, but is deprecated. Its help page,
and the vignette, map each of its arguments to the new functions.
Version 2.0.0 also fixes several estimation bugs; see
[NEWS](https://samleebyu.github.io/tseLCA/NEWS.md).
