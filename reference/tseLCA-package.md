# tseLCA: Three-Step Estimation for Latent Class Analysis

tseLCA relates latent classes to covariates and distal outcomes by
bias-adjusted three-step estimation. The latent class measurement model
is estimated first and held fixed, so the structural variables cannot
change the meaning of the classes; the structural estimates are
corrected for the classification error of the class assignments (BCH and
ML estimators), and their standard errors account for the uncertainty of
the measurement model. Measurement models are estimated with multilevLCA
(Lyrvall et al., 2025).

## The three steps

1.  **Measurement model**
    ([`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)):
    class sizes and class-conditional item-response probabilities,
    estimated from the indicators alone. With several numbers of
    classes, a class-enumeration table (AIC, BIC, SABIC, entropy) for
    choosing the number of classes.

2.  **Classification**
    ([`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)):
    posterior class probabilities, modal or proportional class
    assignments, and the classification-error probabilities \\P(W = s
    \mid X = t)\\.

3.  **Structural model**: a multinomial logit of class membership on
    covariates
    ([`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md)),
    and/or class-specific distributions of a distal outcome
    ([`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)),
    with the ML (Vermunt 2010; Bakk, Tekle & Vermunt 2013) or BCH
    (Bolck, Croon & Hagenaars 2004) correction.

[`tseLCA()`](https://samleebyu.github.io/tseLCA/reference/tseLCA.md)
runs all three steps from one formula,
`indicators ~ covariates | distal outcome`;
[`measurement()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`classification()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`covariate()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
and
[`distal()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)
extract the components of a fitted model.

## Estimators and standard errors

- `method = "ML"` (default):

  Vermunt's (2010) maximum likelihood correction, treating the assigned
  class as an indicator of the true class with known
  classification-error probabilities.

- `method = "BCH"`:

  The Bolck-Croon-Hagenaars correction, reweighting the assignments by
  the inverse classification-error matrix. Reliable when classes are
  well separated.

- `method = "none"`:

  The uncorrected three-step estimator, for comparison.

- `se = "corrected"` (default):

  Sandwich standard errors plus the propagated uncertainty of the Step-1
  measurement model (Bakk, Oberski & Vermunt 2014), and, for combined
  models, of the covariate model.

- `se = "robust"`:

  Sandwich standard errors of Step 3 only.

The two-step estimator of Bakk & Kuha (2018) is available with
[`tse_twostep()`](https://samleebyu.github.io/tseLCA/reference/tse_twostep.md).

## Features

- Binary and polytomous indicators, coded as factors, logicals,
  characters, or numbers; full-information maximum likelihood for
  missing indicator values (`missing = "fiml"`).

- Covariate formulas with factors, interactions, and transformations;
  Wald tests by term
  ([`anova.tseLCA_covariate()`](https://samleebyu.github.io/tseLCA/reference/anova.tseLCA_covariate.md));
  predicted class probabilities
  ([`predict.tseLCA_covariate()`](https://samleebyu.github.io/tseLCA/reference/predict.tseLCA_covariate.md));
  any reference class.

- Gaussian, Poisson, binomial, and multinomial distal outcomes, with an
  omnibus test of equality across classes
  ([`omnibus_test()`](https://samleebyu.github.io/tseLCA/reference/omnibus_test.md)).

- Measurement models estimated on one sample and applied to another
  ([`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)
  with `newdata`).

- Standard methods for fitted models:
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`update()`](https://rdrr.io/r/stats/update.html).

- Simulation from the design of Bakk & Kuha (2018)
  ([`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md)).

The 1.x function
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
is deprecated; its help page maps each of its arguments to the current
interface.

## Getting started

    vignette("tseLCA-workflow", package = "tseLCA")

## References

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, 22(4), 520–540.
<https://www.jstor.org/stable/24573086>

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, 12(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025).
Multilevel latent class analysis: State-of-the-art methodologies and
their implementation in the R package multilevLCA. *Multivariate
Behavioral Research*, 60(4), 731–747.
[doi:10.1080/00273171.2025.2473935](https://doi.org/10.1080/00273171.2025.2473935)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

## See also

Useful links:

- <https://samleebyu.github.io/tseLCA/>

- <https://github.com/SamLeeBYU/tseLCA>

- Report bugs at <https://github.com/SamLeeBYU/tseLCA/issues>

## Author

Sam Lee <samlee@arizona.edu>, Jay Goodliffe <goodliffe@byu.edu>
