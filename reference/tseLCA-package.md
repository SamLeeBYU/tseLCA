# tseLCA: Three-Step Estimation for Latent Class Analysis

**tseLCA** implements bias-adjusted three-step estimators for structural
latent class models with covariates and distal outcomes. Building on the
efficient measurement-model estimation procedures in multilevLCA,
**tseLCA** extends existing functionality through modern three-step
estimators, classification-error corrections, and analytic sandwich
variance estimation that propagates uncertainty from the latent class
measurement stage through to the final structural parameter estimates.

## The three-step approach

Three-step methods separate the latent class model into three stages,
preserving the measurement structure established during class formation.
In contrast to one-step approaches (e.g. poLCA), where including
covariates may alter the underlying latent class definitions, three-step
methods fix the measurement model before estimating structural
relationships and adjust for classification error at the final stage.

- Step 1:

  Estimate the LCA measurement model – class-conditional item-response
  probabilities and class prevalences estimated with multilevLCA.

- Step 2:

  Assign posterior class probabilities and compute the T x T
  misclassification matrix.

- Step 3:

  Estimate the structural model using bias-adjusted weights (ML or BCH
  correction) with sandwich variance estimation.

## Main functions

- [`three_step`](https://samleebyu.github.io/tseLCA/reference/three_step.md):

  Full three-step estimation pipeline. Accepts covariates (`Zp.names`),
  distal outcomes (`Zo.name`), or both. Handles measurement,
  classification, and structural estimation in a single call.

- [`lca_step1`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md):

  Step-1 measurement model via multilevLCA. Also computes two-step
  covariate initialisation via
  [`fitZ_from_fit0`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md)
  when `Zp.names` is supplied.

- [`fitZ_from_fit0`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md):

  Two-step covariate estimation by fixing measurement parameters at
  their Step-1 values and estimating multinomial logit coefficients via
  EM. Returns starting values for Step 3.

- [`fitZ_from_multiLCA`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md):

  Two-step covariate estimation via `multiLCA(fixedpars = 1)`, returning
  multilevLCA's bias-corrected standard errors. Called when
  `get.twostep.vcov = TRUE` in
  [`three_step`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- [`generate_data`](https://samleebyu.github.io/tseLCA/reference/generate_data.md):

  Simulate data from the Bakk & Kuha (2018) three-class LCA design with
  covariate or distal outcome scenarios across three separation levels.

- [`generate_all_conditions`](https://samleebyu.github.io/tseLCA/reference/generate_all_conditions.md):

  Batch simulation across all scenarios, separation levels, and sample
  sizes for replication studies.

## Estimators

- ML (default):

  The Vermunt (2010) ML correction uses weighted pseudo-likelihood with
  the misclassification matrix as a bias adjustment. Set
  `use.bch = FALSE` (default).

- BCH:

  The Bolck, Croon & Hagenaars (2004) correction inverts the
  misclassification matrix to obtain direct class weights. Set
  `use.bch = TRUE`. Works well under high separation but may produce an
  ill-conditioned Hessian with low separation.

## Variance estimation

- Simple (robust):

  Sandwich SEs from Step 3 only (`use.simple.cov = TRUE`). Efficient
  when separation is high and measurement uncertainty is negligible.

- Full correction:

  Analytic propagation of measurement uncertainty from Step 1 through
  the classification-error correction (`use.simple.cov = FALSE`,
  default). Uses an analytic Jacobian computed from soft (proportional)
  posteriors regardless of the `use.modal.assignment` setting.

## Supported features

- Binary and polytomous indicators (following multilevLCA coding
  conventions: consecutive integers from 0).

- Gaussian, Poisson, and binomial distal outcome families.

- Full-information maximum likelihood (FIML) for partially observed
  indicator patterns (`incomplete = TRUE`).

- Flexible measurement and structural samples: the measurement model can
  be estimated on a reference sample and applied to a different analysis
  sample via the `step1` argument.

- Arbitrary reference class via the `rebase` argument.

- S3 methods: `print`, `summary`, `coef`, `vcov`, `plot` for all four
  return subclasses (`tseLCA_measurement`, `tseLCA_covariate`,
  `tseLCA_distal`, `tseLCA_both`).

## Getting started

    # Install from GitHub
    pak::pak("SamLeeBYU/tseLCA")

    # Introductory vignette
    vignette("tseLCA-workflow", package = "tseLCA")

## References

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Bakk, Z., Pohle, M. J., & Kuha, J. (2025). Bias-adjusted three-step
estimation of structural models for latent classes. *Multivariate
Behavioral Research*.
[doi:10.1080/00273171.2025.2473935](https://doi.org/10.1080/00273171.2025.2473935)

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

Sam Lee <samlee@arizona.edu>
