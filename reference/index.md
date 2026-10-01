# Package index

## One-call interface

Fit a complete three-step model from one formula.

- [`tseLCA()`](https://samleebyu.github.io/tseLCA/reference/tseLCA.md) :
  Three-step latent class analysis in one call
- [`measurement()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)
  [`classification()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)
  [`covariate()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)
  [`distal()`](https://samleebyu.github.io/tseLCA/reference/measurement.md)
  : Components of a fitted tseLCA model

## Step 1: measurement model

Estimate the latent classes and choose the number of classes.

- [`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
  : Fit a latent class measurement model (Step 1)
- [`as_tse_lca()`](https://samleebyu.github.io/tseLCA/reference/as_tse_lca.md)
  : Use a measurement model with given parameters (Step 1)
- [`best_model()`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
  [`` `[[`( ``*`<tseLCA_select>`*`)`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
  [`as.data.frame(`*`<tseLCA_select>`*`)`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
  [`print(`*`<tseLCA_select>`*`)`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
  [`plot(`*`<tseLCA_select>`*`)`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
  : Class enumeration results
- [`predict(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/predict.tseLCA_measurement.md)
  [`fitted(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/predict.tseLCA_measurement.md)
  : Class membership predictions from a measurement model
- [`class_sizes()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
  [`item_probs()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
  : Class sizes and item-response probabilities of the measurement model

## Step 2: classification

Assign observations to classes and estimate the classification error.

- [`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)
  [`print(`*`<tseLCA_classify>`*`)`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)
  : Assign observations to latent classes (Step 2)
- [`posterior()`](https://samleebyu.github.io/tseLCA/reference/posterior.md)
  [`classes()`](https://samleebyu.github.io/tseLCA/reference/posterior.md)
  : Posterior class-membership probabilities and modal class assignments

## Step 3: structural models

Relate the classes to covariates and distal outcomes.

- [`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md)
  : Relate latent classes to covariates (Step 3)
- [`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)
  : Relate latent classes to a distal outcome (Step 3)
- [`tse_twostep()`](https://samleebyu.github.io/tseLCA/reference/tse_twostep.md)
  : Two-step estimates of covariate effects
- [`predict(`*`<tseLCA_covariate>`*`)`](https://samleebyu.github.io/tseLCA/reference/predict.tseLCA_covariate.md)
  : Class-membership probabilities from a covariate model
- [`relevel(`*`<tseLCA_covariate>`*`)`](https://samleebyu.github.io/tseLCA/reference/relevel.tseLCA_covariate.md)
  : Change the reference class of a covariate model
- [`anova(`*`<tseLCA_covariate>`*`)`](https://samleebyu.github.io/tseLCA/reference/anova.tseLCA_covariate.md)
  : Wald tests of covariate terms
- [`omnibus_test()`](https://samleebyu.github.io/tseLCA/reference/omnibus_test.md)
  : Omnibus Wald test of class equality for a distal outcome

## Methods for fitted models

- [`summary(`*`<tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`coef(`*`<summary.tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`print(`*`<summary.tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`print(`*`<tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`summary(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`print(`*`<summary.tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  [`print(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/summary.tseLCA_structural.md)
  : Summarize a fitted tseLCA model
- [`coef(`*`<tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/coef.tseLCA_structural.md)
  [`coef(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/coef.tseLCA_structural.md)
  : Coefficients of a fitted tseLCA model
- [`vcov(`*`<tseLCA_structural>`*`)`](https://samleebyu.github.io/tseLCA/reference/vcov.tseLCA_structural.md)
  [`vcov(`*`<tseLCA_measurement>`*`)`](https://samleebyu.github.io/tseLCA/reference/vcov.tseLCA_structural.md)
  : Variance-covariance matrix of a fitted tseLCA model
- [`logLik(`*`<tseLCA>`*`)`](https://samleebyu.github.io/tseLCA/reference/logLik.tseLCA.md)
  [`nobs(`*`<tseLCA>`*`)`](https://samleebyu.github.io/tseLCA/reference/logLik.tseLCA.md)
  : Log-likelihood, number of observations, and information criteria
- [`plot(`*`<tseLCA>`*`)`](https://samleebyu.github.io/tseLCA/reference/plot.tseLCA.md)
  : Plot item-response probability profiles for a tseLCA model

## Estimation settings

- [`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md)
  : Estimation settings for tseLCA models

## Simulation

Data from the design of Bakk and Kuha (2018), used in the simulation
study.

- [`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md)
  : Generate one dataset following the Bakk & Kuha (2018) simulation
  design
- [`generate_all_conditions()`](https://samleebyu.github.io/tseLCA/reference/generate_all_conditions.md)
  : Generate datasets for all 18 conditions in the simulation design
- [`bk2018_params`](https://samleebyu.github.io/tseLCA/reference/bk2018_params.md)
  : Default population parameters for the Bakk & Kuha (2018) simulation
- [`draw_Zo()`](https://samleebyu.github.io/tseLCA/reference/draw_Zo.md)
  : Draw a continuous distal outcome given true class memberships
  (scenario "distal")
- [`draw_Zp()`](https://samleebyu.github.io/tseLCA/reference/draw_Zp.md)
  : Draw the covariate Zp ~ Uniform{1, 2, 3, 4, 5}
- [`draw_classes()`](https://samleebyu.github.io/tseLCA/reference/draw_classes.md)
  : Draw latent class memberships from their marginal distribution
- [`draw_classes_given_Zp()`](https://samleebyu.github.io/tseLCA/reference/draw_classes_given_Zp.md)
  : Draw latent classes conditional on the covariate (scenario
  "covariate")
- [`draw_indicators()`](https://samleebyu.github.io/tseLCA/reference/draw_indicators.md)
  : Draw binary indicators given true class memberships
- [`make_rho()`](https://samleebyu.github.io/tseLCA/reference/make_rho.md)
  : Build the item-response probability matrix for the simulation
- [`mnl_probs()`](https://samleebyu.github.io/tseLCA/reference/mnl_probs.md)
  : Compute multinomial logistic class probabilities given covariates

## Deprecated

The tseLCA 1.x interface. These functions keep working but will be
removed in a future version.

- [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
  : Three-step LCA estimation with covariates and/or distal outcomes
- [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  : Fit the LCA measurement model (Step 1)
- [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md)
  : Estimate covariate effects with measurement parameters fixed
  (two-step EM)
- [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md)
  : Estimate two-step covariate model with multilevLCA (optional
  reference path)
