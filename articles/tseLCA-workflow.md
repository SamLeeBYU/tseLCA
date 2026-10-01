# tseLCA Workflow

``` r

library(tseLCA)
```

## Overview

Latent class analysis (LCA) groups observations into unobserved classes
from a set of categorical indicators. Researchers usually also want to
know how the classes relate to other variables: *covariates* that
predict class membership, and *distal outcomes* that the classes
predict.

`tseLCA` does this with **three-step estimation**:

1.  **Measurement model.** Estimate the latent classes from the
    indicators alone
    ([`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)).
2.  **Classification.** Assign observations to classes and quantify the
    classification error
    ([`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)).
3.  **Structural model.** Relate the classes to covariates
    ([`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md))
    or distal outcomes
    ([`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)),
    correcting for the classification error of Step 2.

Because the measurement model is fixed before any structural variable
enters, the covariates and outcomes cannot change what the classes mean.
That is the main reason to prefer three-step over one-step estimation,
in which indicators, covariates, and outcomes are modeled jointly and
the class solution can shift with every change in the structural
specification. Assigning observations to classes and then analyzing the
assignments, as if they were the true classes, biases the structural
estimates toward zero; the bias-adjusted estimators of Bolck, Croon, and
Hagenaars (2004; “BCH”) and Vermunt (2010; “ML”) remove that bias, and
`tseLCA` adds standard errors that account for the uncertainty of the
Step-1 measurement model (Bakk, Oberski, and Vermunt 2014).

Each step returns an object that can be inspected with the usual R tools
([`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)) before moving
to the next. The one-call interface
[`tseLCA()`](https://samleebyu.github.io/tseLCA/reference/tseLCA.md)
runs all steps at once.

## Example data

[`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md)
simulates data from the design of Bakk and Kuha (2018): three classes,
six binary indicators `Y1`–`Y6`, and either a covariate `Zp`
(`scenario = "covariate"`) or a continuous distal outcome `Zo`
(`scenario = "distal"`). `separation` controls how well the indicators
separate the classes; the true class `X` is included for reference.

``` r

d <- generate_data(n = 1000, separation = "high", scenario = "covariate", seed = 1)
d$Zo <- draw_Zo(d$X, bk2018_params$distal_params) # add a distal outcome
head(d)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp           Zo
#> 1  1  1  1  0  0  0 2  1 -0.230182622
#> 2  0  1  1  0  0  0 3  4  0.769758959
#> 3  1  1  0  0  0  0 2  1  0.175955765
#> 4  0  1  1  0  1  0 2  2  2.670139514
#> 5  0  0  1  0  0  0 3  5  0.005199824
#> 6  1  0  1  1  1  1 1  3 -1.356208871
```

## Step 1: the measurement model

### Choosing the number of classes

The number of classes is chosen from the measurement model alone, before
any covariates or outcomes are considered. With several values of
`nclass`,
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
returns a class-enumeration table of fit statistics (see Nylund,
Asparouhov, and Muthén 2007, and Masyn 2013, for guidance on using
them).

``` r

f_items <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
sel <- tse_lca(f_items, data = d, nclass = 1:4)
sel
#> Latent class enumeration (measurement model)
#> 
#>  nclass logLik npar      AIC      BIC    SABIC entropy.R2 min.class
#>       1  -3989    6 7989.70  8019.15  8000.09          NA   1.00000
#>       2  -3232   13 6490.52  6554.33  6513.04      0.8651   0.38942
#>       3  -2960   20 5959.92  6058.07* 5994.55*     0.8720   0.31573
#>       4  -2953   27 5959.59* 6092.10  6006.35      0.8672   0.01781
#> 
#> * smallest value of each criterion. N = 1000
```

``` r

plot(sel)
```

![](tseLCA-workflow_files/figure-html/enumeration-plot-1.png)

The BIC and the sample-size adjusted BIC (SABIC) favor three classes.
The AIC, which penalizes additional parameters less, marginally prefers
four, but the fourth class holds less than 2% of the sample, a common
sign of over-extraction; the three-class model is also the more
interpretable. In practice the choice combines these criteria with class
sizes, separation (entropy), and substantive interpretation.
[`best_model()`](https://samleebyu.github.io/tseLCA/reference/best_model.md)
extracts the model selected by a criterion; `sel[[k]]` extracts the
`k`-class model.

``` r

m <- best_model(sel, criterion = "BIC")
```

### Inspecting the measurement model

``` r

summary(m)
#> Latent class measurement model
#>   Classes: 3   N: 1000
#>   Log-lik: -2959.9591 (df = 20)   AIC: 5959.92   BIC: 6058.07
#>   Entropy R²: 0.8720
#> 
#> Class sizes:
#>     C1     C2     C3 
#> 0.3501 0.3157 0.3341 
#> 
#> Item-response probabilities:
#>             C1     C2      C3
#> P(Y1|C) 0.8938 0.8915 0.07963
#> P(Y2|C) 0.9045 0.8671 0.13936
#> P(Y3|C) 0.8907 0.8505 0.07364
#> P(Y4|C) 0.8969 0.1062 0.10002
#> P(Y5|C) 0.9218 0.1032 0.10981
#> P(Y6|C) 0.9034 0.1081 0.09452
```

``` r

plot(m)
```

![](tseLCA-workflow_files/figure-html/measurement-plot-1.png)

[`class_sizes()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
and
[`item_probs()`](https://samleebyu.github.io/tseLCA/reference/class_sizes.md)
give the parameters on the probability scale;
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) give them on the
unconstrained log-ratio scale in which they are estimated.

``` r

class_sizes(m)
#>        C1        C2        C3 
#> 0.3501193 0.3157328 0.3341479
item_probs(m)
#>                C1        C2         C3
#> P(Y1|C) 0.8937779 0.8914562 0.07963285
#> P(Y2|C) 0.9044773 0.8670541 0.13936231
#> P(Y3|C) 0.8906931 0.8505388 0.07364458
#> P(Y4|C) 0.8969491 0.1062298 0.10001735
#> P(Y5|C) 0.9217944 0.1031584 0.10980622
#> P(Y6|C) 0.9034094 0.1080517 0.09451956
```

Latent class models can have several local optima. By default
[`tse_lca()`](https://samleebyu.github.io/tseLCA/reference/tse_lca.md)
uses the k-means initialization of **multilevLCA**, with additional
random starts when the entropy is low;
`control = tse_control(n_init = 20)` fits the model from 20 random
starts, and `start` fixes a starting classification.

## Step 2: classification

[`tse_classify()`](https://samleebyu.github.io/tseLCA/reference/tse_classify.md)
assigns each observation to its most likely class (modal assignment) or,
with `assignment = "proportional"`, to every class with its posterior
probability as weight. It reports the classification-error probabilities
$`P(W = s \mid X = t)`$ between the true class $`X`$ and the assigned
class $`W`$, which the Step-3 estimators correct for.

``` r

cl <- tse_classify(m)
cl
#> Latent class assignment (Step 2)
#>   Classes: 3   Assignment: modal   N: 1000   Entropy R²: 0.8721
#> 
#> Classification error probabilities P(W = s | X = t)
#> (rows: true class X; columns: assigned class W)
#>       W=C1  W=C2  W=C3
#> X=C1 0.973 0.023 0.004
#> X=C2 0.031 0.924 0.046
#> X=C3 0.004 0.025 0.971
#> 
#> Class proportions: estimated (measurement model) and assigned
#>              C1    C2    C3
#> estimated 0.350 0.316 0.334
#> assigned  0.352 0.308 0.340
```

With these well-separated classes, 92–97% of the members of each class
are assigned to it.
[`posterior()`](https://samleebyu.github.io/tseLCA/reference/posterior.md)
and
[`classes()`](https://samleebyu.github.io/tseLCA/reference/posterior.md)
return the posterior probabilities and the modal classes.

## Step 3: covariates

[`tse_covariate()`](https://samleebyu.github.io/tseLCA/reference/tse_covariate.md)
estimates a multinomial logistic regression of class membership on the
covariates. The default estimator is ML with standard errors corrected
for the Step-1 uncertainty.

``` r

fc <- tse_covariate(cl, ~ Zp)
summary(fc)
#> Three-step latent class model: covariates
#>   Classes: 3   Estimator: ML   N: 1000
#>   Log-lik: -2717.7982 (df = 22)   AIC: 5479.60   BIC: 5587.57
#>   Entropy R² (covariate-adjusted): 0.8861
#> 
#> Covariate effects on class membership (multinomial logit):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):C2  1.96000    0.24865   7.882 3.21e-15 ***
#> Zp:C2          -0.90913    0.11546  -7.874 3.43e-15 ***
#> (Intercept):C3 -3.79821    0.37232 -10.201  < 2e-16 ***
#> Zp:C3           1.03333    0.09398  10.996  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The coefficients are named `covariate:class` and the first class is the
reference. The generic tools work as usual:

``` r

confint(fc)
#>                     2.5 %     97.5 %
#> (Intercept):C2  1.4726472  2.4473489
#> Zp:C2          -1.1354206 -0.6828305
#> (Intercept):C3 -4.5279423 -3.0684755
#> Zp:C3           0.8491445  1.2175239
AIC(fc)
#> [1] 5479.596
anova(fc) # Wald test of each covariate term, across all classes
#> Wald tests of covariate terms (all class contrasts)
#> 
#>    Df Chisq Pr(>Chisq)    
#> Zp  2 183.4  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`predict()`](https://rdrr.io/r/stats/predict.html) gives the
class-membership probabilities at given covariate values:

``` r

predict(fc, newdata = data.frame(Zp = 1:5))
#>          C1         C2         C3
#> 1 0.2548985 0.72904693 0.01605453
#> 2 0.4293142 0.49469228 0.07599355
#> 3 0.5097606 0.23664534 0.25359408
#> 4 0.3868254 0.07234665 0.54082798
#> 5 0.1998141 0.01505572 0.78513019
```

### Estimators and standard errors

The BCH estimator and the uncorrected three-step estimator (which
analyzes the assigned classes as if they were the true classes) are
available for comparison; so is the two-step estimator of Bakk and Kuha
(2018).

``` r

fc_bch <- tse_covariate(cl, ~ Zp, method = "BCH")
fc_raw <- tse_covariate(cl, ~ Zp, method = "none")
ft <- tse_twostep(m, ~ Zp)
round(cbind(
  ML = coef(fc), BCH = coef(fc_bch), uncorrected = coef(fc_raw), two.step = coef(ft)
), 3)
#>                    ML    BCH uncorrected two.step
#> (Intercept):C2  1.960  1.911       1.521    1.856
#> Zp:C2          -0.909 -0.873      -0.703   -0.866
#> (Intercept):C3 -3.798 -3.901      -3.156   -3.729
#> Zp:C3           1.033  1.055       0.875    1.020
```

The uncorrected estimates are attenuated toward zero (the true slopes
are $`-1`$ and $`1`$); the bias-adjusted estimates are not.

`se = "robust"` omits the Step-1 correction. With well-separated classes
the two are close; the correction matters more when classes are less
distinct.

``` r

round(cbind(
  corrected = sqrt(diag(vcov(fc))),
  robust = sqrt(diag(vcov(tse_covariate(cl, ~ Zp, se = "robust"))))
), 4)
#>                corrected robust
#> (Intercept):C2    0.2487 0.2417
#> Zp:C2             0.1155 0.1093
#> (Intercept):C3    0.3723 0.3608
#> Zp:C3             0.0940 0.0923
```

Under low separation, proportional assignment
(`tse_classify(m, assignment = "proportional")`) with the ML estimator
is generally the most reliable choice.

### Reference class and covariate formulas

`ref` (or [`relevel()`](https://rdrr.io/r/stats/relevel.html) on a
fitted model) changes the reference class. Covariates follow the usual
formula syntax, including factors, interactions, and transformations. A
variable created after classification is supplied with `data`, which
must hold the classified rows (it may add columns).

``` r

coef(relevel(fc, ref = "C3"), matrix = TRUE)
#>                    C1        C2
#> (Intercept)  3.798209  5.758207
#> Zp          -1.033334 -1.942460

d$group <- factor(ifelse(d$Zp > 3, "high", "low"))
anova(tse_covariate(cl, ~ Zp + group, data = d))
#> Wald tests of covariate terms (all class contrasts)
#> 
#>       Df   Chisq Pr(>Chisq)    
#> Zp     2 98.3038     <2e-16 ***
#> group  2  1.5912     0.4513    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Step 3: distal outcomes

[`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)
estimates the distribution of a distal outcome in each class. Outcomes
can be `"gaussian"` (class means with a common variance), `"poisson"`,
`"binomial"`, or `"multinomial"` (nominal).

``` r

fd <- tse_distal(cl, Zo ~ 1)
summary(fd)
#> Three-step latent class model: distal outcome
#>   Classes: 3   Estimator: ML   Family: gaussian   N: 1000
#>   Log-lik: -4429.4154 (df = 24)   AIC: 8906.83   BIC: 9024.62
#> 
#> Distal outcome means by class:
#>       Estimate Std. Error z value Pr(>|z|)    
#> mu_C1 -0.88911    0.05593 -15.896   <2e-16 ***
#> mu_C2  1.04803    0.06688  15.669   <2e-16 ***
#> mu_C3 -0.10743    0.06087  -1.765   0.0776 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`omnibus_test()`](https://samleebyu.github.io/tseLCA/reference/omnibus_test.md)
tests whether the outcome’s distribution differs across classes:

``` r

omnibus_test(fd)
#> 
#>  Wald test of equal distal outcome distributions across latent classes
#> 
#> data:  gaussian distal outcome, 3 classes
#> W = 510.3, df = 2, p-value < 2.2e-16
```

For a nominal outcome, `coef(fit, matrix = TRUE)` gives the
class-by-category probability matrix:

``` r

d$Zcat <- cut(d$Zo, c(-Inf, -0.5, 0.5, Inf), labels = c("low", "mid", "high"))
fm <- tse_distal(cl, Zcat ~ 1, family = "multinomial", data = d)
round(coef(fm, matrix = TRUE), 3)
#>      low   mid  high
#> C1 0.627 0.312 0.061
#> C2 0.086 0.218 0.696
#> C3 0.361 0.350 0.289
omnibus_test(fm)
#> 
#>  Wald test of equal distal outcome distributions across latent classes
#> 
#> data:  multinomial distal outcome, 3 classes
#> W = 417.4, df = 4, p-value < 2.2e-16
```

### Covariates and a distal outcome

Passing a covariate model to
[`tse_distal()`](https://samleebyu.github.io/tseLCA/reference/tse_distal.md)
fits both parts. The class prior then depends on the covariates, and the
uncertainty of the covariate model is propagated to the distal
estimates.

``` r

fb <- tse_distal(fc, Zo ~ 1)
fb
#> Three-step latent class model: covariates and distal outcome
#>   Classes: 3   Estimator: ML   Family: gaussian   N: 1000
#>   Log-lik: -4176.2103 (df = 26)   AIC: 8404.42   BIC: 8532.02
#> 
#> Covariate effects on class membership (multinomial logit):
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):C2  1.96000    0.24865   7.882 3.21e-15 ***
#> Zp:C2          -0.90913    0.11546  -7.874 3.43e-15 ***
#> (Intercept):C3 -3.79821    0.37232 -10.201  < 2e-16 ***
#> Zp:C3           1.03333    0.09398  10.996  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Distal outcome means by class:
#>       Estimate Std. Error z value Pr(>|z|)    
#> mu_C1 -0.89181    0.05550 -16.070   <2e-16 ***
#> mu_C2  1.05856    0.06580  16.087   <2e-16 ***
#> mu_C3 -0.07797    0.05580  -1.397    0.162    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## All steps in one call

[`tseLCA()`](https://samleebyu.github.io/tseLCA/reference/tseLCA.md)
runs the three steps from one formula,
`indicators ~ covariates | distal outcome`. The components remain
available through
[`measurement()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`classification()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
[`covariate()`](https://samleebyu.github.io/tseLCA/reference/measurement.md),
and
[`distal()`](https://samleebyu.github.io/tseLCA/reference/measurement.md).

``` r

fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
all.equal(coef(fit), coef(fb))
#> [1] TRUE
round(classification(fit)$D, 3)
#>       W=C1  W=C2  W=C3
#> X=C1 0.973 0.023 0.004
#> X=C2 0.031 0.924 0.046
#> X=C3 0.004 0.025 0.971
```

## A measurement model from another sample

Because the measurement model is fixed in Step 1, it can be estimated on
one sample and applied to another, for example when covariates are
observed only in a subsample. `tse_classify(m, newdata = ...)`
classifies the new sample with the existing measurement model; the
Step-1 uncertainty in later steps is that of the sample the model was
estimated on.

``` r

sub <- d[1:300, ]
fc_sub <- tse_covariate(tse_classify(m, newdata = sub), ~ Zp)
coef(fc_sub)
#> (Intercept):C2          Zp:C2 (Intercept):C3          Zp:C3 
#>      2.1769002     -1.0400321     -3.4247725      0.9448009
nobs(fc_sub)
#> [1] 300
```

## Missing data and indicator coding

Indicators can be factors, logicals, character variables, or numeric
codes in any coding; their categories are stored with the model and
reused on new data. With `missing = "fiml"`, observations with some
missing indicators are kept (full-information maximum likelihood); the
default drops them. Rows with a missing covariate or distal outcome are
dropped from that Step-3 model only.

``` r

d_miss <- d
set.seed(2)
d_miss$Y1[sample(nrow(d), 100)] <- NA
d_miss$Y2 <- factor(d_miss$Y2, labels = c("no", "yes"))
m_fiml <- tse_lca(f_items, data = d_miss, nclass = 3, missing = "fiml")
nobs(m_fiml)
#> [1] 1000
nobs(tse_lca(f_items, data = d_miss, nclass = 3))
#> [1] 900
```

## Estimation settings

[`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md)
collects the numerical settings: iteration limits and tolerances for
Step 1 and Step 3, random starts, the boundary tolerance for the Step-1
variance, and the information matrix used for Step-3 standard errors.

``` r

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

## Migrating from tseLCA 1.x

[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
still works (with the same estimates) but is deprecated. The table maps
its arguments to the new interface.

| [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md) | tseLCA 2.0 |
|----|----|
| `Y.names`, `n_classes` | `tse_lca(cbind(...) ~ 1, nclass = )` |
| `Zp.names` | `tse_covariate(, ~ ...)` |
| `Zo.name`, `family` | `tse_distal(, outcome ~ 1, family = )` |
| `step1` (measurement model from another sample) | `tse_classify(, newdata = )` |
| `startval`, `n_init` | `tse_lca(start = )`, `tse_control(n_init = )` |
| `use.modal.assignment` | `tse_classify(assignment = )` |
| `use.bch`, `use.simple.cov` | `method = "BCH"`, `se = "robust"` |
| `rebase` | `ref`, or [`relevel()`](https://rdrr.io/r/stats/relevel.html) |
| `incomplete` | `tse_lca(missing = "fiml")` |
| other tuning arguments | [`tse_control()`](https://samleebyu.github.io/tseLCA/reference/tse_control.md) |
| `get.twostep.vcov` | `tse_twostep(se = TRUE)` |

[`coef()`](https://rdrr.io/r/stats/coef.html) now returns a named vector
matching [`vcov()`](https://rdrr.io/r/stats/vcov.html) (use
`coef(fit, matrix = TRUE)` for the coefficient matrix), so
[`confint()`](https://rdrr.io/r/stats/confint.html) works. See `NEWS.md`
for all changes, including bug fixes affecting results for polytomous
indicators, Gaussian distal outcomes, and proportional-assignment ML
distal models.

## References

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, 22(4), 520–540.

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311.

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, 12(1), 3–27.

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*,
Vol. 2, 551–611. Oxford University Press.

Nylund, K. L., Asparouhov, T., & Muthén, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, 14(4),
535–569.

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
