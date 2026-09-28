# tests/testthat/test-bugfixes-2.0.R
#
# Regression tests for bugs in tseLCA 1.1.1 fixed in 2.0:
#   1. Measurement-only posteriors/classifications were not in data-row order
#      (taken from multilevLCA's fit0$mU, which is sorted by response pattern).
#   2. Polytomous items were decoded from fit0$mU as category 0 on the
#      listwise path, making the Step-1 information matrix singular: the
#      measurement vcov and the corrected Step-3 vcov were all NA.
#   3. AIC/BIC parameter counts counted one-hot indicator columns (covariate
#      models) and ignored covariate coefficients (combined models).

dl <- v1_data()

fit_pair <- function(d, ...) {
  set.seed(1L)
  m <- three_step(d, v1_items, 3L, ...)
  set.seed(1L)
  f <- suppressWarnings(three_step(
    d, v1_items, 3L, Zp.names = "Zp", step1 = m, use.simple.cov = TRUE, ...
  ))
  list(m = m, f = f)
}

test_that("measurement-only posteriors are in data-row order", {
  cases <- list(
    binary = list(d = dl$cov_high),
    polytomous = list(d = dl$poly),
    fiml = list(d = dl$miss, incomplete = TRUE)
  )
  for (nm in names(cases)) {
    args <- cases[[nm]]
    p <- do.call(fit_pair, c(list(args$d), args[-1]))
    # Step-2 posteriors (computed independently in lca_step2) for the same rows
    expect_equal(posterior(p$m), posterior(p$f), tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
    expect_equal(classes(p$m), classes(p$f), info = nm)
  }
})

test_that("measurement-only posteriors follow a rebased reference class", {
  set.seed(1L)
  m1 <- three_step(dl$cov_high, v1_items, 3L)
  set.seed(1L)
  m2 <- three_step(dl$cov_high, v1_items, 3L, rebase = "C2")
  expect_equal(posterior(m2), posterior(m1)[, c(2L, 1L, 3L)], ignore_attr = TRUE)
})

test_that("modal classes recover the true classes under high separation", {
  set.seed(1L)
  m <- three_step(dl$cov_high, v1_items, 3L)
  tab <- table(classes(m), dl$cov_high$X)
  # agreement up to label switching
  expect_gt(sum(apply(tab, 1L, max)) / nrow(dl$cov_high), 0.85)
})

test_that("polytomous measurement vcov and corrected Step-3 vcov are finite", {
  for (nm in c("poly", "sparse")) {
    set.seed(1L)
    m <- three_step(dl[[nm]], v1_items, 3L)
    expect_false(anyNA(vcov(m)), info = nm)
    set.seed(1L)
    f <- three_step(dl[[nm]], v1_items, 3L, Zp.names = "Zp", step1 = m)
    set.seed(1L)
    fs <- three_step(dl[[nm]], v1_items, 3L, Zp.names = "Zp", step1 = m,
                     use.simple.cov = TRUE)
    V <- vcov(f)
    expect_false(anyNA(V), info = nm)
    # correcting for Step-1 uncertainty cannot shrink the variances
    expect_true(all(diag(V) >= diag(vcov(fs)) - 1e-10), info = nm)
  }
})

test_that("the sparse data set has item probabilities at the boundary", {
  set.seed(1L)
  m <- three_step(dl$sparse, v1_items, 3L)
  expect_lt(min(item_probs(m)), 1e-3)
})

test_that("fit0$mU fallback decodes both multilevLCA polytomous codings", {
  d <- dl$poly
  d_miss <- d
  set.seed(9L)
  d_miss$Y1[sample(nrow(d), 40L)] <- NA
  for (case in list(list(d = d, inc = FALSE), list(d = d_miss, inc = TRUE))) {
    set.seed(1L)
    m <- suppressWarnings(three_step(case$d, v1_items, 3L, incomplete = case$inc))
    s1 <- m$measurement_model
    decoded <- step1_sample(list(fit0 = s1$fit0), s1$ivItemcat)
    # mU is sorted by response pattern: compare the multisets of rows
    key <- function(M) sort(apply(M, 1L, paste, collapse = ""))
    Yd <- decoded$Y.exp
    Ys <- s1$Y.exp
    if (!is.null(decoded$mDesign)) Yd[decoded$mDesign == 0] <- NA
    if (!is.null(s1$mDesign.exp)) Ys[s1$mDesign.exp == 0] <- NA
    expect_equal(key(Yd), key(Ys), info = paste("incomplete =", case$inc))
  }
})

test_that("a Step-1 model reused through `step1` keeps its own sample", {
  set.seed(1L)
  m <- three_step(dl$cov_high, v1_items, 3L)
  expect_equal(nrow(m$measurement_model$Y.exp), nrow(dl$cov_high))
  small <- dl$cov_mid[1:200, ]
  set.seed(1L)
  f <- three_step(small, v1_items, 3L, Zp.names = "Zp", step1 = m)
  expect_equal(nrow(f$measurement_model$Y.exp), nrow(dl$cov_high))
  expect_equal(nrow(posterior(f)), 200L)
  expect_false(anyNA(vcov(f)))
})

test_that("parameter counts behind logLik/AIC/BIC", {
  n_items <- 6L
  cfg <- function(nm) v1_fit(v1_configs[[nm]], dl)
  # measurement: (T-1) class sizes + T * sum(K-1) item parameters
  expect_equal(attr(logLik(cfg("meas_high")), "df"), 2 + 3 * n_items)
  expect_equal(attr(logLik(cfg("meas_poly")), "df"), 2 + 3 * 2 * n_items)
  # covariate: item parameters + Q * (T-1) logit coefficients (Q = 2)
  expect_equal(attr(logLik(cfg("cov_ml_modal")), "df"), 3 * n_items + 2 * 2)
  expect_equal(attr(logLik(cfg("cov_poly")), "df"), 3 * 2 * n_items + 2 * 2)
  # distal: class sizes + item parameters + T class means (+ sigma2 for
  # the gaussian family)
  expect_equal(attr(logLik(cfg("dis_gauss_ml")), "df"), 2 + 3 * n_items + 3 + 1)
  expect_equal(attr(logLik(cfg("dis_poisson")), "df"), 2 + 3 * n_items + 3)
  expect_equal(attr(logLik(cfg("dis_multinomial")), "df"), 2 + 3 * n_items + 3 * 2)
  # combined: covariate coefficients replace class sizes
  expect_equal(attr(logLik(cfg("both_ml_prop")), "df"), 2 * 2 + 3 * n_items + 3 + 1)
  fit <- cfg("cov_ml_modal")
  expect_equal(AIC(fit), -2 * fit$llik + 2 * 22)
  expect_equal(BIC(fit), -2 * fit$llik + 22 * log(nobs(fit)))
})

test_that("covariate models without an intercept can be fitted", {
  # 1.1.1 labeled the coefficient rows c("Intercept", Zp.names) regardless of
  # include.intercept and failed with a dimnames error.
  d <- generate_data(300L, "high", "covariate", seed = 1L)
  set.seed(1L)
  f <- three_step(d, v1_items, 3L, Zp.names = "Zp",
                  include.intercept = FALSE, use.simple.cov = TRUE)
  expect_equal(dimnames(coef(f, matrix = TRUE)), list("Zp", c("C2", "C3")))
  expect_equal(names(coef(f)), rownames(vcov(f)))
})
