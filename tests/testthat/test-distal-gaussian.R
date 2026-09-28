# tests/testthat/test-distal-gaussian.R
#
# Three-step ML for a gaussian distal outcome estimates the common
# within-class variance sigma2 jointly with the class means (tseLCA 2.0; it
# was fixed at 1 in 1.x, which biased the means whenever the outcome's
# variance differed from 1). These tests check the analytic score, Hessian,
# and Step-1/Step-2 uncertainty propagation against numerical derivatives,
# and the estimator's scale equivariance.

# Synthetic Step-3 inputs, independent of any Step-1 fit.
gauss_inputs <- function(N = 400L, iT = 3L, seed = 1L) {
  set.seed(seed)
  X <- sample.int(iT, N, replace = TRUE, prob = c(.3, .3, .4))
  pwx <- matrix(c(.8, .1, .1, .15, .7, .15, .1, .2, .7), iT) # columns sum to 1
  W <- vapply(X, function(x) sample.int(iT, 1L, prob = pwx[, x]), integer(1))
  w.is <- diag(iT)[W, ]
  z <- rnorm(N, c(-1, 1, 0)[X], 1.5)
  Zc <- cbind(1, rnorm(N))
  list(N = N, iT = iT, pwx = pwx, w.is = w.is, z = z, Zc = Zc,
       gamma = c(0.2, -0.3, 0.4, 0.5), theta = c(-0.8, 0.9, 0.1, 2.1))
}

p_xz <- function(gamma, Zc) {
  eta <- cbind(0, Zc %*% matrix(gamma, ncol = 2L))
  e <- exp(eta - apply(eta, 1L, max))
  e / rowSums(e)
}

# column-softmax parameterization of P(W = s | X = t): off-diagonal
# log-ratios log(pwx[s, t] / pwx[t, t]), ordered t0 then s0 (as in
# lca_vcov_distal)
theta2_of <- function(pwx) {
  iT <- ncol(pwx)
  unlist(lapply(seq_len(iT), function(t) log(pwx[-t, t] / pwx[t, t])))
}
pwx_of <- function(theta2, iT) {
  M <- matrix(0, iT, iT)
  k <- 0L
  for (t in seq_len(iT)) {
    M[-t, t] <- theta2[k + seq_len(iT - 1L)]
    k <- k + iT - 1L
  }
  E <- exp(M)
  sweep(E, 2L, colSums(E), "/")
}

num_jac <- function(f, x, eps = 1e-6) {
  f0 <- f(x)
  vapply(seq_along(x), function(j) {
    h <- eps * max(abs(x[j]), 1)
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(length(f0)))
}

nll_gauss <- function(theta, pi_mat, pwx, z, w.is, iT) {
  mu <- theta[1:iT]
  s2 <- theta[iT + 1L]
  lf <- -0.5 * outer(z, mu, "-")^2 / s2 - 0.5 * log(2 * pi * s2)
  -sum(log(rowSums(pi_mat * (w.is %*% pwx) * exp(lf))))
}

test_that("analytic score and Hessian match numerical derivatives", {
  g <- gauss_inputs()
  pi_mat <- p_xz(g$gamma, g$Zc)
  f <- function(th) nll_gauss(th, pi_mat, g$pwx, g$z, g$w.is, g$iT)
  score <- colSums(ml_score_distal_gaussian(g$theta, pi_mat, g$pwx, g$z, g$w.is, g$iT))
  expect_equal(score, -as.vector(num_jac(f, g$theta)), tolerance = 1e-6)
  H <- ml_hessian_distal_gaussian(g$theta, pi_mat, g$pwx, g$z, g$w.is, g$iT)
  H_num <- num_jac(function(th) {
    colSums(-ml_score_distal_gaussian(th, pi_mat, g$pwx, g$z, g$w.is, g$iT))
  }, g$theta)
  expect_equal(H, H_num, tolerance = 1e-6)
})

test_that("Step-1 and Step-2 uncertainty terms match numerical propagation", {
  g <- gauss_inputs()
  iT <- g$iT
  pi_mat <- p_xz(g$gamma, g$Zc)
  theta <- g$theta
  Psi <- function(th, pwx, gamma) {
    colSums(ml_score_distal_gaussian(th, p_xz(gamma, g$Zc), pwx, g$z, g$w.is, iT))
  }
  H.3.inv <- solve(ml_hessian_distal_gaussian(theta, pi_mat, g$pwx, g$z, g$w.is, iT))
  J <- ml_score_distal_gaussian(theta, pi_mat, g$pwx, g$z, g$w.is, iT)
  meat <- crossprod(J)

  # Arbitrary positive-definite input variances; J.2 = identity maps the
  # Step-1 variance straight onto the classification-error parameters.
  set.seed(2L)
  A <- matrix(rnorm(36), 6)
  Sigma.1 <- crossprod(A) / 500
  B <- matrix(rnorm(16), 4)
  Sigma.3 <- crossprod(B) / 500

  C1_num <- num_jac(function(t2) Psi(theta, pwx_of(t2, iT), g$gamma), theta2_of(g$pwx))
  C_num <- num_jac(function(gm) Psi(theta, g$pwx, gm), g$gamma)
  # lca_vcov_distal treats the Step-3 score as d(-neg.ll); C terms enter the
  # sandwich only through C Sigma C', so their sign convention is irrelevant.
  V_num <- H.3.inv %*%
    (meat + C1_num %*% Sigma.1 %*% t(C1_num) + C_num %*% Sigma.3 %*% t(C_num)) %*%
    H.3.inv

  p.zx <- function(params) {
    s2 <- params[iT + 1L]
    -0.5 * outer(g$z, params[1:iT], "-")^2 / s2 - 0.5 * log(2 * pi * s2)
  }
  V <- lca_vcov_distal(
    mu_hat = theta,
    three_step.score = function(th) {
      ml_score_distal_gaussian(th, pi_mat, g$pwx, g$z, g$w.is, iT)
    },
    pi_adj = pi_mat, w.is = g$w.is, p.wx_mat = g$pwx, p.zx = p.zx,
    family = "gaussian", H.3.inv = H.3.inv, Sigma.1 = Sigma.1,
    s2 = list(J.2 = diag(6)), Sigma.3 = Sigma.3, s3.par = g$gamma,
    p.xz.cov = function(params) p_xz(as.vector(params), g$Zc), Z_mat_cov = g$Zc,
    iT = iT, use.simple.cov = FALSE, use.bch = FALSE,
    unit_scores = function(th) unit_scores_distal_gaussian(th, g$z, iT)
  )
  expect_equal(unname(V), unname(V_num), tolerance = 1e-5)
  expect_equal(rownames(V), c("mu_C1", "mu_C2", "mu_C3", "sigma2"))
})

test_that("gaussian distal estimates are scale equivariant (ML and BCH)", {
  d <- generate_data(600L, "mid", "distal", seed = 7L)
  d$Zo10 <- 10 * d$Zo
  set.seed(1L)
  m <- three_step(d, paste0("Y", 1:6), 3L)
  for (bch in c(FALSE, TRUE)) {
    for (modal in c(TRUE, FALSE)) {
      f1 <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo", step1 = m,
                       use.bch = bch, use.modal.assignment = modal)
      f10 <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo10", step1 = m,
                        use.bch = bch, use.modal.assignment = modal)
      info <- paste("bch =", bch, "modal =", modal)
      expect_equal(coef(f10), 10 * coef(f1), tolerance = 1e-4, info = info)
      expect_equal(sqrt(diag(vcov(f10))), 10 * sqrt(diag(vcov(f1))),
                   tolerance = 1e-4, info = info)
      expect_equal(f10$sigma2[["estimate"]], 100 * f1$sigma2[["estimate"]],
                   tolerance = 1e-4, info = info)
    }
  }
})

test_that("ML with modal assignment recovers means and variance of a wide outcome", {
  d <- generate_data(3000L, "high", "distal", seed = 11L)
  mu <- c(-1, 1, 0)[d$X]
  d$Zo3 <- mu + 3 * (d$Zo - mu) # same class means, residual SD 3
  set.seed(1L)
  m <- three_step(d, paste0("Y", 1:6), 3L)
  f <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo3", step1 = m)
  perm <- apply(table(classes(m), d$X), 2L, which.max)
  expect_equal(unname(coef(f)[perm]), c(-1, 1, 0), tolerance = 0.25)
  expect_equal(f$sigma2[["estimate"]], 9, tolerance = 0.1)
  expect_true(is.finite(f$sigma2[["se"]]))
})
