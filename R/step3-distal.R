# tseLCA/R/step3-distal.R
#
# Step 3 for distal outcomes: class-specific outcome parameters (BCH or ML).
# The ML likelihood machinery shared by all families is in R/distal-ml.R.

#' Step 3 (distal, multinomial): estimate class-conditional category
#' probabilities for a nominal categorical distal outcome
#'
#' Estimates the T x C matrix `pi_hat[t, c] = P(Zo = c | X = t)` for a
#' saturated (nominal) multinomial distal outcome with C categories, using the
#' closed-form weighted-proportion estimator
#' `pi_hat[t, c] = sum_i w_it * 1(y_i = c) / sum_i w_it`, using either the
#' BCH weight matrix (fixed given Step 1/2) or, for ML, EM-updated posterior
#' responsibilities (both give the same closed-form M-step; only the E-step
#' weights differ). Returns the same contract as `lca_step3.distal()`
#' (`res$par`, `H.3.inv`, `three_step.score`), with `res$par` the
#' column-major-flattened `pi_hat` (`matrix(par, nrow = iT, ncol = C)`
#' recovers it) so downstream sandwich-variance code is unchanged.
#'
#' For BCH, the score is the moment/estimating-equation form the weighted
#' proportion solves, `s_i,(t,c) = w_it * (1(y_i=c) - pi_hat[t,c])`; its
#' bread is exactly `diag(1 / colSums(w_it))` (repeated across categories)
#' since `w_it` does not depend on `pi_hat`. For ML, the same estimating
#' equation is used with `w_it` replaced by the (theta-dependent) E-step
#' weights `lambda_it = sum_s w_is R_ist` of the expanded-data likelihood
#' (R/distal-ml.R); its bread is `solve(-Jacobian(Psi))`, with the Jacobian
#' given in closed form by `distal_multinomial_jacobian()`.
#' @noRd
lca_step3.distal.multinomial <- function(
  Y_cat,
  C,
  iT,
  covariate.tol,
  use.bch = FALSE,
  w.is_cc = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  vPi = NULL,
  pi_mat = NULL,
  verbose = FALSE
) {
  N <- length(Y_cat)
  onehot_y <- matrix(0, N, C)
  onehot_y[cbind(seq_len(N), Y_cat)] <- 1

  pi_s <- if (!is.null(pi_mat)) {
    pi_mat
  } else {
    matrix(vPi, ncol = iT, nrow = N, byrow = TRUE)
  }

  # T x C closed-form weighted-proportion M-step for any N x T weight matrix.
  weighted_props <- function(w) {
    sweep(t(w) %*% onehot_y, 1, colSums(w), "/")
  }

  # N x (iT*C) estimating-equation matrix, column-major in (t, c): columns
  # 1:iT are category 1 for classes 1..iT, columns (iT+1):(2*iT) are
  # category 2, etc. -- matches matrix(theta, nrow = iT, ncol = C).
  score_matrix <- function(w, pi_hat) {
    J <- matrix(0, N, iT * C)
    for (c in seq_len(C)) {
      idx <- ((c - 1L) * iT + 1L):(c * iT)
      J[, idx] <- w * onehot_y[, c] - sweep(w, 2, pi_hat[, c], "*")
    }
    J
  }

  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx)
    w_colsums <- colSums(w.it)
    if (any(w_colsums <= 0)) {
      stop(
        "BCH weights have non-positive column sums for at least one class. ",
        "The variance-covariance matrix will not be positive semi-definite. ",
        "Consider use.bch = FALSE.",
        call. = FALSE
      )
    }
    pi_hat <- weighted_props(w.it)
    log_pzx_bch <- log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300))

    three_step.score <- function(params) {
      score_matrix(w.it, matrix(params, nrow = iT, ncol = C))
    }

    H.3.inv <- diag(rep(1 / w_colsums, times = C))
    res <- list(
      par = as.vector(pi_hat),
      value = -sum(w.it * log_pzx_bch),
      convergence = 0L
    )
  } else {
    # EM: E-step (posterior responsibilities) / M-step (weighted
    # proportions, closed form). Initialize from a smoothed crosstab of the
    # modal Step-2 assignment against the observed category.
    init_class <- factor(max.col(w.is_cc), levels = seq_len(iT))
    init_cat <- factor(Y_cat, levels = seq_len(C))
    pi_hat <- unclass(table(init_class, init_cat)) + 0.5
    pi_hat <- pi_hat / rowSums(pi_hat)
    storage.mode(pi_hat) <- "double"

    # Per-record posteriors (see R/distal-ml.R) at a T x C probability matrix
    records_at <- function(pi_hat) {
      distal_records(
        log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300)), # N x T
        pi_s,
        pwx
      )
    }

    ll_prev <- -Inf
    for (iter in seq_len(em.maxIter)) {
      rec <- records_at(pi_hat)
      pi_hat_new <- weighted_props(distal_lambda(rec$R, w.is_cc))
      pi_hat_new <- pmax(pmin(pi_hat_new, 1 - 1e-10), 1e-10)
      pi_hat_new <- pi_hat_new / rowSums(pi_hat_new)

      ll_new <- distal_loglik(rec$logM, w.is_cc)
      if (iter > 1L && abs(ll_new - ll_prev) < covariate.tol) {
        pi_hat <- pi_hat_new
        if (verbose) {
          message(sprintf(
            "Multinomial ML EM converged in %d iterations.",
            iter
          ))
        }
        break
      }
      pi_hat <- pi_hat_new
      ll_prev <- ll_new
      if (iter == em.maxIter) {
        warning("Multinomial ML EM reached maximum iterations.")
      }
    }

    three_step.score <- function(params) {
      pi_hat_p <- matrix(params, nrow = iT, ncol = C)
      score_matrix(distal_lambda(records_at(pi_hat_p)$R, w.is_cc), pi_hat_p)
    }

    theta_hat <- as.vector(pi_hat)
    rec <- records_at(pi_hat)
    Jac <- distal_multinomial_jacobian(pi_hat, rec, w.is_cc, Y_cat)

    H.3.inv <- tryCatch(
      qr.solve(-Jac),
      error = function(e) {
        warning("Hessian inversion failed. SEs will be NA.")
        matrix(NA_real_, iT * C, iT * C)
      }
    )
    res <- list(
      par = theta_hat,
      value = -distal_loglik(rec$logM, w.is_cc),
      convergence = 0L
    )
  }

  list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  )
}

#' Step 3 (distal): estimate class-specific distal outcome parameters
#'
#' Estimates mu = (mu_1, ..., mu_T) for Gaussian (means), Poisson (log-rates),
#' or Binomial (logits) distal outcomes with either BCH (closed-form or Newton-Rhapson) or
#' ML EM. Returns the parameter estimates, the inverted Hessian H.3.inv, and
#' the case-wise score function three_step.score for sandwich variance
#' propagation in lca_vcov_distal.
#' @noRd
lca_step3.distal <- function(
  neg.ll,
  beta_init,
  iT,
  covariate.tol,
  use.bch = FALSE,
  Zo_cc = NULL,
  w.is_cc = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  family = "gaussian",
  p.zx = NULL,
  vPi = NULL,
  pi_mat = NULL,
  verbose = FALSE
) {
  N <- length(Zo_cc)
  #use covariate-adjusted pi if provided, otherwise flat vPi
  pi_s <- if (!is.null(pi_mat)) {
    pi_mat
  } else {
    matrix(vPi, ncol = iT, nrow = N, byrow = TRUE)
  }

  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx) # N x T

    score_nt_bch <- function(mu) {
      if (family == "gaussian") {
        resid <- outer(Zo_cc, mu, "-")
        w.it * resid
      } else if (family == "poisson") {
        mu_val <- exp(mu)
        w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
      } else if (family == "binomial") {
        mu_val <- 1 / (1 + exp(-mu))
        w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
      }
    }

    w_colsums <- colSums(w.it)

    if (any(w_colsums < 0)) {
      stop(
        "BCH weights have negative column sums for at least one class. ",
        "The variance-covariance matrix will not be positive semi-definite. ",
        "Consider use.bch = FALSE."
      )
    }

    if (family == "gaussian") {
      beta <- colSums(w.it * Zo_cc) / w_colsums # closed-form weighted mean
      resid <- outer(Zo_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:iT]
        resid <- outer(Zo_cc, mu, "-")
        w.it * resid / sigma2
      }

      H.3.inv <- diag(sigma2 / w_colsums)

      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    } else {
      beta <- beta_init
      for (nr in seq_len(em.maxIter)) {
        grad_vec <- colSums(score_nt_bch(beta))

        if (family == "gaussian") {
          H_diag <- w_colsums
        } else if (family == "poisson") {
          H_diag <- exp(beta) * w_colsums
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-beta))
          H_diag <- mu_val * (1 - mu_val) * w_colsums
        }

        direction <- grad_vec / H_diag

        step <- 1.0
        ll_cur <- -neg.ll(beta)
        for (ls in seq_len(20L)) {
          beta_new <- beta + step * direction
          ll_new <- tryCatch(-neg.ll(beta_new), error = function(e) -Inf)
          if (is.finite(ll_new) && ll_new > ll_cur) {
            break
          }
          step <- step * 0.5
        }

        delta <- step * direction
        beta <- beta + delta

        if (max(abs(delta)) < covariate.tol) {
          if (verbose) {
            message(sprintf("BCH NR converged in %d iterations.", nr))
          }
          break
        }
        if (nr == em.maxIter) warning("BCH NR reached maximum iterations.")
      }

      resid <- outer(Zo_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:iT]
        if (family == "gaussian") {
          resid <- outer(Zo_cc, mu, "-")
          w.it * resid
        } else if (family == "poisson") {
          mu_val <- exp(mu)
          w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-mu))
          w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
        }
      }

      H.3.inv <- tryCatch(
        diag(sigma2 / w_colsums),
        error = function(e) {
          warning("Hessian inversion failed. SEs will be NA.")
          matrix(NA_real_, iT, iT)
        }
      )
      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    }
  } else {
    if (is.null(p.zx)) {
      stop("p.zx must be provided for ML distal outcome estimation.")
    }

    # ML score helper
    # Three-step ML over the expanded data (R/distal-ml.R): the E-step
    # weights lambda_it = sum_s w_is R_ist average the per-record posteriors,
    # so proportional-assignment weights stay outside the log.
    #
    # For the gaussian family the common within-class variance sigma2 is
    # estimated jointly with the class means (Bakk, Tekle & Vermunt 2013:
    # normal distal outcome with constant error variance). Unlike ordinary
    # regression, sigma2 does not factor out of the mean estimates here: it
    # enters the posterior weights P(X = t | W_i, Zo_i), so fixing it would
    # bias the means whenever the true variance differs.
    gaussian <- family == "gaussian"
    theta_of <- function(mu, sigma2) if (gaussian) c(mu, sigma2) else mu
    records_at <- function(theta) distal_records(p.zx(theta), pi_s, pwx)
    derivs_at <- function(theta) distal_unit_derivs(theta, Zo_cc, iT, family)

    beta <- beta_init
    sigma2 <- if (gaussian) stats::var(Zo_cc) else NULL
    Z_long <- rep(Zo_cc, iT)
    X_long <- factor(rep(seq_len(iT), each = N))

    for (iter in seq_len(em.maxIter)) {
      theta <- theta_of(beta, sigma2)
      lambda <- distal_lambda(records_at(theta)$R, w.is_cc)

      fit <- glm(Z_long ~ X_long - 1, family = family, weights = as.vector(lambda))
      beta_new <- coef(fit)
      sigma2_new <- if (gaussian) {
        sum(lambda * outer(Zo_cc, beta_new, "-")^2) / sum(lambda)
      } else {
        NULL
      }

      converged <- abs(neg.ll(theta_of(beta_new, sigma2_new)) - neg.ll(theta)) <
        covariate.tol
      beta <- beta_new
      sigma2 <- sigma2_new
      if (converged) {
        break
      }
      if (iter == em.maxIter) {
        warning("ML distal EM reached maximum iterations.")
      }
    }
    theta_hat <- theta_of(beta, sigma2)
    P <- length(theta_hat)

    three_step.score <- function(params) {
      distal_score(distal_lambda(records_at(params)$R, w.is_cc), derivs_at(params)$G)
    }

    H.3.inv <- tryCatch(
      qr.solve(distal_neg_hessian(records_at(theta_hat), w.is_cc, derivs_at(theta_hat))),
      error = function(e) {
        warning("Hessian inversion failed. SEs will be NA.")
        matrix(NA_real_, P, P)
      }
    )
    res <- list(
      par = beta,
      theta = theta_hat,
      value = neg.ll(theta_hat),
      convergence = 0L,
      sigma2 = sigma2
    )
  }

  return(list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  ))
}
