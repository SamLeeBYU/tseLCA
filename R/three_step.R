#' One-hot expand an integer response matrix
#'
#' Converts an N x H matrix of 0-based integer category values into an
#' N x sum(ivItemcat) binary indicator matrix, one column per category per item.
#' @noRd
expand_Y <- function(mY_int, ivItemcat) {
  # mY_int: N x H matrix of integer category values (0-based)
  # ivItemcat: length-H vector of number of categories per item
  N <- nrow(mY_int)
  H <- ncol(mY_int)
  out <- matrix(0, N, sum(ivItemcat))
  col_start <- 1L
  for (h in seq_len(H)) {
    K_h <- ivItemcat[h]
    for (k in seq_len(K_h)) {
      out[, col_start + k - 1L] <- as.integer(mY_int[, h] == (k - 1L))
    }
    col_start <- col_start + K_h
  }
  out
}

#' Expand a compact mPhi to a full item-probability matrix
#'
#' Converts multilevLCA's storage convention (one row per dichotomous item,
#' K rows per polytomous item) into a sum(ivItemcat) x T expanded matrix where
#' each item block contains all K category probabilities including the reference.
#' @noRd
expand_Phi <- function(phi_mat, ivItemcat) {
  dichotomous <- ivItemcat == 2L
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    if (dichotomous[h]) {
      result[[h]] <- rbind(1 - phi_mat[h_phi, ], phi_mat[h_phi, ])
      h_phi <- h_phi + 1L
    } else {
      K_h <- ivItemcat[h]
      result[[h]] <- phi_mat[h_phi:(h_phi + K_h - 1L), , drop = FALSE]
      h_phi <- h_phi + K_h
    }
  }
  do.call(rbind, result)
}

#' Expand a free-parameter phi matrix to full category probabilities
#'
#' Inverse of the simplex constraint: given (K-1) free rows per polytomous item
#' and 1 row per dichotomous item, prepends the reference P(Y=0|C) row for each
#' polytomous item so that the result aligns with expand_Y output.
#' @noRd
expand_Phi_free <- function(phi_free, ivItemcat) {
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    if (K_h == 2L) {
      result[[h]] <- phi_free[h_phi, , drop = FALSE]
      h_phi <- h_phi + 1L
    } else {
      rows_h <- phi_free[h_phi:(h_phi + K_h - 2L), , drop = FALSE]
      result[[h]] <- rbind(1 - colSums(rows_h), rows_h)
      h_phi <- h_phi + K_h - 1L
    }
  }
  do.call(rbind, result)
}

#' Per-observation class log-likelihood matrix
#'
#' Returns an N x T matrix where entry `[i, t]` is the conditional log-likelihood
#' log P(Y_i | X=t) under the expanded item-probability matrix mPhi.
#' mDesign masks missing indicators (0 = missing, 1 = observed).
#' @noRd
log_lik_matrix <- function(Y, mPhi, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }
  (mDesign * Y) %*% log(mPhi)
}

#' Joint observed-data log-likelihood with covariates
#'
#' Computes sum_i log P(Y_i, Z_i) = sum_i log sum_t P(Y_i|X=t) P(X=t|Z_i)
#' under a multinomial logit structural model with coefficient matrix gamma.coefs
#' (Q x (T-1), reference class absorbed into the intercept column of Z).
#' @noRd
joint_log_lik <- function(Y, Z, mPhi, gamma.coefs, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }

  log_P_Y_given_X <- log_lik_matrix(Y, mPhi, mDesign)

  eta <- Z %*% gamma.coefs
  eta_full <- cbind(0, eta)

  row_maxes_eta <- apply(eta_full, 1, max)
  log_denom_eta <- row_maxes_eta + log(rowSums(exp(eta_full - row_maxes_eta)))
  log_P_X_given_Z <- eta_full - log_denom_eta

  log_joint_prob <- log_P_Y_given_X + log_P_X_given_Z

  row_maxes_joint <- apply(log_joint_prob, 1, max)
  log_marg_prob <- row_maxes_joint +
    log(rowSums(exp(log_joint_prob - row_maxes_joint)))

  sum(log_marg_prob)
}

#' Joint log-likelihood for the distal outcome model
#'
#' Computes sum_i log( sum_t P(X=t|Zp_i) * P(Zo_i|X=t) * P(Y_i|X=t) )
#' When Zp is absent, P(X=t|Zp) = vPi (flat prevalences).
#'
#' @param Y       N x K_total expanded one-hot response matrix.
#' @param Zo      Length-N distal outcome vector.
#' @param mPhi    K_total x T expanded item-response probability matrix.
#' @param p.zx    N x T matrix of log-densities log P(Zo_i|X=t).
#' @param pi_mat  N x T matrix of class priors P(X=t|Zp_i). If NULL, uses
#'   flat prevalences from the row means of p.zx (not used; vPi supplied).
#' @param vPi     Length-T flat prevalences, used when pi_mat is NULL.
#' @param mDesign N x K_total design matrix (NULL for complete data).
#' @noRd
joint_log_lik_distal <- function(
  Y,
  mPhi,
  log_pZo_t,
  pi_mat = NULL,
  vPi = NULL,
  mDesign = NULL
) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }

  # log P(Y_i | X=t): N x T
  log_P_Y_t <- log_lik_matrix(Y, mPhi, mDesign)

  # log P(X=t | Zp_i): N x T
  if (!is.null(pi_mat)) {
    log_P_X_t <- log(pmax(pi_mat, 1e-300))
  } else {
    # flat prevalences: broadcast vPi across rows
    log_P_X_t <- matrix(
      log(pmax(vPi, 1e-300)),
      nrow(Y),
      length(vPi),
      byrow = TRUE
    )
  }

  # log P(Zo_i | X=t): N x T  (passed in as log_pZo_t)
  log_joint <- log_P_X_t + log_pZo_t + log_P_Y_t # N x T

  row_max <- apply(log_joint, 1L, max)
  log_marg <- row_max + log(rowSums(exp(log_joint - row_max)))
  sum(log_marg)
}

#' Compute posterior class probabilities from the unconstrained theta1 vector
#'
#' Reconstructs vPi and phi from the stacked parameter vector theta1 =
#' `c(vPi[-1], phi_free)` used internally by lca_step2, then returns the
#' N x T soft posterior matrix P(X=t|Y_i).
#' @noRd
compute_posteriors <- function(Y, mDesign, theta1, ivItemcat, iT) {
  vPi_free <- theta1[1:(iT - 1L)]
  vPi <- c(1 - sum(vPi_free), vPi_free)
  n_free <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat - 1L))
  phi_free <- matrix(
    theta1[iT:(iT + n_free * iT - 1L)],
    nrow = n_free,
    ncol = iT
  )

  mPhi <- expand_Phi(expand_Phi_free(phi_free, ivItemcat), ivItemcat)

  log_joint <- sweep(log_lik_matrix(Y, mPhi, mDesign), 2, log(vPi), "+")
  row_max <- apply(log_joint, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_joint - row_max)))
  exp(log_joint - log_denom)
}

#' Posterior class probabilities from a Step-1 fit, in data-row order
#'
#' P(X = t | Y_i) for the rows of `Y.exp` (one-hot expanded, as returned by
#' clean_data()), using the item-response probabilities and class sizes in
#' `fit0` as they stand (i.e. after any rebase permutation).
#' @noRd
step1_posteriors <- function(Y.exp, mDesign, fit0, ivItemcat) {
  log_joint <- sweep(
    log_lik_matrix(Y.exp, expand_Phi(fit0$mPhi, ivItemcat), mDesign),
    2L,
    log(fit0$vPi),
    "+"
  )
  row_max <- apply(log_joint, 1L, max)
  post <- exp(log_joint - (row_max + log(rowSums(exp(log_joint - row_max)))))
  colnames(post) <- paste0("C", seq_len(ncol(post)))
  post
}

#' Step-1 sample used for the measurement-model variance (Sigma.1)
#'
#' Prefers the data stored with the measurement model at fit time
#' (`s1$Y.exp`, `s1$mDesign.exp`, in data-row order). Falls back to decoding
#' multilevLCA's `fit0$mU` for measurement models that do not carry their
#' data (e.g. raw lca_step1() output); `mU` is sorted by response pattern,
#' which is harmless here because its rows are used consistently. Returns
#' NULL if neither is available.
#' @noRd
step1_sample <- function(s1, ivItemcat, ref_idx = 1L) {
  if (!is.null(s1$Y.exp)) {
    return(list(
      Y.exp = s1$Y.exp,
      mDesign = s1$mDesign.exp,
      ivItemcat = ivItemcat,
      u_post = NULL
    ))
  }
  fit0 <- s1$fit0
  if (is.null(fit0$mU)) {
    return(NULL)
  }
  raw <- extract_Y_from_mU(fit0, ivItemcat)
  if (ref_idx != 1L) {
    iT <- length(fit0$vPi)
    raw$u_post <- raw$u_post[, c(ref_idx, seq_len(iT)[-ref_idx]), drop = FALSE]
  }
  raw
}

#' Compute classification-error matrix with optional covariate-adjusted prior
#'
#' Returns posteriors, modal/soft assignments (w.is), and the T x T
#' classification-error probability matrix p.wx_mat = P(W=s|X=t).
#' When pi_adj (N x T) is supplied, uses person-specific class priors from the
#' covariate model; otherwise falls back to the flat vPi from fit0.
#' @noRd
compute_pwx_adj <- function(
  Y.obs,
  fit0,
  ivItemcat,
  mDesign = NULL,
  use.modal.assignment = TRUE,
  pi_adj = NULL # N x T covariate-adjusted class probs, or NULL for flat vPi
) {
  N <- nrow(Y.obs)
  iT <- ncol(fit0$mPhi)

  mPhi_exp <- expand_Phi(fit0$mPhi, ivItemcat)
  phi_clamped <- pmax(pmin(mPhi_exp, 1 - 1e-10), 1e-10)

  log_p_it <- if (is.null(mDesign)) {
    Y.obs %*% log(phi_clamped)
  } else {
    (mDesign * Y.obs) %*% log(phi_clamped)
  }

  # use adjusted or flat priors
  if (!is.null(pi_adj)) {
    log_prior <- log(pi_adj) # N x T, person-specific
  } else {
    log_prior <- matrix(log(fit0$vPi), nrow = N, ncol = iT, byrow = TRUE)
  }

  log_joint <- log_p_it + log_prior # N x T
  row_max <- apply(log_joint, 1, max)
  post <- exp(log_joint - row_max - log(rowSums(exp(log_joint - row_max)))) # N x T posteriors

  w.is <- if (use.modal.assignment) {
    w <- matrix(0L, N, iT)
    w[cbind(seq_len(N), max.col(post))] <- 1L
    w
  } else {
    post
  }

  p.wx_joint <- (t(w.is) %*% post) / N
  p.wx_mat <- sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")

  list(
    post = post,
    w.is = w.is,
    p.wx_mat = p.wx_mat
  )
}

# -- Step 2: Posteriors and classification-error matrix -----------------------

#' Step 2: posteriors, classification-error matrix, and Jacobian closure
#'
#' Computes all Step-2 quantities needed for Step 3 and variance propagation:
#' theta1 and theta2 (constrained and unconstrained parameterizations of the
#' classification-error matrix), w.is (modal or soft assignments), p.wx_mat,
#' and optionally a closure compute_J_unc for the analytic Jacobian
#' d theta2 / d u used in the measurement-uncertainty correction.
#' Returns NULL for compute_J_unc when use.simple.cov = TRUE.
#' @noRd
lca_step2 <- function(
  Y.obs,
  fit0,
  n_classes,
  use.modal.assignment,
  boundary.tol,
  use.simple.cov,
  ivItemcat,
  mDesign = NULL
) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y.obs), ncol(Y.obs))
  }

  N <- nrow(Y.obs)
  iT <- n_classes
  K <- ncol(Y.obs) #sum(K_h)

  # Number of free rows per item in mPhi
  n_free <- ivItemcat - 1L

  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )

  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))

  phi_free <- fit0$mPhi[free_idx, ]

  theta1 <- c(
    fit0$vPi[2:iT],
    phi_free
  )

  p.xy <- compute_posteriors(Y.obs, mDesign, theta1, ivItemcat, iT)
  assignment <- max.col(p.xy)

  make_w <- function(posteriors) {
    if (use.modal.assignment) {
      w <- matrix(0L, nrow = N, ncol = iT)
      w[cbind(seq_len(N), max.col(posteriors))] <- 1L
      w
    } else {
      posteriors
    }
  }

  w.is <- make_w(p.xy)

  compute_pwx <- function(t1) {
    post <- compute_posteriors(Y.obs, mDesign, t1, ivItemcat, iT)
    w_local <- make_w(post)
    p.wx_joint <- (t(w_local) %*% post) / N
    sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")
  }

  theta2_from_theta1 <- function(th1) {
    # rho <- c(1 - sum(th1[1:(iT - 1)]), th1[1:(iT - 1)])
    # phi <- matrix(th1[iT:length(th1)], nrow = K, ncol = iT)
    p.wx_mat <- compute_pwx(th1)
    log_ref <- log(diag(p.wx_mat))
    gamma_mat <- sweep(log(p.wx_mat), 2, log_ref, "-")
    gamma_mat[row(gamma_mat) != col(gamma_mat)]
  }

  gamma_vec_to_pwx <- function(gamma_vec) {
    gamma_mat <- matrix(0, nrow = iT, ncol = iT)
    gamma_mat[row(gamma_mat) != col(gamma_mat)] <- gamma_vec
    exp_mat <- exp(gamma_mat)
    sweep(exp_mat, 2, colSums(exp_mat), "/")
  }

  theta2 <- theta2_from_theta1(theta1)
  p.wx_mat <- gamma_vec_to_pwx(theta2)

  if (!use.simple.cov) {
    compute_J_unc_analytical <- function(
      p_ik,
      Y_obs,
      mDes,
      th1,
      ivItemcat,
      T_classes
    ) {
      N <- nrow(p_ik)

      A <- t(p_ik) %*% p_ik
      A[A < 1e-12] <- 1e-12

      #Extract item probabilities to match the free parameter structure
      phi_mat <- matrix(
        th1[T_classes:length(th1)],
        nrow = sum(ivItemcat - 1L),
        ncol = T_classes
      )

      L_rho <- T_classes - 1L
      L_phi <- sum((ivItemcat - 1L) * T_classes)
      L <- L_rho + L_phi

      J <- matrix(0, nrow = T_classes * (T_classes - 1L), ncol = L)

      #Offsets for locating items and categories in the expanded matrices
      starts_Y <- c(1L, cumsum(ivItemcat)[-length(ivItemcat)] + 1L)
      item_offsets <- c(
        0L,
        cumsum((ivItemcat - 1L) * T_classes)[-length(ivItemcat)]
      )

      # Map (s, t) to the correct row in J (matching gamma_mat[row != col] column-major)
      st_idx <- 1L
      st_map <- matrix(0L, nrow = T_classes, ncol = T_classes)
      for (t in seq_len(T_classes)) {
        for (s in seq_len(T_classes)) {
          if (s != t) {
            st_map[s, t] <- st_idx
            st_idx <- st_idx + 1L
          }
        }
      }

      for (t in seq_len(T_classes)) {
        P_tt <- (p_ik[, t]^2) / A[t, t]

        for (s in seq_len(T_classes)) {
          if (s == t) {
            next
          }
          row_J <- st_map[s, t]
          P_st <- (p_ik[, s] * p_ik[, t]) / A[s, t]

          for (c_prime in seq_len(T_classes)) {
            I_s <- if (s == c_prime) 1.0 else 0.0
            I_t <- if (t == c_prime) 1.0 else 0.0

            #shared derivative component for class c_prime
            Q_stc <- P_st *
              (I_s + I_t - 2 * p_ik[, c_prime]) -
              2 * P_tt * (I_t - p_ik[, c_prime])

            sum_Q <- sum(Q_stc)

            #Derivative for class prevalence u^rho (only c' >= 2)
            if (c_prime >= 2L) {
              col_rho <- c_prime - 1L
              J[row_J, col_rho] <- sum_Q
            }

            #Derivative for item response u^phi
            for (h in seq_along(ivItemcat)) {
              K_h <- ivItemcat[h]
              n_free <- K_h - 1L
              Y_cols <- (starts_Y[h] + 1L):(starts_Y[h] + K_h - 1L)

              phi_row_start <- if (h == 1L) {
                1L
              } else {
                sum(ivItemcat[1:(h - 1L)] - 1L) + 1L
              }
              phi_vals <- phi_mat[
                phi_row_start:(phi_row_start + n_free - 1L),
                c_prime
              ]

              col_J_start <- L_rho + item_offsets[h] + (c_prime - 1L) * n_free

              for (k in seq_len(n_free)) {
                Y_col <- Y_cols[k]

                val <- sum(Q_stc * Y_obs[, Y_col]) -
                  phi_vals[k] * sum(Q_stc * mDes[, Y_col])
                J[row_J, col_J_start + k] <- val
              }
            }
          }
        }
      }
      J
    }
  }

  list(
    theta1 = theta1,
    theta2 = theta2,
    w.is = w.is,
    p.wx_mat = p.wx_mat,
    gamma_vec_to_pwx = gamma_vec_to_pwx,
    theta2_from_theta1 = theta2_from_theta1,
    p.xy = p.xy,
    compute_J_unc = if (!use.simple.cov) compute_J_unc_analytical else NULL
  )
}

#' BCH classification-error-corrected weight matrix
#'
#' Computes the N x T BCH weight matrix used throughout the BCH estimators:
#' \code{w.it = w.is \%*\% t(pwx)^-1}, where \code{pwx[s, t] = P(W = s | X =
#' t)} is the column-stochastic classification-error matrix from
#' \code{compute_pwx_adj()}/\code{lca_step2()} (\code{colSums(pwx) == 1}).
#' @noRd
bch_weight_matrix <- function(w.is, pwx) {
  w.is %*% t(qr.solve(pwx))
}

#' Moore-Penrose pseudo-inverse and numerical rank with SVD
#'
#' Used by the omnibus class-equality Wald test (`omnibus_test()`), whose
#' contrast covariance is rank-deficient for the multinomial family (each
#' class's C-vector of category probabilities sums to 1, so a difference of
#' two classes' full probability vectors always sums to 0 across
#' categories) and may be for other families too under boundary/near-
#' collinear fits. `qr()`-based rank/solve is avoided because it is less
#' numerically stable than SVD for a covariance matrix that is exactly
#' singular by construction, not just ill-conditioned.
#' @noRd
pinv_rank <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  keep <- s$d > (tol * max(s$d))
  d_inv <- ifelse(keep, 1 / s$d, 0)
  list(pinv = s$v %*% (d_inv * t(s$u)), rank = sum(keep))
}

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

#' Step 3 (covariate): estimate multinomial logit gamma with either BCH or ML EM
#'
#' Optimizes the Q x (T-1) coefficient matrix gamma for P(X=t|Z_i) with
#' Newton-Raphson (BCH) or EM with an inner NR M-step (ML). Returns the
#' parameter vector, the inverted Hessian H.3.inv (or NA matrix on failure),
#' used by lca_vcov for sandwich variance propagation.
#' @noRd
lca_step3 <- function(
  neg.ll,
  gamma_init,
  Q,
  iT,
  covariate.tol,
  use.bch = FALSE,
  gradient = NULL,
  Z_mat_cc = NULL,
  w.is_cc = NULL,
  p.xz = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  verbose = FALSE,
  correct.spec = FALSE
) {
  N <- nrow(Z_mat_cc)
  beta <- matrix(gamma_init, nrow = Q, ncol = iT - 1)
  ll_prev <- -neg.ll(c(beta))
  H <- NULL
  # print(ll_prev)
  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx) # N x T
    w.it_plus <- rowSums(w.it)

    for (nr in seq_len(em.maxIter)) {
      # print(beta)
      # print(-neg.ll(c(beta)))
      grad_vec <- -gradient(c(beta)) # gradient of pos. ll: Q*(T-1) vector

      H <- matrix(0, Q * (iT - 1), Q * (iT - 1))
      pi_ <- p.xz(beta)
      for (k in seq_len(iT - 1)) {
        for (l in k:(iT - 1)) {
          w_kl <- w.it_plus * pi_[, k + 1L] * ((k == l) - pi_[, l + 1L])
          idx_k <- ((k - 1) * Q + 1):(k * Q)
          idx_l <- ((l - 1) * Q + 1):(l * Q)
          block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
          H[idx_k, idx_l] <- block
          if (k != l) H[idx_l, idx_k] <- t(block) # Clairaut: H symmetric
        }
      }

      if (
        nr >= max(em.maxIter / 5, 1) && #Wait a little bit before testing for PSD
          inherits(
            tryCatch(chol(H), error = function(e) e),
            "error"
          )
      ) {
        stop(
          sprintf(
            "BCH Newton-Raphson failed after %d iterations: Hessian is not positive semi-definite. ",
            nr
          ),
          "This typically occurs under low class separation. ",
          "Try use.bch = FALSE to use the ML estimator instead. Or you can try increasing em.maxIter...",
          call. = FALSE
        )
      }

      direction <- tryCatch(
        qr.solve(-H, grad_vec),
        error = function(e) rep(0, Q * (iT - 1))
      )

      alpha <- 1
      current_ll <- -neg.ll(c(beta))

      beta_vec <- c(beta)

      while (alpha > (covariate.tol / 2)) {
        trial <- beta_vec + alpha * direction
        trial_ll <- tryCatch(-neg.ll(trial), error = function(e) -Inf)
        if (is.finite(trial_ll) && trial_ll > current_ll) {
          break
        }
        alpha <- alpha / 2
      }

      delta <- alpha * direction
      #print(max(abs(delta)))
      beta <- beta + matrix(delta, nrow = Q, ncol = iT - 1)
      if (max(abs(delta)) < covariate.tol) {
        if (verbose) {
          message(sprintf("BCH NR converged in %d iterations.", nr))
        }
        break
      }
      if (nr == em.maxIter) stop("BCH NR reached maximum iterations.")
    }
    #print(H)
  } else {
    for (iter in seq_len(em.maxIter)) {
      #print(beta)
      #print(ll_prev)

      # E-step ######################################################################
      p <- p.xz(beta)

      q <- p %*% t(pwx)

      gamma <- matrix(0, nrow = N, ncol = iT)
      for (t in seq_len(iT)) {
        gamma[, t] <- p[, t] * rowSums(w.is_cc * outer(rep(1, N), pwx[, t]) / q)
      }
      ###############################################################################

      #M-step ########################################################################
      gamma_plus <- rowSums(gamma)
      Gamma_nr <- gamma[, -1, drop = FALSE]

      for (nr in seq_len(10L)) {
        p_nr <- p.xz(beta)
        p_nr1 <- p_nr[, -1, drop = FALSE]

        grad <- t(Z_mat_cc) %*% (Gamma_nr - p_nr1 * gamma_plus)

        H <- matrix(0, Q * (iT - 1), Q * (iT - 1))
        for (k in seq_len(iT - 1)) {
          for (l in k:(iT - 1)) {
            w_kl <- gamma_plus * p_nr1[, k] * ((k == l) - p_nr1[, l])
            idx_k <- ((k - 1) * Q + 1):(k * Q)
            idx_l <- ((l - 1) * Q + 1):(l * Q)
            block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
            H[idx_k, idx_l] <- block
            # Clairaut: H symmetric
            if (k != l) H[idx_l, idx_k] <- t(block)
          }
        }

        delta <- tryCatch(
          qr.solve(-H, as.vector(grad)),
          error = function(e) rep(0, Q * (iT - 1))
        )
        beta <- beta + matrix(delta, nrow = Q, ncol = iT - 1)
        if (max(abs(delta)) < covariate.tol) break
      }
      if (nr == 10L && max(abs(delta)) >= covariate.tol) {
        warning(sprintf(
          "M-step in EM algorithm did not converge at iteration %d",
          iter
        ))
      }

      #Alternatively, fit a multinomial logistic regression model ##########################
      # class_exp <- rep(seq_len(iT), each = N)
      # Z_exp <- Z_mat_cc[rep(seq_len(N), iT), , drop = FALSE]
      # w_exp <- as.vector(gamma)

      # fit <- nnet::multinom(
      #   class_exp ~ Z_exp - 1,
      #   weights = w_exp,
      #   trace = FALSE,
      #   maxit = 500L
      # )

      # beta <- t(coef(fit))
      ######################################################################################

      ll_curr <- -neg.ll(c(beta))
      if (abs(ll_curr - ll_prev) < covariate.tol && iter > 1L) {
        if (verbose) {
          message(sprintf("EM converged in %d iterations.", iter))
        }
        break
      }
      ll_prev <- ll_curr

      if (iter == em.maxIter) {
        warning("EM reached maximum iterations without converging.")
      }
    }
  }
  res <- list(par = c(beta), value = neg.ll(c(beta)), convergence = 0L)

  H.3.inv <- tryCatch(
    {
      if (!use.bch) {
        if (correct.spec) {
          matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
        } else {
          # -- Analytic observed-data Hessian of neg.ll (checked with sympy)--------------------------
          # neg.ll = -sum_i sum_s w_{is} * log(q_{is})
          # q_{is} = sum_t pi_{it}(beta) * pwx[s,t]
          # r_{is} = w_{is} / q_{is}
          #
          # H_{(q,k),(p,l)} = sum_i z_{iq}*z_{ip} * [
          #   pi_{i,k+1}*(I(k==l)-pi_{i,l+1}) * F_k
          # - pi_{i,k+1} * pi_{i,l+1} * G_{kl}
          # ]
          # where:
          #   F_k  = sum_s w_{is}*pwx[s,k+1]/q_{is} - sum_s w_{is}
          #   G_{kl} = sum_s w_{is}*pwx[s,k+1]*(pwx[s,l+1]-q_{is})/q_{is}^2

          p_ <- p.xz(beta) # N x T
          q_mat <- p_ %*% t(pwx) # N x T: q[i,s]
          r_mat <- w.is_cc / q_mat # N x T: r[i,s]

          # F_k for each non-reference class k (N x (T-1))
          F_mat <- matrix(0, N, iT - 1L)
          for (k in seq_len(iT - 1L)) {
            F_mat[, k] <- rowSums(r_mat * pwx[, k + 1L][col(r_mat)]) -
              rowSums(w.is_cc)
          }

          # G_{kl} for each pair (k,l)
          H_obs <- matrix(0, Q * (iT - 1L), Q * (iT - 1L))
          for (k in seq_len(iT - 1L)) {
            for (l in k:(iT - 1L)) {
              idx_k <- ((k - 1L) * Q + 1L):(k * Q)
              idx_l <- ((l - 1L) * Q + 1L):(l * Q)
              tA <- p_[, k + 1L] * ((k == l) - p_[, l + 1L]) * F_mat[, k]
              G_kl <- rowSums(
                w.is_cc *
                  pwx[, k + 1L][col(w.is_cc)] *
                  (pwx[, l + 1L][col(w.is_cc)] - q_mat) /
                  q_mat^2
              )
              tB <- -p_[, k + 1L] * p_[, l + 1L] * G_kl

              block <- -t(Z_mat_cc) %*% ((tA + tB) * Z_mat_cc) # Q x Q
              H_obs[idx_k, idx_l] <- block
              if (k != l) {
                H_obs[idx_l, idx_k] <- t(block)
              } # Clairaut: H is symmetric
            }
          }
          qr.solve(H_obs)
        }
      } else {
        qr.solve(-H)
      }
    },
    error = function(e) {
      if (!use.bch) {
        warning(
          "Exact Hessian inversion failed. Falling back to Hessian estimated from the cross-product of the case-wise log-likelihood gradients. This assumes a correct third-stage model specification."
        )
        matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
      } else {
        warning(
          "Hessian inversion failed after BCH. Try again with use.bch = FALSE."
        )
        matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
      }
    }
  )
  return(list(res = res, H.3.inv = H.3.inv))
}


# -- Variance estimation (Bakk et al., 2014) ----------------------------------

#' Individual-level BHHH variance matrix for binary and polytomous LCA
#'
#' Computes the outer-product (BHHH) information matrix and variance-covariance
#' matrix for LCA measurement model parameters in the unconstrained
#' (logit/log-ratio) space, matching \pkg{multilevLCA}'s \code{$Varmat}.
#'
#' The score in unconstrained space is
#' \eqn{s_{it} = u_{it}(y_i - d_i \circ p_{it})},
#' where \eqn{d_i} is the missing-data design indicator matrix.
#'
#' Assumes \code{fit0$mPhi} follows the \pkg{multilevLCA} storage convention:
#' \itemize{
#'   \item Dichotomous item h (\code{ivItemcat[h] == 2}): 1 row =
#'     \eqn{P(Y=1|C)}; the base level \eqn{P(Y=0|C)} is excluded.
#'   \item Polytomous item h (\code{ivItemcat[h] > 2}): \code{K_h} rows =
#'     \eqn{P(Y=0|C), \ldots, P(Y=K_h-1|C)}; the base level is included.
#' }
#' \code{expand_Y} produces one-hot columns in the same order so that
#' \code{expand_Phi(fit0$mPhi, ivItemcat)} aligns column-wise with
#' \code{expand_Y(mY, ivItemcat)}.  Free (estimable) parameters per item are
#' the single \eqn{P(Y=1|C)} row for dichotomous items, and rows 2 through
#' \eqn{K_h} for polytomous items (row 1, \eqn{P(Y=0|C)}, is the reference).
#' Boundary parameters (within \code{boundary.tol} of 0 or 1) are treated as
#' fixed: their score columns are zeroed and they do not contribute to the
#' information matrix.
#'
#' @param Y.exp       N x sum(K_h) expanded one-hot indicator matrix.
#' @param mDesign.exp Expanded design matrix (same dimensions as \code{Y.exp}),
#'   or \code{NULL} for complete data.
#' @param fit0        Step-1 fit object with \code{$vPi} and \code{$mPhi}.
#' @param ivItemcat   Integer vector of category counts per item.
#' @param boundary.tol Scalar tolerance for boundary detection. Default
#'   \code{1e-2}.
#' @param use.freq    Logical. Collapse duplicate score rows before computing
#'   the cross-product, weighting by frequency. Default \code{TRUE}.
#' @param u_post      Optional N x T matrix of posterior class probabilities.
#'   When supplied (e.g. extracted from \code{fit0$mU} with
#'   \code{extract_Y_from_mU}), \code{compute_posteriors} is skipped.
#'   Default \code{NULL}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`Infomat`}{Square BHHH information matrix of dimension p x p,
#'       where p = (iT-1) + sum(ivItemcat - 1) * iT is the total number of free
#'       parameters. Boundary parameters have zero rows and columns.}
#'     \item{`Varmat`}{Inverse of \code{Infomat} divided by N, giving the
#'       asymptotic variance-covariance matrix on the same scale as
#'       \pkg{multilevLCA}'s \code{$Varmat}. Boundary parameters have zero
#'       rows and columns.}
#'     \item{`SEs`}{Numeric vector of length p. Square root of the diagonal of
#'       \code{Varmat}; zero for boundary parameters.}
#'     \item{`mScore`}{N x p matrix of individual score contributions in the
#'       unconstrained parameterization, used for sandwich variance propagation
#'       in \code{lca_vcov} and \code{lca_vcov_distal}.}
#'   }
#' @keywords internal
lca_indiv_varmat <- function(
  Y.exp,
  mDesign.exp,
  fit0,
  ivItemcat,
  boundary.tol = 1e-2,
  use.freq = TRUE,
  u_post = NULL
) {
  pi_ <- fit0$vPi
  phi <- fit0$mPhi
  iT <- length(pi_)
  N <- nrow(Y.exp)

  if (is.null(mDesign.exp)) {
    mDesign.exp <- matrix(1L, N, ncol(Y.exp))
  }

  #Protect against parameter estimates on the boundary of the support (zero out their score contributions)
  pi_bdry <- pi_ <= boundary.tol | pi_ >= (1 - boundary.tol)
  phi_bdry <- phi <= boundary.tol | phi >= (1 - boundary.tol)

  #Clamp boundary parameters before computing posteriors
  pi_[pi_bdry] <- pmax(pmin(pi_[pi_bdry], 1 - 1e-6), 1e-6)
  phi[phi_bdry] <- pmax(pmin(phi[phi_bdry], 1 - 1e-6), 1e-6)

  # ---- Build theta1 and compute posteriors -----------------------------------
  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )
  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))
  phi_free <- phi[free_idx, , drop = FALSE]

  # Use pre-computed posteriors when available (e.g. extracted from fit0$mU),
  # otherwise compute them from theta1.
  if (is.null(u_post)) {
    theta1 <- c(pi_[-1L], phi_free)
    u_post <- compute_posteriors(Y.exp, mDesign.exp, theta1, ivItemcat, iT)
  }

  # ---- Expand phi for residual computation -----------------------------------
  phi_exp <- expand_Phi(phi, ivItemcat) # K_total x T

  # ---- Build free_cols (columns of phi_exp for free parameters) -------------
  # free_cols: indices into phi_exp rows for the free categories.
  # phi_exp has K_total rows; free categories are the non-reference rows
  # within each item block. The mapping from free_idx (mPhi rows) to
  # phi_exp rows differs between dichotomous and polytomous items:
  #
  #   Dichotomous (K=2): phi_exp block = [P(Y=0), P(Y=1)], 2 rows.
  #     free_idx points to the single mPhi row = P(Y=1) = phi_exp row 2
  #     within the block (col_start + 1).
  #
  #   Polytomous (K>2): phi_exp block = [P(Y=0)..P(Y=K-1)], K rows.
  #     free_idx points to mPhi rows 2..K within the block (P(Y=1)..P(Y=K-1))
  #     = phi_exp rows col_start+1 .. col_start+K-1.
  #
  # In both cases the free phi_exp rows are exactly col_start+1..col_start+K-1
  # (dropping col_start = reference P(Y=0) for binary and polytomous alike).

  free_cols <- integer(0L)
  col_start <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    free_cols <- c(free_cols, (col_start + 1L):(col_start + K_h - 1L))
    col_start <- col_start + K_h
  }
  n_free_phi <- length(free_cols) # = sum(ivItemcat - 1) = nrow(phi_free) for all items

  # phi_bdry_free: n_free_phi x T, boundary flags for free phi parameters.
  # free_idx already selects the free mPhi rows in item order, so:
  phi_bdry_free <- phi_bdry[free_idx, , drop = FALSE]

  # ---- Pi scores (T-1 columns, one per non-reference class t=2..T) -----------
  s_u_pi <- sweep(u_post[, -1L, drop = FALSE], 2L, pi_[-1L], "-")
  pi_free_bdry <- pi_bdry[-1L] # drop reference class t=1
  if (any(pi_free_bdry)) {
    s_u_pi[, pi_free_bdry] <- 0
  }

  # ---- Phi scores (n_free_phi * T columns, class-major) ----------------------
  Y_free <- Y.exp[, free_cols, drop = FALSE]
  D_free <- mDesign.exp[, free_cols, drop = FALSE]

  s_u_phi <- matrix(0, N, n_free_phi * iT)

  for (t in seq_len(iT)) {
    idx <- ((t - 1L) * n_free_phi + 1L):(t * n_free_phi)
    resid <- Y_free -
      D_free *
        matrix(
          phi_exp[free_cols, t],
          nrow = N,
          ncol = n_free_phi,
          byrow = TRUE
        )
    s_col <- u_post[, t] * resid

    #Zero boundary free parameters for this class
    bdry_t <- phi_bdry_free[, t]
    if (any(bdry_t)) {
      s_col[, bdry_t] <- 0
    }

    s_u_phi[, idx] <- s_col
  }

  S <- cbind(s_u_pi, s_u_phi) # N x p

  # ---- Identify active (non-boundary) columns --------------------------------
  # Boundary parameters have all-zero score columns to avoid rank deficiency, then restore zero rows/cols after
  active <- which(colSums(S != 0) > 0L)
  p_full <- ncol(S)

  # ---- BHHH information matrix (on active columns only) ----------------------
  S_active <- S[, active, drop = FALSE]

  if (use.freq) {
    S_char <- apply(S_active, 1L, paste, collapse = "\r")
    uniq <- !duplicated(S_char)
    freq <- tabulate(match(S_char, S_char[uniq]))
    Infomat_active <- crossprod(S_active[uniq, , drop = FALSE] * sqrt(freq)) / N
  } else {
    Infomat_active <- crossprod(S_active) / N
  }

  Varmat_active <- tryCatch(
    qr.solve(Infomat_active) / N,
    error = function(e) {
      warning(
        "lca_indiv_varmat: Infomat is singular even after removing boundary ",
        "parameters; returning NA matrix. Check for near-empty classes.",
        call. = FALSE
      )
      matrix(NA_real_, length(active), length(active))
    }
  )

  # ---- Restore full-size Infomat and Varmat ----------------------------------
  Infomat <- matrix(0, p_full, p_full)
  Infomat[active, active] <- Infomat_active

  Varmat <- matrix(0, p_full, p_full)
  Varmat[active, active] <- Varmat_active

  list(
    Infomat = Infomat,
    Varmat = Varmat,
    SEs = sqrt(diag(Varmat)),
    mScore = S
  )
}

#' Covariate model variance-covariance with measurement-uncertainty correction
#'
#' Assembles the sandwich variance matrix for the Step-3 gamma estimates,
#' optionally propagating Step-1 measurement uncertainty through the analytic
#' `C_mat = d/d(theta2) [sum_i score_3_i]` and Jacobian `J.2 = d(theta2)/d(u)`.
#' When use.simple.cov = TRUE, returns the plain robust sandwich H^{-1} S'S H^{-1}.
#' @noRd
lca_vcov <- function(
  coefs,
  three_step.score,
  H.3.inv,
  Sigma.1,
  theta2,
  J.2,
  p.wx_mat,
  w.is,
  Z_mat,
  n_classes,
  p.xz,
  s2,
  use.simple.cov
) {
  J.3 <- three_step.score(c(coefs))

  Sigma.3.robust <- H.3.inv %*% crossprod(J.3) %*% H.3.inv
  if (!use.simple.cov) {
    # -- Analytic C_mat = d/d theta2 [colSums(score_3)] (checked with sympy) ----------------------------

    T_ <- n_classes
    pwx <- p.wx_mat
    pi_ <- p.xz(matrix(coefs, ncol = T_ - 1L)) # N x T
    N <- nrow(Z_mat)
    Q_ <- ncol(Z_mat)

    q <- pi_ %*% t(pwx) # N x T: q[i,s]
    V <- s2$w.is / q # N x T: w[i,s]/q[i,s]
    U <- s2$w.is / q^2 # N x T: w[i,s]/q[i,s]^2

    # VP[i, t0] = sum_s V[i,s] * pwx[s, t0]
    VP <- V %*% pwx

    n_theta2 <- T_ * (T_ - 1L)
    n_coef <- Q_ * (T_ - 1L)
    C_mat <- matrix(0, nrow = n_coef, ncol = n_theta2)

    # Map (s0, t0) mapping to the correct column index in C_mat
    theta2_idx <- matrix(0, T_, T_)
    theta2_idx[row(theta2_idx) != col(theta2_idx)] <- seq_len(n_theta2)

    for (t0 in seq_len(T_)) {
      for (s0 in seq_len(T_)) {
        if (s0 == t0) {
          next
        }
        idx_theta2 <- theta2_idx[s0, t0]

        # Store derivatives for all classes k for this specific (s0, t0)
        d_score_mat <- matrix(0, nrow = N, ncol = T_ - 1L)

        for (k in seq_len(T_ - 1L)) {
          # UP[i] = sum_s U[i,s] * pwx[s, k+1] * pwx[s, t0]
          UP_k_t0 <- rowSums(
            U *
              matrix(
                pwx[, k + 1] * pwx[, t0],
                nrow = N,
                ncol = T_,
                byrow = TRUE
              )
          )

          term1 <- 0
          if (k + 1L == t0) {
            term1 <- pi_[, k + 1] * (V[, s0] - VP[, t0])
          }

          term2 <- pi_[, k + 1] *
            pi_[, t0] *
            (U[, s0] * pwx[s0, k + 1] - UP_k_t0)

          d_score_mat[, k] <- pwx[s0, t0] * (term1 - term2)
        }

        C_mat[, idx_theta2] <- as.vector(crossprod(Z_mat, d_score_mat))
      }
    }

    step1.uncertainty <- C_mat %*% J.2 %*% Sigma.1 %*% t(J.2) %*% t(C_mat)
    Sigma.3 <- H.3.inv %*% (crossprod(J.3) + step1.uncertainty) %*% H.3.inv
  } else {
    Sigma.3 <- Sigma.3.robust
    step1.uncertainty <- NULL
  }

  param_names <- as.vector(outer(
    rownames(coefs),
    colnames(coefs),
    paste,
    sep = ":"
  ))
  rownames(Sigma.3) <- param_names
  colnames(Sigma.3) <- param_names

  Sigma.3
}

#' Distal outcome variance-covariance with full uncertainty propagation
#'
#' Assembles the T x T sandwich variance matrix for the distal outcome mu
#' estimates, propagating Step-1 measurement uncertainty
#' and, when both a covariate and distal model are fitted, Step-3 covariate
#' uncertainty. Skips steps 1/2 uncertainty corrections when
#' use.bch = TRUE or use.simple.cov = TRUE.
#' @noRd
lca_vcov_distal <- function(
  mu_hat,
  three_step.score,
  pi_adj,
  w.is,
  p.wx_mat,
  p.zx,
  family,
  H.3.inv,
  Sigma.1,
  s2,
  Sigma.3 = NULL,
  s3.par = NULL,
  p.xz.cov = NULL,
  Z_mat_cov = NULL,
  iT,
  use.simple.cov,
  use.bch,
  unit_scores = NULL
) {
  # mu_hat is the full Step-3 parameter vector: the T class parameters, plus
  # sigma2 for the gaussian ML model. `unit_scores(mu_hat)`, when supplied,
  # returns a list with one N x T matrix per parameter p holding
  # d log f(z_i | X = t) / d theta_p; otherwise each class parameter only
  # moves its own class and the unit scores are recovered from the score.
  J.3 <- three_step.score(mu_hat)
  meat <- crossprod(J.3)
  P <- ncol(J.3)
  par_names <- c(paste0("mu_C", seq_len(iT)), if (P > iT) "sigma2")

  if (use.bch || use.simple.cov) {
    result <- H.3.inv %*% meat %*% H.3.inv
    dimnames(result) <- list(par_names, par_names)
    return(result)
  }

  if (is.null(unit_scores)) {
    stop("lca_vcov_distal(): `unit_scores` is required for ML uncertainty propagation.")
  }
  # Per-record posteriors and unit scores at the estimates; the Step-1 (C1)
  # and Step-2 (C_mat) cross-derivatives are derived in R/distal-ml.R.
  rec <- distal_records(p.zx(mu_hat), pi_adj, p.wx_mat)
  cross <- distal_cross_derivs(
    rec,
    w.is,
    p.wx_mat,
    unit_scores(mu_hat),
    Z_mat = if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
      Z_mat_cov
    } else {
      NULL
    }
  )

  step1.uncertainty <- cross$C1 %*%
    s2$J.2 %*%
    Sigma.1 %*%
    t(s2$J.2) %*%
    t(cross$C1)

  step2.uncertainty <- if (!is.null(cross$C_mat)) {
    cross$C_mat %*% Sigma.3 %*% t(cross$C_mat)
  } else {
    matrix(0, P, P)
  }

  result <- H.3.inv %*%
    (meat + step1.uncertainty + step2.uncertainty) %*%
    H.3.inv
  dimnames(result) <- list(par_names, par_names)
  result
}

#' Distal outcome variance-covariance for the multinomial family
#'
#' Multinomial analog of \code{lca_vcov_distal()}: assembles the
#' \code{(iT*C) x (iT*C)} sandwich variance matrix for the flattened T x C
#' class-conditional probability matrix \code{pi_hat}
#' (\code{matrix(theta_hat, nrow = iT, ncol = C)} recovers it), propagating
#' Step-1 measurement uncertainty (\code{C1_mat}/\code{step1.uncertainty})
#' and, when both a covariate and distal model are fitted, Step-3 covariate
#' uncertainty (\code{C_mat}/\code{step2.uncertainty}) when
#' \code{use.bch = FALSE} and \code{use.simple.cov = FALSE}. The
#' generalization from scalar \code{mu_t} (one parameter per class, as in
#' \code{lca_vcov_distal()}) to \code{pi_hat[t, ]} (C parameters per class)
#' only touches the "unit score" \code{g_it}: instead of dividing the
#' length-\code{iT} score by \code{r_it} once, the length-\code{iT*C} score
#' is divided by \code{r_it} replicated across the C categories, since
#' neither chain-rule term (\code{dr}, through \code{d ae/d theta2}; nor
#' \code{inner_t}, through \code{d r_it/d gamma}) depends on the category
#' dimension at all -- only on which class \code{t} a given column belongs
#' to. Both terms were cross-validated against independent numerical
#' differentiation of the case-wise estimating equation (see the "full
#' propagation" tests in test-integration.R).
#' @noRd
lca_vcov_distal_multinomial <- function(
  theta_hat,
  three_step.score,
  pi_adj,
  w.is,
  p.wx_mat,
  Y_cat,
  C,
  H.3.inv,
  Sigma.1,
  s2,
  iT,
  use.simple.cov,
  use.bch,
  Sigma.3 = NULL,
  s3.par = NULL,
  p.xz.cov = NULL,
  Z_mat_cov = NULL
) {
  J.3 <- three_step.score(theta_hat) # N x (iT*C)
  meat <- crossprod(J.3)

  if (use.bch || use.simple.cov) {
    return(H.3.inv %*% meat %*% H.3.inv)
  }

  pi_hat <- matrix(theta_hat, nrow = iT, ncol = C)
  n_par <- iT * C
  rec <- distal_records(
    log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300)),
    pi_adj,
    p.wx_mat
  )
  cross <- distal_cross_derivs(
    rec,
    w.is,
    p.wx_mat,
    distal_unit_derivs(theta_hat, Y_cat, iT, "multinomial", C = C)$G,
    Z_mat = if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
      Z_mat_cov
    } else {
      NULL
    }
  )

  step1.uncertainty <- cross$C1 %*% s2$J.2 %*% Sigma.1 %*% t(s2$J.2) %*% t(cross$C1)
  step2.uncertainty <- if (!is.null(cross$C_mat)) {
    cross$C_mat %*% Sigma.3 %*% t(cross$C_mat)
  } else {
    matrix(0, n_par, n_par)
  }

  H.3.inv %*% (meat + step1.uncertainty + step2.uncertainty) %*% H.3.inv
}


#' Three-step LCA estimation with covariates and/or distal outcomes
#'
#' Fits a three-step latent class model through the following steps:
#' \enumerate{
#'   \item \strong{Measurement model}: estimates latent class parameters
#'     (\eqn{\pi}, \eqn{\phi}) using \pkg{multilevLCA}
#'     (Lyrvall et al., 2025).
#'   \item \strong{Classification-error matrix}: computes posterior class
#'     probabilities and the T x T misclassification probability matrix
#'     \eqn{P(W = s \mid X = t)}, with standard errors corrected for
#'     classification-error propagation (Bakk, Oberski & Vermunt, 2014).
#'   \item \strong{Structural model}: estimates covariate effects using
#'     two-step starting values (Bakk & Kuha, 2018) and/or distal outcome
#'     means following Bakk, Tekle & Vermunt (2013), with the ML correction (Vermunt, 2010) or BCH correction
#'     (Bolck, Croon & Hagenaars, 2004). See
#'     \code{vignette("tseLCA", package = "tseLCA")} for a worked example.
#' }
#'
#' @param data A data.frame containing all columns referenced by \code{Y.names},
#'   \code{Zp.names}, and \code{Zo.name}.
#' @param Y.names Character vector of indicator column names. Need to be coded as consecutive integers with base level starting at `0`.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or \code{NULL}
#'   for a measurement-only fit. Default \code{NULL}.
#' @param Zo.name Single character name of the distal outcome column, or
#'   \code{NULL}. Default \code{NULL}.
#' @param step1 Pre-fitted Step-1 object (output of [tseLCA::lca_step1()] or a
#'   prior \code{three_step()} call), or \code{NULL} to run Step 1 internally.
#'   Default \code{NULL}.
#' @param startval Optional starting classification for the Step-1
#'   measurement model, either an integer vector of length \code{nrow(data)}
#'   (a class assignment \code{1..n_classes} for every row) or a numeric
#'   matrix of conditional item-response probabilities
#'   \eqn{P(Y_h = k \mid X = t)} (one row per item-category pair in
#'   \code{Y.names} order, one column per class) from which a classification
#'   is derived internally. See [lca_step1_startval()] for the full
#'   description of both forms and typical sources (an external solver run
#'   with many random starts, or a published item-response table).
#'   \pkg{multilevLCA}'s default initialization (k-means on principal
#'   components) is deterministic given the data and can converge to a local
#'   optimum of the Step-1 log-likelihood; supplying \code{startval} bypasses
#'   it entirely (\code{kmea = FALSE} with the classification injected as
#'   multilevLCA's \code{startval}). Mutually exclusive with \code{step1} and
#'   \code{n_init}. Default \code{NULL}.
#' @param n_init Optional positive integer. If supplied, fits the Step-1
#'   measurement model \code{n_init} times from independent uniform-random
#'   classifications (\code{kmea = FALSE}, not multilevLCA's k-means-on-PCA
#'   path) and keeps the fit with the highest log-likelihood -- the
#'   unconditional multi-start analog of \code{n_init} in \pkg{StepMix} or
#'   \code{nrep} in \pkg{poLCA}. This is distinct from
#'   \code{iter.measurement}, which instead reruns multilevLCA's own k-means
#'   initialization, and only when the entropy R\eqn{^2} of the default fit
#'   is below \code{R2.threshold}; \code{n_init} restarts always run.
#'   Mutually exclusive with \code{step1} and \code{startval}. Default
#'   \code{NULL}.
#' @param use.two.step Logical. Initialize Step-3 from two-step estimates.
#'   Default \code{TRUE}.
#' @param use.modal.assignment Logical. Use modal (hard) class assignments in
#'   Step 2 and 3. \code{FALSE} uses soft posterior weights. Default \code{TRUE}.
#' @param include.intercept Logical. Prepend an intercept column to the
#'   covariate design matrix. Default \code{TRUE}.
#' @param use.simple.cov Logical. Skip the Step-1 measurement-uncertainty
#'   correction and return only the robust sandwich variance. Faster but
#'   underestimates standard errors when class separation is low. Default
#'   \code{FALSE}.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}.
#'   Default \code{FALSE}.
#' @param boundary.tol Scalar. Parameters within this tolerance of 0 or 1 are
#'   treated as fixed when computing the Step-1 variance matrix for numerical stability. Default
#'   \code{1e-2}.
#' @param maxIter.measurement Integer. Maximum EM iterations for Step 1.
#'   Default \code{5000L}.
#' @param measurement.tol Scalar. Convergence tolerance for the Step-1 EM
#'   algorithm. Default \code{1e-8}.
#' @param covariate.tol Scalar. Convergence tolerance for the Step-3
#'   Newton-Raphson or EM algorithm. Default \code{1e-6}.
#' @param iter.measurement Integer. Number of random restarts triggered when
#'   the Step-1 entropy R\eqn{^2} falls below \code{R2.threshold}. Default
#'   \code{10L}.
#' @param R2.threshold Scalar. Entropy R\eqn{^2} threshold below which Step-1
#'   random restarts are triggered. Default \code{0.70}.
#' @param use.bch Logical. Use BCH-corrected weights instead of the ML
#'   estimator in Step 3. May error if BCH weights induce a non-positive semi-definite Hessian in the third step (common in cases of low separation). Default \code{FALSE}.
#' @param em.maxIter Integer. Maximum EM iterations for the Step-3 covariate
#'   or distal outcome model. Default \code{200L}.
#' @param get.twostep.vcov Logical. If \code{TRUE}, obtain \pkg{multilevLCA}'s
#'   bias-corrected variance-covariance matrix for the two-step gamma estimates
#'   and store it in \code{$two_step_vcov}. If the \code{fitZ} object passed
#'   through \code{step1} already contains a \code{Varmat_cor} (from a prior
#'   [fitZ_from_multiLCA()] or plain \code{multiLCA} call), it is attached
#'   automatically even when \code{get.twostep.vcov = FALSE}. Default
#'   \code{FALSE}.
#' @param rebase Character (e.g. \code{"C1"}, \code{"C2"}) or integer
#'   specifying which latent class to use as the reference category in the
#'   multinomial logit. The measurement model is permuted so this class becomes
#'   column 1 before any structural estimation. Default \code{"C1"}.
#' @param family Character. Distal outcome family: one of \code{"gaussian"}
#'   (class means), \code{"poisson"} (log-rates), \code{"binomial"}
#'   (logits), or \code{"multinomial"} (a saturated model for a nominal
#'   categorical outcome with 2 or more categories -- \code{Zo.name} may be
#'   a factor, character, or integer column; categories are taken from
#'   \code{sort(unique(data[[Zo.name]]))} with \code{factor()}). For
#'   \code{"multinomial"}, \code{coef()} returns a \code{T x C} matrix of
#'   class-conditional category probabilities
#'   \eqn{\hat\pi_{tc} = P(Zo = c \mid X = t)} (rows sum to 1) instead of a
#'   length-\code{T} vector, and \code{vcov()} returns its
#'   \code{(T*C) x (T*C)} sandwich covariance (necessarily singular, since
#'   each class's row sums to 1 -- see \code{\link{omnibus_test}()} for a
#'   Wald test that accounts for this). Unlike \code{"binomial"}, whose
#'   \code{coef()}/\code{vcov()} are on the logit scale, \code{"multinomial"}
#'   reports \code{coef()}/\code{vcov()} directly on the probability scale,
#'   so \code{Std.Error} is directly interpretable without a delta-method
#'   back-transform -- but a symmetric interval
#'   \code{Estimate +/- 1.96*Std.Error} can fall outside \eqn{[0, 1]} for a
#'   probability near a boundary, the same well-known limitation as a naive
#'   Wald interval for any sample proportion. The \code{z.value}/\code{p.value}
#'   columns \code{summary()}/\code{print()} show for this family test each
#'   probability against 0, which is rarely the question of interest;
#'   \code{\link{omnibus_test}()} is the intended, boundary-safe test of
#'   whether the outcome's distribution differs across classes. Combining
#'   \code{family = "multinomial"} with both \code{Zp.names} and
#'   \code{Zo.name} fully
#'   propagates both Step-1 measurement and Step-3 covariate uncertainty
#'   under \code{use.simple.cov = FALSE}, the same as the other families.
#'   Default \code{"gaussian"}.
#' @param correct.spec Logical. Use the model-robust outer-product Hessian for
#'   Step-3 standard errors rather than the observed-data Hessian. Not appropriate
#'   when the Step-3 model may be misspecified. Default \code{FALSE}.
#' @param verbose Logical. Print convergence messages. Default \code{FALSE}.
#'
#' @return An S3 object of class \code{tseLCA}. The subclass depends on which
#'   models were estimated:
#'   \describe{
#'     \item{`tseLCA_measurement`}{Returned when neither \code{Zp.names} nor
#'       \code{Zo.name} is supplied. Contains the following elements:
#'       \describe{
#'         \item{`measurement_model`}{Step-1 output list from [tseLCA::lca_step1()].}
#'         \item{`llik`}{Final Step-1 log-likelihood.}
#'         \item{`AIC`, `BIC`}{Information criteria from the measurement model.}
#'         \item{`R2entr`}{Entropy R\eqn{^2} of the measurement model.}
#'         \item{`n_classes`}{Number of latent classes.}
#'         \item{`posteriors`}{N x T matrix of soft posterior class probabilities.}
#'         \item{`classifications`}{Length-N integer vector of modal class assignments.}
#'       }
#'     }
#'     \item{`tseLCA_covariate`}{Returned when \code{Zp.names} is supplied and
#'       \code{Zo.name} is \code{NULL}. Contains all elements of
#'       \code{tseLCA_measurement} plus:
#'       \describe{
#'         \item{`three_step`}{Q x (T-1) matrix of Step-3 gamma coefficients.}
#'         \item{`three_step_vcov`}{Q(T-1) x Q(T-1) variance-covariance matrix
#'           for \code{three_step}, with measurement-uncertainty correction
#'           unless \code{use.simple.cov = TRUE}.}
#'         \item{`two_step`}{Q x (T-1) matrix of two-step starting values, or
#'           \code{NULL} if \code{use.two.step = FALSE}.}
#'         \item{`two_step_vcov`}{\pkg{multilevLCA} bias-corrected vcov for the
#'           two-step estimates, or \code{NULL}.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`entropy.R2`}{Covariate-adjusted entropy R\eqn{^2}.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|Z_{p,i};\hat{\gamma}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\phi}} held fixed. By construction
#'           smaller than the equivalent one-step MLE likelihood.}
#'       }
#'     }
#'     \item{`tseLCA_distal`}{Returned when \code{Zo.name} is supplied and
#'       \code{Zp.names} is \code{NULL}. Contains:
#'       \describe{
#'         \item{`three_step`}{Named length-T vector of Step-3 distal outcome
#'           parameters (means, log-rates, or logits depending on
#'           \code{family}) -- or, for \code{family = "multinomial"}, a
#'           \code{T x C} matrix of class-conditional category probabilities
#'           (rows sum to 1).}
#'         \item{`three_step_vcov`}{T x T variance-covariance matrix for
#'           \code{three_step}, named \code{mu_C1} through \code{mu_CT} --
#'           or, for \code{family = "multinomial"}, a \code{(T*C) x (T*C)}
#'           (necessarily rank-deficient) matrix named \code{"C{t}:{category}"}.}
#'         \item{`three_step.llik`}{Step-3 distal log-likelihood
#'           \eqn{\log P(Z_o|X=t)} at converged estimates.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|\hat{\pi}) P(Z_{o,i}|X=t;\hat{\mu}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\pi}, \hat{\phi}} held fixed.
#'           By construction smaller than the equivalent one-step MLE likelihood.}
#'         \item{`AIC`}{Akaike information criterion based on \code{llik}.}
#'         \item{`BIC`}{Bayesian information criterion based on \code{llik},
#'           using the number of distal-complete observations.}
#'         \item{`family`}{Character. The distal outcome family used.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`posteriors`}{N x T soft posterior matrix.}
#'         \item{`classifications`}{Length-N modal class assignment vector.}
#'       }
#'     }
#'     \item{`tseLCA_both`}{Returned when both \code{Zp.names} and
#'       \code{Zo.name} are supplied. Contains:
#'       \describe{
#'         \item{`covariate`}{A \code{tseLCA_covariate}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{entropy.R2}.}
#'         \item{`distal`}{A \code{tseLCA_distal}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{three_step.llik}.}
#'         \item{`family`, `n_classes`, `estimator`}{Shared top-level fields.}
#'         \item{`posteriors`, `classifications`}{Shared N x T posterior
#'           matrix and length-N modal class vector.}
#'       }
#'     }
#'   }
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#'   between latent class membership and external variables using bias-adjusted
#'   three-step approaches. \emph{Sociological Methodology}, 43(1), 272--311.
#'   \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' Bakk, Z., Pohle, M. J., & Kuha, J. (2025). Bias-adjusted three-step
#'   estimation of structural models for latent classes. \emph{Multivariate
#'   Behavioral Research}. \doi{10.1080/00273171.2025.2473935}
#'
#' @seealso \code{vignette("tseLCA", package = "tseLCA")} for a full worked
#'   example; [tseLCA::lca_step1()] for standalone Step-1 estimation
#'   (including from an externally supplied starting classification, with its
#'   own `startval` argument); [fitZ_from_fit0()] and [fitZ_from_multiLCA()]
#'   for two-step covariate estimation.
#'
#' @examples
#' d <- generate_data(n = 200, separation = "high",
#'                    scenario = "covariate", seed = 1)
#'
#' # Measurement model only
#' fit_m <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#'
#' # ML three-step with simple SEs (fast)
#' fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#'
#' # Full measurement-uncertainty correction (see vignette for interpretation)
#' fit_cor <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.simple.cov = FALSE,
#'                       use.modal.assignment = FALSE)
#' summary(fit_cor)
#'
#' # BCH estimator
#' fit_bch <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.bch = TRUE,
#'                       use.simple.cov = TRUE)
#' summary(fit_bch)
#'
#' # Change reference class
#' fit_c2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                      Zp.names = "Zp", use.simple.cov = TRUE,
#'                      rebase = "C2")
#' summary(fit_c2)
#'
#' # Gaussian distal outcome
#' d2 <- generate_data(200, "high", "distal", seed = 2)
#' fit_dis <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zo", family = "gaussian",
#'                       use.simple.cov = TRUE)
#' summary(fit_dis)
#'
#' # Nominal categorical distal outcome (3+ categories): coef() returns a
#' # T x C matrix of class-conditional category probabilities; omnibus_test()
#' # gives a single Wald test of whether the category distribution differs
#' # across classes at all.
#' d2$Zcat <- factor(sample(c("low", "mid", "high"), nrow(d2), replace = TRUE))
#' fit_cat <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zcat", family = "multinomial",
#'                       use.simple.cov = TRUE)
#' coef(fit_cat)
#' omnibus_test(fit_cat)
#'
#' # Pass a pre-fitted measurement model to skip Step 1
#' fit_step1 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' fit2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                    Zp.names = "Zp", step1 = fit_step1,
#'                    use.simple.cov = TRUE)
#' summary(fit2)
#'
#' # Supply an external starting classification for Step 1 (bypasses
#' # multilevLCA's k-means-on-PCA initialization; here we use the DGP's own
#' # true classes as a stand-in for e.g. a StepMix solution with many
#' # random starts)
#' fit_ext <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       startval = d$X, use.simple.cov = TRUE)
#' summary(fit_ext)
#'
#' # Many random-classification restarts for Step 1, keeping the best
#' # (analogous to n_init in StepMix or nrep in poLCA)
#' fit_ninit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                         n_init = 20L, use.simple.cov = TRUE)
#' summary(fit_ninit)
#'
#' # Plot item-response profiles from the measurement model
#' plot(fit)
#'
#' @export
three_step <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  Zo.name = NULL,
  step1 = NULL,
  startval = NULL,
  n_init = NULL,
  use.two.step = TRUE,
  use.modal.assignment = TRUE,
  include.intercept = TRUE,
  use.simple.cov = FALSE,
  incomplete = FALSE,
  boundary.tol = 1e-2,
  maxIter.measurement = 5000,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.bch = FALSE,
  em.maxIter = 200L,
  get.twostep.vcov = FALSE,
  rebase = "C1",
  family = "gaussian",
  correct.spec = FALSE,
  verbose = FALSE
) {
  # -- Step 1: measurement model -----------------------------------------------
  ref_idx <- parse_rebase(rebase, n_classes)

  n_step1_inputs <- sum(!is.null(step1), !is.null(startval), !is.null(n_init))
  if (n_step1_inputs > 1L) {
    stop(
      "`step1`, `startval`, and `n_init` are mutually exclusive ways of ",
      "controlling Step 1: supply a pre-fitted measurement model with ",
      "`step1`, a starting classification (or item-response probability ",
      "matrix) to fit one with `startval`, or a number of random restarts ",
      "with `n_init`, not more than one of these.",
      call. = FALSE
    )
  }

  if (!is.null(step1)) {
    # Normalize: accept raw lca_step1() list or any tseLCA object
    s1 <- if (inherits(step1, "tseLCA")) step1$measurement_model else step1
    # Apply rebase permutation so the desired reference class is column 1.
    s1$fit0 <- permute_fit0_classes(s1$fit0, ref_idx)
    if (!is.null(s1$fitZ)) {
      s1$fitZ <- normalize_fitZ_names(
        s1$fitZ,
        n_classes = n_classes
      )
      s1$fitZ <- permute_fitZ_classes(s1$fitZ, ref_idx)
    }
  } else {
    s1 <- lca_step1(
      data,
      Y.names,
      n_classes,
      Zp.names,
      maxIter.measurement,
      measurement.tol,
      covariate.tol,
      iter.measurement,
      R2.threshold,
      get.twostep.vcov,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      startval = startval,
      n_init = n_init,
      verbose = verbose
    )
  }
  fit0 <- s1$fit0

  # Compute two-step coefficients with EM algorithm holding step1 fixed.
  # This runs when use.two.step = TRUE and no fitZ was already computed
  # (e.g. when step1 was passed in from a measurement-only fit).
  if (use.two.step && is.null(s1$fitZ) && !is.null(Zp.names)) {
    s1$fitZ <- fitZ_from_fit0(
      fit0 = s1$fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = covariate.tol,
      maxIter = em.maxIter,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  }
  fitZ <- s1$fitZ

  # -- Multinomial distal outcome: encode categories as 1..C -------------------
  zo_levels <- NULL
  if (!is.null(Zo.name) && family == "multinomial") {
    zo_factor <- factor(data[[Zo.name]])
    zo_levels <- levels(zo_factor)
    if (length(zo_levels) < 2L) {
      stop(
        "`Zo.name` must have at least 2 distinct categories for ",
        "family = \"multinomial\".",
        call. = FALSE
      )
    }
    data[[Zo.name]] <- as.integer(zo_factor)
  }

  # -- Data preparation --------------------------------------------------------
  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    Zo.name = Zo.name,
    incomplete = incomplete,
    include.intercept = include.intercept,
    verbose = verbose
  )
  Y.obs <- cd$Y.obs
  mDesign <- cd$mDesign
  ivItemcat <- cd$ivItemcat
  keep_Y <- cd$keep_Y
  Z_mat <- cd$Z_mat
  keep_step3_Z_in_Y <- cd$keep_step3_Z_in_Y
  Zo_mat <- cd$Zo_mat
  keep_step3_Zo_in_Y <- cd$keep_step3_Zo_in_Y
  keep_step3_Zo <- cd$keep_step3_Zo

  s1$Y.names <- Y.names
  s1$ivItemcat <- ivItemcat
  s1$ref_idx <- ref_idx
  # Keep the Step-1 sample, in data-row order, with the measurement model: it
  # is needed for posteriors and for the Step-1 variance when this model is
  # reused (possibly on another sample) through `step1`. multilevLCA's
  # fit0$mU is not used as the data source because it is sorted by response
  # pattern and its polytomous coding differs between the listwise and FIML
  # paths.
  if (is.null(step1)) {
    s1$Y.exp <- Y.obs
    s1$mDesign.exp <- mDesign
  }

  # -- Early return if no covariates -------------------------------------------
  if (is.null(Zp.names) && is.null(Zo.name)) {
    posts <- step1_posteriors(Y.obs, mDesign, s1$fit0, ivItemcat)
    return(structure(
      list(
        measurement_model = s1,
        llik = s1$fit0$LLKSeries[nrow(s1$fit0$LLKSeries)],
        AIC = s1$fit0$AIC,
        BIC = s1$fit0$BIC,
        R2entr = s1$fit0$R2entr,
        n_classes = n_classes,
        npar = (n_classes - 1L) + n_classes * sum(ivItemcat - 1L),
        nobs = nrow(Y.obs),
        posteriors = posts,
        classifications = max.col(posts)
      ),
      class = c("tseLCA_measurement", "tseLCA")
    ))
  }

  # -- Step 2 ------------------------------------------------------------------
  s2 <- lca_step2(
    Y.obs,
    fit0,
    n_classes,
    use.modal.assignment,
    boundary.tol,
    use.simple.cov || use.bch,
    ivItemcat = ivItemcat,
    mDesign = mDesign
  )

  Q <- if (!is.null(Z_mat)) ncol(Z_mat) else 0L
  iT <- n_classes

  # Subset s2 outputs to the rows used in each Step 3 model.
  # s2 was estimated on all keep_Y rows; Step 3 models only use complete-Z rows.
  s2_for_cov <- if (!is.null(Z_mat)) {
    J.2_cov <- if (!is.null(s2$compute_J_unc)) {
      Y_cov <- Y.obs[keep_step3_Z_in_Y, , drop = FALSE]
      mDes_cov <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Z_in_Y, , drop = FALSE]
      } else {
        matrix(1L, nrow(Y_cov), ncol(Y_cov))
      }
      p_xy_cov <- s2$p.xy[keep_step3_Z_in_Y, , drop = FALSE]
      s2$compute_J_unc(
        p_xy_cov,
        Y_cov,
        mDes_cov,
        s2$theta1,
        ivItemcat,
        iT
      )
    } else {
      NULL
    }

    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = J.2_cov,
      w.is = s2$w.is[keep_step3_Z_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Z_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  s2_for_dis <- if (!is.null(Zo_mat)) {
    J.2_dis <- if (!is.null(s2$compute_J_unc)) {
      Y_dis <- Y.obs[keep_step3_Zo_in_Y, , drop = FALSE]
      mDes_dis <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        matrix(1L, nrow(Y_dis), ncol(Y_dis))
      }
      p_xy_dis <- s2$p.xy[keep_step3_Zo_in_Y, , drop = FALSE]
      s2$compute_J_unc(
        p_xy_dis,
        Y_dis,
        mDes_dis,
        s2$theta1,
        ivItemcat,
        iT
      )
    } else {
      NULL
    }

    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = J.2_dis,
      w.is = s2$w.is[keep_step3_Zo_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  # Step-1 sample inputs for Sigma.1.
  step1_Y <- step1_sample(s1, ivItemcat, ref_idx)
  if (is.null(step1_Y)) {
    step1_Y <- list(
      Y.exp = Y.obs,
      mDesign = mDesign,
      ivItemcat = ivItemcat,
      u_post = NULL
    )
  }

  #For covariate estimation
  if (!is.null(Z_mat)) {
    p.xz <- function(params) {
      eta_full <- cbind(0, Z_mat %*% params)
      row_max <- apply(eta_full, 1, max)
      exp_eta <- exp(eta_full - row_max)
      exp_eta / rowSums(exp_eta)
    }

    if (use.bch) {
      w.it <- bch_weight_matrix(s2_for_cov$w.is, s2_for_cov$p.wx_mat)

      .ll_bch <- function(params, pwx = NULL) {
        beta.cur <- matrix(params, ncol = iT - 1)
        rowSums(w.it * log(p.xz(beta.cur)))
      }

      neg.ll <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        -sum(w.it * log(pmax(p.xz(beta.cur), 1e-6)))
      }

      .grad_bch <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        -as.vector(t(Z_mat) %*% resid)
      }

      .score_bch <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        resid[, rep(seq_len(iT - 1L), each = Q)] *
          Z_mat[, rep(seq_len(Q), iT - 1L)]
      }

      three_step.ll <- .ll_bch
      three_step.grad <- .grad_bch
      three_step.score <- .score_bch
    } else {
      .ll_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        probs <- p.xz(matrix(params, ncol = iT - 1))
        rowSums(s2_for_cov$w.is * log(probs %*% t(pwx)))
      }

      .grad_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        beta <- matrix(params, ncol = iT - 1)
        p <- p.xz(beta)
        q <- p %*% t(pwx)
        r <- s2_for_cov$w.is / q
        grad <- matrix(0, nrow = iT - 1, ncol = Q)
        for (k in seq_len(iT - 1L)) {
          score_i <- p[, k + 1L] *
            (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
          grad[k, ] <- t(Z_mat) %*% score_i
        }
        as.vector(t(grad))
      }

      .score_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        beta <- matrix(params, ncol = iT - 1)
        p <- p.xz(beta)
        q <- p %*% t(pwx)
        r <- s2_for_cov$w.is / q
        score_ik <- matrix(0, nrow = nrow(Z_mat), ncol = iT - 1)
        for (k in seq_len(iT - 1L)) {
          score_ik[, k] <- p[, k + 1L] *
            (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
        }
        score_ik[, rep(seq_len(iT - 1L), each = Q)] *
          Z_mat[, rep(seq_len(Q), iT - 1L)]
      }

      three_step.ll <- .ll_ml
      three_step.grad <- .grad_ml
      three_step.score <- .score_ml
      neg.ll <- function(params) -sum(three_step.ll(params))
    }

    gamma_init <- if (!is.null(fitZ$mGamma) && use.two.step) {
      c(fitZ$mGamma)
    } else {
      rep(0, Q * (iT - 1))
    }

    # -- Step 3 ------------------------------------------------------------------
    s3 <- lca_step3(
      neg.ll,
      gamma_init,
      Q,
      iT,
      covariate.tol,
      gradient = three_step.grad,
      use.bch = use.bch,
      Z_mat_cc = Z_mat,
      w.is_cc = s2_for_cov$w.is,
      p.xz = p.xz,
      pwx = s2_for_cov$p.wx_mat,
      em.maxIter = em.maxIter,
      verbose = verbose,
      correct.spec = correct.spec
    )
    if (
      (correct.spec && !use.bch) ||
        is.null(s3$H.3.inv) ||
        any(is.na(s3$H.3.inv))
    ) {
      s3$H.3.inv <- qr.solve(crossprod(three_step.score(s3$res$par)))
    }

    coefs <- matrix(s3$res$par, ncol = iT - 1)
    ref_idx <- parse_rebase(rebase, iT)
    non_ref_classes <- seq_len(iT)[-ref_idx]
    colnames(coefs) <- paste0("C", non_ref_classes)
    rownames(coefs) <- c("Intercept", Zp.names)

    # -- Variance -----------------------------------------------------------------
    Sigma.3 <- lca_vcov(
      coefs = coefs,
      three_step.score = three_step.score,
      H.3.inv = s3$H.3.inv,
      Sigma.1 = if (use.simple.cov || use.bch) {
        NULL
      } else {
        lca_indiv_varmat(
          step1_Y$Y.exp,
          step1_Y$mDesign,
          fit0,
          step1_Y$ivItemcat,
          boundary.tol = boundary.tol,
          u_post = step1_Y$u_post
        )$Varmat
      },
      J.2 = s2_for_cov$J.2,
      p.wx_mat = s2_for_cov$p.wx_mat,
      w.is = s2_for_cov$w.is,
      Z_mat = Z_mat,
      n_classes = n_classes,
      p.xz = p.xz,
      s2 = s2_for_cov,
      use.simple.cov = use.simple.cov || use.bch
    )
    Sigma.3.covariate <- Sigma.3

    # -- Model fit ----------------------------------------------------------------
    Y_cc <- Y.obs[keep_step3_Z_in_Y, , drop = FALSE]
    mDes_cc <- if (!is.null(mDesign)) {
      mDesign[keep_step3_Z_in_Y, , drop = FALSE]
    } else {
      NULL
    }
    total.llik <- joint_log_lik(
      Y_cc,
      Z_mat,
      expand_Phi(fit0$mPhi, ivItemcat),
      coefs,
      mDes_cc
    )
    # Free parameters of the joint model: item-response log-ratios plus the
    # multinomial-logit coefficients (whose intercepts replace class sizes).
    total.k <- iT * sum(ivItemcat - 1L) + Q * (iT - 1L)
    total.AIC <- -2 * total.llik + 2 * total.k
    total.BIC <- -2 * total.llik + total.k * log(nrow(Y_cc))

    # -- Optional two-step vcov from multiLCA ------------------------------------
    # get.twostep.vcov = TRUE requests multilevLCA's bias-corrected SEs.
    # We skip re-estimation if a Varmat_cor is already attached to fitZ --
    # this handles three cases:
    #   (a) fitZ_from_multiLCA output: Varmat_cor lives at fitZ$raw_fit$Varmat_cor
    #   (b) Plain multiLCA output passed directly: fitZ$Varmat_cor
    #   (c) fitZ_from_fit0 output: no Varmat_cor anywhere -> must re-estimate
    #
    # If a Varmat_cor is already present on fitZ (regardless of get.twostep.vcov),
    # we always attach it to the output -- the user shouldn't lose it just because
    # they didn't set get.twostep.vcov = TRUE.

    .extract_varmat <- function(fZ) {
      if (is.null(fZ)) {
        return(NULL)
      }
      if (!is.null(fZ$Varmat_cor)) {
        return(fZ$Varmat_cor)
      }
      if (!is.null(fZ$raw_fit$Varmat_cor)) {
        return(fZ$raw_fit$Varmat_cor)
      }
      if (!is.null(fZ$raw_fit$SEs_cor_gamma)) {
        return(diag(as.vector(fZ$raw_fit$SEs_cor_gamma)^2))
      }
      NULL
    }

    .name_varmat <- function(V, fZ) {
      if (is.null(V) || is.null(fZ$mGamma)) {
        return(V)
      }
      param_names <- as.vector(outer(
        rownames(fZ$mGamma),
        colnames(fZ$mGamma),
        paste,
        sep = ":"
      ))
      rownames(V) <- param_names
      colnames(V) <- param_names
      V
    }

    existing_varmat <- .extract_varmat(fitZ)

    two_step_vcov <- if (!is.null(existing_varmat)) {
      # Already have it
      .name_varmat(existing_varmat, fitZ)
    } else if (get.twostep.vcov) {
      # No existing vcov
      fZ_ml <- fitZ_from_multiLCA(
        data = data,
        Y.names = Y.names,
        n_classes = n_classes,
        Zp.names = Zp.names,
        maxIter.measurement = maxIter.measurement,
        measurement.tol = measurement.tol,
        covariate.tol = covariate.tol,
        iter.measurement = iter.measurement,
        R2.threshold = R2.threshold,
        incomplete = incomplete,
        rebase = rebase,
        startval = startval,
        n_init = n_init,
        verbose = verbose
      )
      if (is.null(fitZ)) {
        fitZ <- fZ_ml
        s1$fitZ <- fZ_ml
      }
      raw_varmat <- .extract_varmat(fZ_ml)
      if (!is.null(raw_varmat)) {
        .name_varmat(raw_varmat, fZ_ml)
      } else {
        warning(
          "get.twostep.vcov: neither Varmat_cor nor SEs_cor_gamma found."
        )
        NULL
      }
    } else {
      NULL
    }

    # -- Covariate-adjusted entropy R^2 ------------------------------------------
    # Measures how much the items reduce classification uncertainty *beyond*
    # what the covariates already explain.
    #   error_prior = H(X|Z):   average entropy of P(X|Z_i) under fitted gamma
    #   error_post  = H(X|Y,Z): average entropy of P(X|Y_i,Z_i), recomputed
    #                            with the covariate-adjusted prior from
    #                            compute_pwx_adj (soft posterior assignment)
    #   R^2 = (H(X|Z) - H(X|Y,Z)) / H(X|Z)
    .h <- function(p) {
      p <- p[p > sqrt(.Machine$double.eps)]
      -sum(p * log(p))
    }
    pi_adj_cov <- p.xz(matrix(s3$res$par, ncol = iT - 1L)) # N x T

    error_prior <- mean(apply(pi_adj_cov, 1L, .h))

    adj_res <- compute_pwx_adj(
      Y.obs = Y_cc,
      fit0 = fit0,
      ivItemcat = ivItemcat,
      mDesign = mDes_cc,
      use.modal.assignment = FALSE,
      pi_adj = pi_adj_cov
    )
    error_post <- mean(apply(adj_res$post, 1L, .h))

    entropy.R2 <- if (error_prior > 1e-8) {
      (error_prior - error_post) / error_prior
    } else {
      1.0 # covariates already explain all class membership
    }

    s3.covariate <- structure(
      list(
        measurement_model = s1,
        two_step = if (!is.null(fitZ)) fitZ$mGamma else NULL,
        two_step_vcov = two_step_vcov,
        three_step = coefs,
        three_step_vcov = Sigma.3,
        three_step.llik = -s3$res$value,
        neg.ll = neg.ll,
        llik = total.llik,
        AIC = total.AIC,
        BIC = total.BIC,
        npar = total.k,
        nobs = nrow(Y_cc),
        n_classes = iT,
        estimator = if (use.bch) "BCH" else "ML",
        entropy.R2 = entropy.R2,
        posteriors = s2$p.xy,
        classifications = max.col(s2$p.xy)
      ),
      class = c("tseLCA_covariate", "tseLCA_structural", "tseLCA")
    )
  }

  if (!is.null(Zo_mat)) {
    if (!(family %in% c("gaussian", "poisson", "binomial", "multinomial"))) {
      message(
        'Provided family is not one of "gaussian", "poisson", "binomial", nor "multinomial". Defaulting to family="gaussain".'
      )
    }
    if (!is.null(Zp.names)) {
      Z_mat_dis <- if (!is.null(Z_mat) && length(keep_step3_Zo) > 0L) {
        Z_full_raw <- if (include.intercept) {
          m <- cbind(1, as.matrix(data[, Zp.names, drop = FALSE]))
          colnames(m) <- c("Intercept", Zp.names)
          m
        } else {
          as.matrix(data[, Zp.names, drop = FALSE])
        }
        Z_full_raw[keep_step3_Zo, , drop = FALSE]
      } else {
        NULL
      }

      if (!is.null(Z_mat_dis)) {
        pi_adj_full <- cbind(1, p.xz(matrix(s3$res$par, ncol = iT - 1)))
        p.xz_dis <- function(params) {
          eta_full <- cbind(0, Z_mat_dis %*% params)
          row_max <- apply(eta_full, 1L, max)
          exp_eta <- exp(eta_full - row_max)
          exp_eta / rowSums(exp_eta)
        }
        pi_adj <- p.xz_dis(matrix(s3$res$par, ncol = iT - 1))
      } else {
        pi_adj <- matrix(
          fit0$vPi,
          nrow = length(keep_step3_Zo_in_Y),
          ncol = iT,
          byrow = TRUE
        )
      }

      res_adj <- compute_pwx_adj(
        Y.obs[keep_step3_Zo_in_Y, , drop = FALSE],
        fit0,
        ivItemcat,
        if (!is.null(mDesign)) {
          mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
        } else {
          NULL
        },
        use.modal.assignment,
        pi_adj = pi_adj
      )
    } else {
      pi_adj <- matrix(
        fit0$vPi,
        nrow = length(keep_step3_Zo_in_Y),
        ncol = iT,
        byrow = TRUE
      )
      res_adj <- list(
        w.is = s2_for_dis$w.is,
        p.wx_mat = s2_for_dis$p.wx_mat
      )
    }

    w.is_dis <- res_adj$w.is

    if (family == "multinomial") {
      Y_cat_dis <- Zo_mat[, 1L] # already 1..C integer-coded, see above
      C <- length(zo_levels)

      s3.distal <- lca_step3.distal.multinomial(
        Y_cat = Y_cat_dis,
        C = C,
        iT = iT,
        covariate.tol = covariate.tol,
        use.bch = use.bch,
        w.is_cc = res_adj$w.is,
        pwx = res_adj$p.wx_mat,
        em.maxIter = em.maxIter,
        vPi = fit0$vPi,
        pi_mat = pi_adj,
        verbose = verbose
      )

      #Variance-covariance
      Sigma.3.distal <- lca_vcov_distal_multinomial(
        theta_hat = s3.distal$res$par,
        three_step.score = s3.distal$three_step.score,
        pi_adj = pi_adj,
        w.is = res_adj$w.is,
        p.wx_mat = res_adj$p.wx_mat,
        Y_cat = Y_cat_dis,
        C = C,
        H.3.inv = s3.distal$H.3.inv,
        Sigma.1 = if (use.simple.cov || use.bch) {
          NULL
        } else {
          lca_indiv_varmat(
            step1_Y$Y.exp,
            step1_Y$mDesign,
            fit0,
            step1_Y$ivItemcat,
            boundary.tol = boundary.tol,
            u_post = step1_Y$u_post
          )$Varmat
        },
        s2 = s2_for_dis,
        iT = iT,
        use.simple.cov = use.simple.cov,
        use.bch = use.bch,
        Sigma.3 = if (!is.null(Zp.names)) Sigma.3 else NULL,
        s3.par = if (!is.null(Zp.names)) s3$res$par else NULL,
        p.xz.cov = if (!is.null(Zp.names) && exists("p.xz_dis")) {
          p.xz_dis
        } else if (!is.null(Zp.names)) {
          p.xz
        } else {
          NULL
        },
        Z_mat_cov = if (!is.null(Zp.names) && !is.null(Z_mat_dis)) {
          Z_mat_dis
        } else if (!is.null(Zp.names)) {
          Z_mat
        } else {
          NULL
        }
      )

      class_labels <- paste0("C", seq_len(iT))
      pi_hat_mat <- matrix(s3.distal$res$par, nrow = iT, ncol = C)
      dimnames(pi_hat_mat) <- list(class_labels, zo_levels)
      param_labels <- as.vector(outer(
        class_labels,
        zo_levels,
        paste,
        sep = ":"
      ))
      dimnames(Sigma.3.distal) <- list(param_labels, param_labels)

      distal_par <- pi_hat_mat

      # -- Distal log-likelihood, AIC, BIC ------------------------------------
      distal.llik <- -s3.distal$res$value

      # log P(Zo_i = y_i | X=t) at converged pi_hat: N_dis x iT
      log_pZo_t <- log(pmax(t(pi_hat_mat[, Y_cat_dis, drop = FALSE]), 1e-300))

      Y_dis <- Y.obs[keep_step3_Zo_in_Y, , drop = FALSE]
      mDes_dis <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        NULL
      }

      total.llik.dis <- joint_log_lik_distal(
        Y = Y_dis,
        mPhi = expand_Phi(fit0$mPhi, ivItemcat),
        log_pZo_t = log_pZo_t,
        pi_mat = pi_adj, # N x T: P(X=t|Zp_i) or flat vPi
        mDesign = mDes_dis
      )

      # Class sizes, or the covariate logit coefficients when the class
      # prior depends on covariates (combined model), plus item parameters.
      n_meas_params <- (if (!is.null(Z_mat)) ncol(Z_mat) * (iT - 1L) else iT - 1L) +
        sum(ivItemcat - 1L) * iT
      n_distal_params <- iT * (C - 1L) # T x (C-1) free simplex parameters
      total.k.dis <- n_meas_params + n_distal_params
      N_dis <- length(keep_step3_Zo)

      s3.distal.list <- list(
        three_step = distal_par,
        three_step_vcov = Sigma.3.distal,
        three_step.llik = distal.llik,
        llik = total.llik.dis,
        AIC = -2 * total.llik.dis + 2 * total.k.dis,
        BIC = -2 * total.llik.dis + total.k.dis * log(N_dis),
        npar = total.k.dis,
        nobs = N_dis,
        zo_levels = zo_levels
      )
    } else {
      # Create p.zx : The function that maps latent indicators to oberved distal outcome Zo (need a choice of likelihood)

      if (family == "poisson") {
        p.zx <- function(params) {
          log_mu <- params[1:iT]
          mu <- exp(log_mu)
          z <- Zo_mat[, 1L]
          outer(z, log_mu, "*") - # N x T: z_i * log(mu_t)
            outer(rep(1, nrow(Zo_mat)), mu, "*") - # N x T: mu_t
            lgamma(z + 1L) # N x 1, recycled
        }
        starting.lm <- glm(
          Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
          family = poisson()
        )
        beta_init <- coef(starting.lm)
      } else if (family == "binomial") {
        # params: logit(mu_t)
        p.zx <- function(params) {
          logit_mu <- params[1:iT]
          mu <- 1 / (1 + exp(-logit_mu))
          z <- Zo_mat[, 1L]
          outer(z, log(mu), "*") + # N x T
            outer(1 - z, log(1 - mu), "*")
        }
        starting.lm <- glm(
          Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
          family = binomial()
        )
        beta_init <- coef(starting.lm) # already on logit scale
      } else {
        #(family == "gaussian")
        # params: class means, optionally followed by the common within-class
        # variance sigma2 (defaults to 1 only when omitted, e.g. for the
        # closed-form BCH means, which do not depend on it)
        p.zx <- function(params) {
          mu <- params[1:iT]
          sigma2 <- if (length(params) > iT) params[iT + 1L] else 1
          resid <- outer(Zo_mat[, 1L], mu, "-")
          -0.5 * resid^2 / sigma2 - 0.5 * log(2 * pi * sigma2)
        }
        starting.lm <- lm(Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)))
        beta_init <- coef(starting.lm)
      }

      if (use.bch) {
        w.it <- bch_weight_matrix(res_adj$w.is, res_adj$p.wx_mat)

        neg.ll <- function(params) {
          -sum(w.it * p.zx(params))
        }
      } else {
        # expanded-data log-likelihood, sum_i sum_s w_is log M_is
        # (R/distal-ml.R)
        neg.ll <- function(params) {
          rec <- distal_records(p.zx(params), pi_adj, res_adj$p.wx_mat)
          -distal_loglik(rec$logM, res_adj$w.is)
        }
      }

      s3.distal <- lca_step3.distal(
        neg.ll = neg.ll,
        em.maxIter = em.maxIter,
        pwx = res_adj$p.wx_mat,
        w.is_cc = res_adj$w.is,
        Zo_cc = Zo_mat[, 1L],
        use.bch = use.bch,
        covariate.tol = covariate.tol,
        iT = iT,
        beta_init = beta_init,
        family = family,
        p.zx = p.zx,
        vPi = fit0$vPi,
        pi_mat = pi_adj
      )

      # Full Step-3 parameter vector: class parameters, plus sigma2 for the
      # gaussian family (estimated jointly under ML; the BCH means do not
      # depend on it and it is estimated from the weighted residuals).
      gaussian_ml <- family == "gaussian" && !use.bch
      theta_zx <- if (family == "gaussian") {
        c(s3.distal$res$par, s3.distal$res$sigma2)
      } else {
        s3.distal$res$par
      }

      #Variance-covariance
      Sigma.3.distal <- lca_vcov_distal(
        mu_hat = if (gaussian_ml) theta_zx else s3.distal$res$par,
        three_step.score = s3.distal$three_step.score,
        pi_adj = pi_adj,
        w.is = res_adj$w.is,
        p.wx_mat = res_adj$p.wx_mat,
        p.zx = p.zx,
        family = family,
        H.3.inv = s3.distal$H.3.inv,
        Sigma.1 = if (use.simple.cov || use.bch) {
          NULL
        } else {
          lca_indiv_varmat(
            step1_Y$Y.exp,
            step1_Y$mDesign,
            fit0,
            step1_Y$ivItemcat,
            boundary.tol = boundary.tol,
            u_post = step1_Y$u_post
          )$Varmat
        },
        s2 = s2_for_dis,
        Sigma.3 = if (!is.null(Zp.names)) Sigma.3 else NULL,
        s3.par = if (!is.null(Zp.names)) s3$res$par else NULL,
        p.xz.cov = if (!is.null(Zp.names) && exists("p.xz_dis")) {
          p.xz_dis
        } else if (!is.null(Zp.names)) {
          p.xz
        } else {
          NULL
        },
        Z_mat_cov = if (!is.null(Zp.names) && !is.null(Z_mat_dis)) {
          Z_mat_dis
        } else if (!is.null(Zp.names)) {
          Z_mat
        } else {
          NULL
        },
        iT = iT,
        use.simple.cov = use.simple.cov,
        use.bch = use.bch,
        unit_scores = function(theta) {
          distal_unit_derivs(theta, Zo_mat[, 1L], iT, family)$G
        }
      )

      # Report the class parameters; keep sigma2 (and its SE under ML)
      # separately.
      sigma2_hat <- NULL
      if (family == "gaussian") {
        sigma2_hat <- c(
          estimate = s3.distal$res$sigma2,
          se = if (gaussian_ml) sqrt(Sigma.3.distal["sigma2", "sigma2"]) else NA_real_
        )
        Sigma.3.distal <- Sigma.3.distal[seq_len(iT), seq_len(iT), drop = FALSE]
      }

      distal_par <- s3.distal$res$par
      names(distal_par) <- paste0("mu_C", seq_len(iT))

      # -- Distal log-likelihood, AIC, BIC ----------------------------------------
      # Step-3 llik: log P(Zo|X=t) weighted by class assignments.
      # Total joint llik: sum_i log[ sum_t P(X=t|Zp_i) P(Zo_i|X=t) P(Y_i|X=t) ]
      distal.llik <- -neg.ll(theta_zx)

      # log P(Zo_i | X=t) at the converged parameters: N_dis x iT
      log_pZo_t <- p.zx(theta_zx)

      Y_dis <- Y.obs[keep_step3_Zo_in_Y, , drop = FALSE]
      mDes_dis <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        NULL
      }

      total.llik.dis <- joint_log_lik_distal(
        Y = Y_dis,
        mPhi = expand_Phi(fit0$mPhi, ivItemcat),
        log_pZo_t = log_pZo_t,
        pi_mat = pi_adj, # N x T: P(X=t|Zp_i) or flat vPi
        mDesign = mDes_dis
      )

      # Class sizes, or the covariate logit coefficients when the class
      # prior depends on covariates (combined model), plus item parameters.
      n_meas_params <- (if (!is.null(Z_mat)) ncol(Z_mat) * (iT - 1L) else iT - 1L) +
        sum(ivItemcat - 1L) * iT
      n_distal_params <- iT + as.integer(family == "gaussian") # + sigma2
      total.k.dis <- n_meas_params + n_distal_params
      N_dis <- length(keep_step3_Zo)

      s3.distal.list <- list(
        three_step = distal_par,
        three_step_vcov = Sigma.3.distal,
        three_step.llik = distal.llik,
        llik = total.llik.dis,
        AIC = -2 * total.llik.dis + 2 * total.k.dis,
        BIC = -2 * total.llik.dis + total.k.dis * log(N_dis),
        npar = total.k.dis,
        nobs = N_dis,
        sigma2 = sigma2_hat
      )
    } # end else (family != "multinomial")
  }

  if (!is.null(Zo_mat) && is.null(Z_mat)) {
    out <- s3.distal.list
    out$measurement_model <- s1
    out$family <- family
    out$n_classes <- iT
    out$estimator <- if (use.bch) "BCH" else "ML"
    out$posteriors <- s2$p.xy
    out$classifications <- max.col(s2$p.xy)
    class(out) <- c("tseLCA_distal", "tseLCA_structural", "tseLCA")
    return(out)
  }
  if (!is.null(Z_mat) && is.null(Zo_mat)) {
    class(s3.covariate) <- c("tseLCA_covariate", "tseLCA_structural", "tseLCA")
    return(s3.covariate)
  }

  out <- list(
    measurement_model = s1,
    covariate = s3.covariate,
    distal = s3.distal.list,
    family = family,
    n_classes = iT,
    estimator = if (use.bch) "BCH" else "ML",
    posteriors = s2$p.xy,
    classifications = max.col(s2$p.xy)
  )
  class(out) <- c("tseLCA_both", "tseLCA_structural", "tseLCA")
  return(out)
}
