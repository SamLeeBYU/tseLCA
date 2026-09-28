# tseLCA/R/step3-covariate.R
#
# Step 3 for explanatory covariates: bias-adjusted multinomial logit of
# latent class on covariates (BCH or ML).

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
