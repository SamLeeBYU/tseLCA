# tseLCA/R/step1.R
#
# Step 1 (measurement model): internal helpers. The multilevLCA-based
# estimation routines themselves live in R/lca_measurement.R.

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
