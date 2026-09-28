# tseLCA/R/methods-tseLCA.R
#
# S3 class hierarchy and methods for fitted tseLCA objects.
#
#   tseLCA                      plot, logLik, nobs, posterior, classes
#   +- tseLCA_measurement       print, summary, coef, vcov, class_sizes, item_probs
#   +- tseLCA_structural        print, summary, coef, vcov  (confint and AIC/BIC
#      |                        work through the stats defaults)
#      +- tseLCA_covariate      Step-3 multinomial logit of class on covariates
#      +- tseLCA_distal         Step-3 class-specific distal outcome model
#      +- tseLCA_both           covariate model + distal outcome model
#
# Every fitted object carries `measurement_model` (the Step-1 fit), `llik`,
# `npar`, `nobs`, `n_classes`, `posteriors`, and `classifications`. A
# tseLCA_both object stores its two structural components in `$covariate`
# and `$distal`; its log-likelihood is that of the distal component, which is
# the full model for (Y, Zo) given Zp.

# -- internal helpers ------------------------------------------------------------

#' Fit-statistic carrier: the object itself, or the distal component of a
#' tseLCA_both object
#' @noRd
.fit_part <- function(x) if (inherits(x, "tseLCA_both")) x$distal else x

#' Name a coefficient matrix's entries "row:col" in column-major order, the
#' same convention used for the rows/columns of the Step-3 vcov matrices
#' @noRd
.flatten_coef <- function(m) {
  rn <- rownames(m)
  if (is.null(rn)) rn <- paste0("V", seq_len(nrow(m)))
  rn[rn == ""] <- "(Intercept)"
  cn <- colnames(m)
  if (is.null(cn)) cn <- paste0("C", seq_len(ncol(m)) + 1L)
  stats::setNames(as.vector(m), as.vector(outer(rn, cn, paste, sep = ":")))
}

#' Reject the tseLCA 1.x `which =` argument of coef()/vcov(), which would
#' otherwise be swallowed by `...` and silently return a different quantity
#' @noRd
.check_old_which <- function(...) {
  if ("which" %in% names(list(...))) {
    stop(
      "`which` was replaced in tseLCA 2.0: use `component = \"covariate\"` ",
      "or `\"distal\"` to select a component, and `step = \"two_step\"` for ",
      "two-step estimates.",
      call. = FALSE
    )
  }
}

#' Structural components of a fitted object, in print order
#' @noRd
.components <- function(x) {
  if (inherits(x, "tseLCA_both")) {
    list(covariate = x$covariate, distal = x$distal)
  } else if (inherits(x, "tseLCA_covariate")) {
    list(covariate = x)
  } else {
    list(distal = x)
  }
}

#' Coefficient vector of one structural component
#' @noRd
.component_coef <- function(part, type, step = "three_step") {
  est <- if (step == "two_step") part$two_step else part$three_step
  if (is.null(est)) {
    stop("No two-step estimates available in this object.", call. = FALSE)
  }
  if (type == "covariate") {
    return(.flatten_coef(est))
  }
  if (is.matrix(est)) {
    # multinomial distal outcome: T x C matrix, flattened column-major to
    # match the (t, c) naming of its vcov
    return(stats::setNames(as.vector(est), rownames(part$three_step_vcov)))
  }
  est
}

#' Variance matrix of one structural component
#' @noRd
.component_vcov <- function(part, step = "three_step") {
  if (step == "two_step") {
    if (is.null(part$two_step_vcov)) {
      stop(
        "No two-step vcov available. Set get.twostep.vcov = TRUE in three_step().",
        call. = FALSE
      )
    }
    return(part$two_step_vcov)
  }
  part$three_step_vcov
}

#' Estimate / Std. Error / z value / Pr(>|z|) matrix, as used by
#' stats::printCoefmat()
#' @noRd
.coefmat <- function(est, V) {
  se <- sqrt(diag(V))
  z <- est / se
  cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))
  )
}

.distal_scale <- function(family) {
  switch(
    family,
    gaussian = "means",
    poisson = "log means",
    binomial = "logits",
    multinomial = "category probabilities",
    "parameters"
  )
}

#' One-line reminder printed after a multinomial distal-outcome table
#'
#' Standard errors are on the probability scale, the per-cell z test is
#' against 0 (rarely of interest), and a symmetric interval can leave [0, 1];
#' omnibus_test() gives the intended test of whether the distribution differs
#' across classes.
#' @noRd
.multinomial_caveat_note <- function() {
  cat(
    "Note: standard errors are on the probability scale; the per-cell z test\n",
    "is against 0 (rarely of interest), and a symmetric interval can fall\n",
    "outside [0, 1] near a boundary. See omnibus_test() for a test of whether\n",
    "the distribution differs across classes.\n",
    sep = ""
  )
}

# -- accessor generics -------------------------------------------------------------

#' Posterior class-membership probabilities and modal class assignments
#'
#' `posterior()` returns the N x T matrix of posterior class-membership
#' probabilities used by a fitted model; `classes()` returns the modal
#' (most likely) class of each observation.
#'
#' @param object A fitted `tseLCA` object.
#' @param ... Further arguments (currently unused).
#' @return `posterior()`: a numeric N x T matrix. `classes()`: an integer
#'   vector of length N with values in `1..T`.
#' @examples
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' head(posterior(fit))
#' table(classes(fit))
#' @export
posterior <- function(object, ...) UseMethod("posterior")

#' @rdname posterior
#' @export
posterior.tseLCA <- function(object, ...) object$posteriors

#' @rdname posterior
#' @export
classes <- function(object, ...) UseMethod("classes")

#' @rdname posterior
#' @export
classes.tseLCA <- function(object, ...) object$classifications

#' Class sizes and item-response probabilities of the measurement model
#'
#' `class_sizes()` returns the estimated class proportions and
#' `item_probs()` the class-conditional item-response probabilities of the
#' Step-1 measurement model underlying any fitted `tseLCA` object.
#'
#' @param object A fitted `tseLCA` object.
#' @param ... Further arguments (currently unused).
#' @return `class_sizes()`: a named numeric vector of length T summing to one.
#'   `item_probs()`: a matrix with one row per item (binary items:
#'   \eqn{P(Y = 1 \mid X = t)}) or per item category (polytomous items:
#'   \eqn{P(Y = k \mid X = t)}) and one column per class.
#' @examples
#' d <- generate_data(200, "high", "covariate", seed = 1)
#' m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' class_sizes(m)
#' item_probs(m)
#' @export
class_sizes <- function(object, ...) UseMethod("class_sizes")

#' @rdname class_sizes
#' @export
class_sizes.tseLCA <- function(object, ...) {
  vPi <- as.vector(object$measurement_model$fit0$vPi)
  stats::setNames(vPi, paste0("C", seq_along(vPi)))
}

#' @rdname class_sizes
#' @export
item_probs <- function(object, ...) UseMethod("item_probs")

#' @rdname class_sizes
#' @export
item_probs.tseLCA <- function(object, ...) {
  object$measurement_model$fit0$mPhi
}

# -- methods shared by all tseLCA objects -------------------------------------------

#' Log-likelihood, number of observations, and information criteria
#'
#' `logLik()` returns the log-likelihood of a fitted model with its number of
#' free parameters (`df`) and observations (`nobs`), so that [stats::AIC()]
#' and [stats::BIC()] work directly. For a measurement model this is the
#' Step-1 log-likelihood. For structural models it is the log-likelihood of
#' the joint model for the indicators and the structural variables evaluated
#' with the Step-1 parameters held fixed; for a `tseLCA_both` object it is the
#' distal-outcome component, which conditions on the covariates.
#'
#' @param object A fitted `tseLCA` object.
#' @param ... Further arguments (currently unused).
#' @return `logLik()`: an object of class `"logLik"`. `nobs()`: an integer.
#' @examples
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' nobs(fit)
#' @export
logLik.tseLCA <- function(object, ...) {
  p <- .fit_part(object)
  structure(p$llik, df = p$npar, nobs = p$nobs, class = "logLik")
}

#' @rdname logLik.tseLCA
#' @export
nobs.tseLCA <- function(object, ...) .fit_part(object)$nobs

#' Plot item-response probability profiles for a tseLCA model
#'
#' Delegates to `plot.multiLCA` from \pkg{multilevLCA}, which draws the
#' class-specific item-response probability profiles of the Step-1
#' measurement model.
#'
#' @param x A fitted `tseLCA` object.
#' @param horiz Logical. If `TRUE`, item labels are drawn horizontally.
#' @param clab Optional character vector of length T giving class labels.
#' @param ... Further arguments passed to `plot.multiLCA`.
#' @return Called for its side effect (a base-graphics plot). Invisibly
#'   returns `NULL`.
#' @examples
#' d     <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' plot(fit_m)
#' plot(fit_m, clab = c("Low", "Mixed", "High"))
#' @export
plot.tseLCA <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(x$measurement_model$fit0, horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

# -- measurement model ---------------------------------------------------------------

#' Measurement-model parameters in the unconstrained log-ratio scale
#'
#' Ordering and names match `vcov.tseLCA_measurement()`: class-size
#' log-ratios `log(pi_t/pi_1)`, then, class by class and item by item,
#' `log(P(Y=k|C_t)/P(Y=0|C_t))` for each non-reference category k.
#' @noRd
.measurement_coef <- function(object) {
  fit0 <- object$measurement_model$fit0
  ivItemcat <- object$measurement_model$ivItemcat
  Y.names <- object$measurement_model$Y.names
  vPi <- as.vector(fit0$vPi)
  iT <- length(vPi)
  H <- length(ivItemcat)
  items <- if (!is.null(Y.names)) Y.names else paste0("Y", seq_len(H))
  cls <- paste0("C", seq_len(iT))
  Phi <- expand_Phi(fit0$mPhi, ivItemcat) # one row per (item, category)
  first <- cumsum(c(1L, utils::head(ivItemcat, -1L)))

  est <- log(vPi[-1L] / vPi[1L])
  nms <- paste0("log(pi_", cls[-1L], "/pi_", cls[1L], ")")
  for (t in seq_len(iT)) {
    for (h in seq_len(H)) {
      for (k in seq_len(ivItemcat[h] - 1L)) {
        est <- c(est, log(Phi[first[h] + k, t] / Phi[first[h], t]))
        nms <- c(nms, sprintf(
          "log(P(%s=%d|%s)/P(%s=0|%s))", items[h], k, cls[t], items[h], cls[t]
        ))
      }
    }
  }
  stats::setNames(est, nms)
}

#' Coefficients of a fitted tseLCA model
#'
#' Returns a named coefficient vector whose names match the rows and columns
#' of [vcov()], so that [stats::confint()] and other generic tools work.
#'
#' * Measurement models: class-size log-ratios \eqn{\log(\pi_t/\pi_1)} and
#'   item-response log-ratios \eqn{\log(P(Y=k \mid t)/P(Y=0 \mid t))}. Use
#'   [class_sizes()] and [item_probs()] for the probability scale.
#' * Covariate models: multinomial-logit coefficients named
#'   `"<covariate>:<class>"`, e.g. `"Zp:C2"`.
#' * Distal-outcome models: class-specific means (`gaussian`), log means
#'   (`poisson`), logits (`binomial`), or category probabilities
#'   (`multinomial`, named `"<class>:<category>"`).
#'
#' @param object A fitted `tseLCA` object.
#' @param component For `tseLCA_both` objects: `"all"` (default; covariate
#'   then distal coefficients), `"covariate"`, or `"distal"`.
#' @param step `"three_step"` (default) or `"two_step"` (the two-step
#'   estimates used to initialize Step 3; covariate models only).
#' @param matrix Logical. If `TRUE`, return the coefficients in their natural
#'   matrix layout instead: Q x (T-1) for covariate models, T x C for
#'   multinomial distal outcomes, a list of both for `tseLCA_both`.
#' @param ... Further arguments (currently unused).
#' @return A named numeric vector (or matrix / list if `matrix = TRUE`).
#' @examples
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' coef(fit)
#' coef(fit, matrix = TRUE)
#' confint(fit)
#' @export
coef.tseLCA_structural <- function(
  object,
  component = c("all", "covariate", "distal"),
  step = c("three_step", "two_step"),
  matrix = FALSE,
  ...
) {
  .check_old_which(...)
  component <- match.arg(component)
  step <- match.arg(step)
  parts <- .components(object)
  if (component != "all") {
    if (is.null(parts[[component]])) {
      stop(sprintf("This model has no %s component.", component), call. = FALSE)
    }
    parts <- parts[component]
  }
  if (step == "two_step" && is.null(parts$covariate)) {
    stop("Two-step estimates exist only for covariate models.", call. = FALSE)
  }
  if (step == "two_step") parts <- parts["covariate"]
  if (matrix) {
    out <- lapply(parts, function(p) if (step == "two_step") p$two_step else p$three_step)
    return(if (length(out) == 1L) out[[1L]] else out)
  }
  unlist(unname(Map(.component_coef, parts, names(parts), step)))
}

#' @rdname coef.tseLCA_structural
#' @export
coef.tseLCA_measurement <- function(object, ...) .measurement_coef(object)

#' Variance-covariance matrix of a fitted tseLCA model
#'
#' Row and column names match [coef()].
#'
#' * Measurement models: the BHHH variance matrix of the log-ratio
#'   parameters (attribute `"parameterization"` records the scale).
#' * Covariate and distal-outcome models: the Step-3 variance matrix, which
#'   includes the correction for Step-1 uncertainty unless the model was
#'   fitted with `use.simple.cov = TRUE`. For `family = "multinomial"` it is
#'   on the probability scale and rank-deficient (each class's probabilities
#'   sum to one).
#' * `tseLCA_both` with `component = "all"`: the covariate and distal blocks
#'   on the diagonal; the cross-covariances between the two sets of
#'   parameters are not computed and are `NA`.
#'
#' @inheritParams coef.tseLCA_structural
#' @param boundary.tol Measurement models only: parameters within this
#'   tolerance of 0 or 1 are treated as fixed. Default `1e-2`.
#' @return A named square matrix.
#' @examples
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' vcov(fit)
#' @export
vcov.tseLCA_structural <- function(
  object,
  component = c("all", "covariate", "distal"),
  step = c("three_step", "two_step"),
  ...
) {
  .check_old_which(...)
  component <- match.arg(component)
  step <- match.arg(step)
  parts <- .components(object)
  if (component != "all") {
    if (is.null(parts[[component]])) {
      stop(sprintf("This model has no %s component.", component), call. = FALSE)
    }
    parts <- parts[component]
  }
  if (step == "two_step") {
    if (is.null(parts$covariate)) {
      stop("Two-step estimates exist only for covariate models.", call. = FALSE)
    }
    parts <- parts["covariate"]
  }
  Vs <- lapply(parts, .component_vcov, step = step)
  if (length(Vs) == 1L) return(Vs[[1L]])
  nms <- unname(unlist(lapply(Vs, rownames)))
  V <- matrix(NA_real_, length(nms), length(nms), dimnames = list(nms, nms))
  at <- 0L
  for (Vi in Vs) {
    idx <- at + seq_len(nrow(Vi))
    V[idx, idx] <- Vi
    at <- at + nrow(Vi)
  }
  V
}

#' @rdname vcov.tseLCA_structural
#' @export
vcov.tseLCA_measurement <- function(object, boundary.tol = 1e-2, ...) {
  s1 <- object$measurement_model
  ref_idx <- if (!is.null(s1$ref_idx)) s1$ref_idx else 1L
  sample1 <- step1_sample(s1, s1$ivItemcat, ref_idx)
  if (is.null(sample1)) {
    stop(
      "The measurement model does not carry its Step-1 data; ",
      "re-estimate it to obtain its vcov.",
      call. = FALSE
    )
  }
  V <- lca_indiv_varmat(
    sample1$Y.exp, sample1$mDesign, s1$fit0, sample1$ivItemcat,
    boundary.tol = boundary.tol, u_post = sample1$u_post
  )$Varmat
  nms <- names(.measurement_coef(object))
  dimnames(V) <- list(nms, nms)
  attr(V, "parameterization") <- "log-ratio (unconstrained); NOT probabilities"
  V
}

# -- summary and print ------------------------------------------------------------

#' Summarize a fitted tseLCA model
#'
#' `summary()` collects fit statistics and coefficient tables; printing the
#' result formats the tables with [stats::printCoefmat()]. The coefficient
#' table (columns `Estimate`, `Std. Error`, `z value`, `Pr(>|z|)`) can be
#' extracted with `coef(summary(fit))`.
#'
#' @param object A fitted `tseLCA` object.
#' @param x A `summary.tseLCA_structural` or `summary.tseLCA_measurement`
#'   object, or a fitted `tseLCA` object (for `print`).
#' @param digits Number of significant digits to print.
#' @param signif.stars Logical; print significance stars?
#' @param ... Further arguments passed to [stats::printCoefmat()].
#' @return `summary()` returns an object of class `"summary.tseLCA_structural"`
#'   or `"summary.tseLCA_measurement"`. Print methods return their argument
#'   invisibly.
#' @examples
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' printCoefmat(coef(summary(fit)))
#' @export
summary.tseLCA_structural <- function(object, ...) {
  parts <- .components(object)
  tables <- Map(function(p, type) {
    .coefmat(.component_coef(p, type), p$three_step_vcov)
  }, parts, names(parts))
  fp <- .fit_part(object)
  structure(
    list(
      type = if (length(parts) == 2L) "both" else names(parts),
      n_classes = object$n_classes,
      estimator = object$estimator,
      family = object$family,
      logLik = logLik(object),
      AIC = fp$AIC,
      BIC = fp$BIC,
      entropy.R2 = if (!is.null(parts$covariate)) parts$covariate$entropy.R2,
      two_step = if (!is.null(parts$covariate)) parts$covariate$two_step,
      coefficients = tables
    ),
    class = "summary.tseLCA_structural"
  )
}

#' @rdname summary.tseLCA_structural
#' @export
coef.summary.tseLCA_structural <- function(object, ...) {
  do.call(rbind, unname(object$coefficients))
}

.structural_header <- function(x) {
  title <- switch(
    x$type,
    covariate = "Three-step latent class model: covariates",
    distal = "Three-step latent class model: distal outcome",
    both = "Three-step latent class model: covariates and distal outcome"
  )
  cat(title, "\n", sep = "")
  info <- sprintf("  Classes: %d   Estimator: %s", x$n_classes, x$estimator)
  if (x$type != "covariate") info <- paste0(info, "   Family: ", x$family)
  cat(info, "\n", sep = "")
  cat(sprintf(
    "  Log-lik: %.4f (df = %d)   AIC: %.2f   BIC: %.2f   N: %d\n",
    as.numeric(x$logLik), attr(x$logLik, "df"), x$AIC, x$BIC,
    attr(x$logLik, "nobs")
  ))
}

.print_tables <- function(x, digits, signif.stars, ...) {
  for (nm in names(x$coefficients)) {
    cat(if (nm == "covariate") {
      "\nCovariate effects on class membership (multinomial logit):\n"
    } else {
      sprintf("\nDistal outcome %s by class:\n", .distal_scale(x$family))
    })
    stats::printCoefmat(
      x$coefficients[[nm]], digits = digits, signif.stars = signif.stars,
      has.Pvalue = TRUE, P.values = TRUE, ...
    )
    if (nm == "distal" && identical(x$family, "multinomial")) {
      .multinomial_caveat_note()
    }
  }
}

#' @rdname summary.tseLCA_structural
#' @export
print.summary.tseLCA_structural <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  .structural_header(x)
  if (!is.null(x$entropy.R2)) {
    cat(sprintf("  Entropy R\u00b2 (covariate-adjusted): %.4f\n", x$entropy.R2))
  }
  if (!is.null(x$two_step)) {
    cat("\nTwo-step estimates (used to initialize Step 3):\n")
    print(x$two_step, digits = digits)
  }
  .print_tables(x, digits, signif.stars, ...)
  invisible(x)
}

#' @rdname summary.tseLCA_structural
#' @export
print.tseLCA_structural <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  s <- summary(x)
  .structural_header(s)
  .print_tables(s, digits, signif.stars, ...)
  invisible(x)
}

#' @rdname summary.tseLCA_structural
#' @export
summary.tseLCA_measurement <- function(object, ...) {
  structure(
    list(
      n_classes = object$n_classes,
      logLik = logLik(object),
      AIC = object$AIC,
      BIC = object$BIC,
      R2entr = object$R2entr,
      class_sizes = class_sizes(object),
      item_probs = item_probs(object)
    ),
    class = "summary.tseLCA_measurement"
  )
}

.measurement_header <- function(x) {
  cat("Latent class measurement model\n")
  cat(sprintf(
    "  Classes: %d   Log-lik: %.4f (df = %d)   AIC: %.2f   BIC: %.2f   N: %d\n",
    x$n_classes, as.numeric(x$logLik), attr(x$logLik, "df"), x$AIC, x$BIC,
    attr(x$logLik, "nobs")
  ))
  if (!is.null(x$R2entr)) cat(sprintf("  Entropy R\u00b2: %.4f\n", x$R2entr))
}

#' @rdname summary.tseLCA_structural
#' @export
print.summary.tseLCA_measurement <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
) {
  .measurement_header(x)
  cat("\nClass sizes:\n")
  print(x$class_sizes, digits = digits)
  cat("\nItem-response probabilities:\n")
  print(x$item_probs, digits = digits)
  invisible(x)
}

#' @rdname summary.tseLCA_structural
#' @export
print.tseLCA_measurement <- function(x, ...) {
  .measurement_header(summary(x))
  invisible(x)
}

# -- omnibus test --------------------------------------------------------------------

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

#' Build a between-class contrast matrix and run a generalized Wald test
#'
#' Tests \eqn{H_0: \theta_1 = \theta_2 = \dots = \theta_T}, where
#' \eqn{\theta_t} is the length-\code{d} distal-outcome parameter vector for
#' class \code{t} (\code{d = 1} for gaussian/poisson/binomial, \code{d = C}
#' -- the number of categories -- for multinomial), against class 1 as the
#' reference. \code{theta} and \code{V} must both use the column-major
#' \code{(t, k)} ordering \code{three_step()} produces (\code{as.vector()}
#' of a T x C matrix, or a plain length-T vector when \code{d = 1}): class
#' index varying fastest within blocks of \code{d}. A Moore-Penrose
#' pseudo-inverse (`pinv_rank()`) is used because the contrast covariance is
#' singular for multinomial (see `pinv_rank()`'s docs); its rank is used as
#' the chi-squared degrees of freedom, which for multinomial recovers the
#' textbook \eqn{(T-1)(C-1)} df of a chi-squared test of homogeneity in a
#' \eqn{T \times C} contingency table.
#' @noRd
wald_class_equality <- function(theta, V, iT) {
  d <- length(theta) %/% iT
  R <- matrix(0, (iT - 1L) * d, iT * d)
  for (k in seq_len(d)) {
    base_col <- (k - 1L) * iT
    base_row <- (k - 1L) * (iT - 1L)
    for (t in seq_len(iT - 1L)) {
      R[base_row + t, base_col + 1L] <- -1
      R[base_row + t, base_col + 1L + t] <- 1
    }
  }

  Rtheta <- R %*% theta
  RVR <- R %*% V %*% t(R)
  pr <- pinv_rank(RVR)

  statistic <- as.numeric(t(Rtheta) %*% pr$pinv %*% Rtheta)
  df <- pr$rank
  list(
    statistic = statistic,
    df = df,
    p.value = pchisq(statistic, df, lower.tail = FALSE)
  )
}

#' Omnibus Wald test of class equality for a distal outcome
#'
#' Tests whether the distal outcome distribution differs across latent
#' classes at all -- \eqn{H_0: \theta_1 = \theta_2 = \dots = \theta_T} for
#' the class-specific distal parameters \eqn{\theta_t} (class means for
#' \code{family = "gaussian"}, log-rates for \code{"poisson"}, logits for
#' \code{"binomial"}, or the full length-C category-probability vector for
#' \code{"multinomial"}) -- using a generalized Wald test with
#' \code{vcov()}. This answers whether the outcome's distribution is
#' associated with class membership at all, before drilling into which
#' classes differ. The degrees of freedom equal the rank of the contrast
#' covariance (`T - 1` for scalar outcomes; `(T - 1) * (C - 1)` for
#' multinomial, i.e. the textbook chi-squared test of homogeneity in a
#' \eqn{T \times C} table), computed with a Moore-Penrose pseudo-inverse so
#' the test remains valid despite \code{multinomial}'s inherently singular
#' covariance (each class's category probabilities sum to 1).
#'
#' @param object A \code{tseLCA_distal} object, or a \code{tseLCA_both}
#'   object (tests its distal component).
#' @param ... Unused; present for S3 method consistency.
#'
#' @return An object of class \code{"tseLCA_omnibus"}: a list with
#'   \code{$statistic} (the Wald chi-squared statistic), \code{$df}, and
#'   \code{$p.value}.
#' @examples
#' \donttest{
#' d <- generate_data(300, "high", "distal", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' omnibus_test(fit)
#' }
#' @export
omnibus_test <- function(object, ...) {
  UseMethod("omnibus_test")
}

#' @rdname omnibus_test
#' @export
omnibus_test.tseLCA_distal <- function(object, ...) {
  .omnibus_test_distal(object, object$n_classes, object$family)
}

#' @rdname omnibus_test
#' @export
omnibus_test.tseLCA_both <- function(object, ...) {
  .omnibus_test_distal(object$distal, object$n_classes, object$family)
}

#' @noRd
.omnibus_test_distal <- function(x, n_classes, family) {
  theta <- as.vector(x$three_step)
  V <- x$three_step_vcov
  iT <- n_classes

  test <- wald_class_equality(theta, V, iT)
  result <- list(
    statistic = test$statistic,
    df = test$df,
    p.value = test$p.value,
    family = if (!is.null(family)) family else "gaussian",
    n_classes = iT
  )
  class(result) <- "tseLCA_omnibus"
  result
}

#' @rdname omnibus_test
#' @param x A \code{tseLCA_omnibus} object.
#' @param digits Integer. Number of decimal places for the test statistic.
#' @export
print.tseLCA_omnibus <- function(x, digits = 4, ...) {
  cat("Omnibus Wald test of class equality (distal outcome)\n")
  cat(sprintf("  Family: %s   Classes: %d\n", x$family, x$n_classes))
  p_str <- if (x$p.value < 0.001) {
    "< 0.001"
  } else {
    sprintf("%.4f", x$p.value)
  }
  cat(sprintf(
    "  W(%d) = %.*f, p %s\n",
    x$df,
    digits,
    x$statistic,
    if (startsWith(p_str, "<")) p_str else paste0("= ", p_str)
  ))
  invisible(x)
}
