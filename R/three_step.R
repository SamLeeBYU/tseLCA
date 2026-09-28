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
