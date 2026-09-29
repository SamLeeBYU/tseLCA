#####################################################
### Simulation study for the tseLCA package
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee (samlee@arizona.edu)
#####################################################
#
# Reproduces the simulation study (Section 2.4) and its tables: bias, and
# coverage with SE/SD ratios, for the covariate and distal-outcome scenarios,
# with the three-step BCH and ML estimators under modal and proportional
# assignment, and the two-step estimator for covariates.
#
# Run in batch mode from a directory of your choice:
#   Rscript --vanilla tseLCA_sim.R
# or: source(system.file("examples", "tseLCA_sim.R", package = "tseLCA"))
#
# Settings (environment variables; defaults in parentheses):
#   TSELCA_SIM_QUICK      TRUE runs 25 replications per condition instead of
#                         500, for a quick check of the pipeline (FALSE).
#   TSELCA_SIM_SCENARIOS  "covariate", "distal", or both ("covariate,distal").
#   TSELCA_SIM_DIR        Directory for the simulated data, Step-1 models, and
#                         results ("tseLCA_sim_output").
#   TSELCA_SIM_CORES      Parallel workers (number of cores - 1).
#
# The computational cost is dominated by Step 1 (500 replications x 18
# conditions, each measurement model fitted from 20 random starts), which
# takes several hours on one core. The data and Step-1 models are saved in
# TSELCA_SIM_DIR and reused when the script is run again; runs are resumed
# condition by condition. Results for each replication are saved as well, and
# the tables are printed from them at the end.
#
# Monte Carlo error: with 500 replications the standard error of a coverage
# rate near .95 is about .01, and that of the bias is reported next to it.

library(tseLCA)
library(parallel)

QUICK <- isTRUE(as.logical(Sys.getenv("TSELCA_SIM_QUICK", "FALSE")))
SCENARIOS <- strsplit(Sys.getenv("TSELCA_SIM_SCENARIOS", "covariate,distal"), ",")[[1]]
OUT_DIR <- Sys.getenv("TSELCA_SIM_DIR", "tseLCA_sim_output")
N_CORES <- as.integer(Sys.getenv("TSELCA_SIM_CORES", max(1L, detectCores() - 1L)))

N_REP <- if (QUICK) 25L else 500L
SEP_LEVELS <- c("low", "mid", "high")
SAMPLE_SIZES <- c("500", "1000", "2000")
F_ITEMS <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
TRUTH <- c(covariate = 1, distal = 0) # Zp:C3 slope; mu_C3
TARGET <- c(covariate = "Zp:C3", distal = "mu_C3")
ALPHA <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf(
  "tseLCA %s simulation: %s, %d replications per condition, %d worker(s), output in %s\n",
  as.character(packageVersion("tseLCA")), paste(SCENARIOS, collapse = " + "),
  N_REP, N_CORES, normalizePath(OUT_DIR)
))

###################################################
### 1. Simulated data
###################################################

# 500 replications of each condition (scenario x separation x n) of the
# Bakk and Kuha (2018) design; QUICK mode uses the first N_REP of them.
data_path <- file.path(OUT_DIR, "sim_datasets.rds")
if (file.exists(data_path)) {
  datasets <- readRDS(data_path)
} else {
  datasets <- generate_all_conditions(
    n_rep = 500L, base_seed = 06262026L, sep_levels = SEP_LEVELS, verbose = FALSE
  )
  saveRDS(datasets, data_path)
}

conditions <- expand.grid(
  scenario = SCENARIOS, separation = SEP_LEVELS, n = SAMPLE_SIZES,
  stringsAsFactors = FALSE
)

###################################################
### 2. Step 1: measurement models
###################################################

# Each measurement model is fitted from 20 random starts (keeping the best).
# Only its parameters are stored; later steps rebuild it with as_tse_lca().
# For the covariate scenario, the two-step estimates and their corrected
# variance are obtained from multilevLCA (Lyrvall et al. 2025), initialized at
# the classes implied by this measurement model and the two-step coefficients.
#
# Models saved by tseLCA 1.x (lists with fit0 = list(vPi, mPhi) and fitZ) are
# read as well.

model_params <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.null(x$class_sizes)) return(x)
  list(
    class_sizes = as.vector(x$fit0$vPi),
    item_probs = x$fit0$mPhi,
    two_step = if (!is.null(x$fitZ)) list(coef = x$fitZ$mGamma, vcov = x$fitZ$Varmat_cor),
    converged = !isFALSE(x$fitZ_converged)
  )
}

fit_step1 <- function(d, sc, seed) {
  set.seed(seed)
  m <- tse_lca(F_ITEMS, data = d, nclass = 3, control = tse_control(n_init = 20))
  out <- list(class_sizes = class_sizes(m), item_probs = item_probs(m), converged = TRUE)
  if (sc == "covariate") {
    # classes implied by the measurement model and the two-step coefficients
    prior <- predict(tse_twostep(m, ~ Zp), newdata = d)
    w <- prior * posterior(m) / matrix(class_sizes(m), nrow(d), 3, byrow = TRUE)
    d$startval <- max.col(w)
    fz <- multilevLCA::multiLCA(
      data = d, Y = paste0("Y", 1:6), iT = 3L, Z = "Zp",
      startval = if (length(unique(d$startval)) == 3L) "startval" else NULL,
      extout = TRUE, verbose = FALSE
    )
    out$two_step <- list(coef = fz$mGamma, vcov = fz$Varmat_cor)
    out$converged <- abs(diff(utils::tail(fz$LLKSeries, 2))) < 1e-8
  }
  out
}

model_path <- file.path(OUT_DIR, "measurement_models.rds")
models <- if (file.exists(model_path)) readRDS(model_path) else list()

todo <- which(vapply(seq_len(nrow(conditions)), function(i) {
  cond <- conditions[i, ]
  have <- models[[cond$scenario]][[cond$separation]][[cond$n]]
  length(have) < N_REP || any(vapply(have[seq_len(N_REP)], is.null, logical(1)))
}, logical(1)))

if (length(todo) > 0L) {
  cat("Fitting Step-1 models for", length(todo), "condition(s)...\n")
  cl <- makePSOCKcluster(min(N_CORES, length(todo)))
  invisible(clusterEvalQ(cl, library(tseLCA)))
  clusterExport(cl, c("F_ITEMS", "fit_step1", "N_REP", "SEP_LEVELS", "datasets", "models"))
  fitted <- parLapply(cl, split(conditions[todo, ], seq_along(todo)), function(cond) {
    sc <- cond$scenario; sep <- cond$separation; nn <- cond$n
    have <- models[[sc]][[sep]][[nn]]
    if (is.null(have)) have <- vector("list", N_REP)
    for (r in seq_len(N_REP)) {
      if (length(have) >= r && !is.null(have[[r]])) next
      seed <- match(sc, c("covariate", "distal")) * 1e6 + match(sep, SEP_LEVELS) * 1e4 +
        as.integer(nn) + r
      have[[r]] <- tryCatch(fit_step1(datasets[[sc]][[sep]][[nn]][[r]], sc, seed),
                            error = function(e) NULL)
    }
    have
  })
  stopCluster(cl)
  for (k in seq_along(todo)) {
    cond <- conditions[todo[k], ]
    models[[cond$scenario]][[cond$separation]][[cond$n]] <- fitted[[k]]
  }
  saveRDS(models, model_path)
}

###################################################
### 3. Step 3: structural models, replication by replication
###################################################

run_replication <- function(d, p, sc) {
  if (is.null(p) || !isTRUE(p$converged)) return(NULL)
  m <- as_tse_lca(F_ITEMS, data = d, class_sizes = p$class_sizes, item_probs = p$item_probs)
  ctl <- tse_control(step3.maxit = 500)
  rows <- list()
  for (a in c("modal", "proportional")) {
    clf <- tse_classify(m, assignment = a)
    for (method in c("ML", "BCH")) {
      fit <- tryCatch(
        if (sc == "covariate") {
          # start from the two-step estimates, as in the original study
          start <- if (!is.null(p$two_step)) unname(p$two_step$coef) else NULL
          tse_covariate(clf, ~ Zp, method = method, start = start, control = ctl)
        } else {
          tse_distal(clf, Zo ~ 1, method = method, control = ctl)
        },
        error = function(e) NULL
      )
      est <- se <- NA_real_
      if (!is.null(fit)) {
        est <- coef(fit)[[TARGET[[sc]]]]
        se <- sqrt(diag(vcov(fit)))[[TARGET[[sc]]]]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        estimator = paste0(if (a == "modal") "modal." else "prop.", tolower(method)),
        estimate = est, se = se
      )
    }
  }
  if (sc == "covariate" && !is.null(p$two_step)) {
    rows[[length(rows) + 1L]] <- data.frame(
      estimator = "two_step",
      estimate = as.vector(p$two_step$coef)[4L],
      se = if (!is.null(p$two_step$vcov)) sqrt(diag(p$two_step$vcov))[4L] else NA_real_
    )
  }
  do.call(rbind, rows)
}

rep_path <- file.path(OUT_DIR, sprintf("sim_replicates%s.rds", if (QUICK) "_quick" else ""))
cat("Estimating Step-3 models...\n")
cl <- makePSOCKcluster(min(N_CORES, nrow(conditions)))
invisible(clusterEvalQ(cl, library(tseLCA)))
clusterExport(cl, c("F_ITEMS", "TARGET", "N_REP", "run_replication", "datasets", "models",
                    "model_params"))
replicates <- parLapply(cl, split(conditions, seq_len(nrow(conditions))), function(cond) {
  sc <- cond$scenario; sep <- cond$separation; nn <- cond$n
  out <- lapply(seq_len(N_REP), function(r) {
    res <- tryCatch(
      run_replication(datasets[[sc]][[sep]][[nn]][[r]],
                      model_params(models[[sc]][[sep]][[nn]][[r]]), sc),
      error = function(e) NULL
    )
    if (!is.null(res)) cbind(scenario = sc, separation = sep, n = nn, rep = r, res)
  })
  do.call(rbind, out)
})
stopCluster(cl)
replicates <- do.call(rbind, replicates)
rownames(replicates) <- NULL
saveRDS(replicates, rep_path)

###################################################
### 4. Tables
###################################################

# Performance of each estimator in each condition: bias and its Monte Carlo
# standard error, root mean squared error, coverage of 95% Wald intervals and
# its Monte Carlo standard error, and the ratio of the mean estimated standard
# error to the standard deviation of the estimates.
summarize_simulation <- function(replicates, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  groups <- split(replicates, replicates[c("scenario", "separation", "n", "estimator")],
                  drop = TRUE)
  rows <- lapply(groups, function(g) {
    ok <- is.finite(g$estimate)
    e <- g$estimate[ok]
    se <- g$se[ok]
    err <- e - TRUTH[[g$scenario[1]]]
    covered <- abs(err) <= z * se
    cov_rate <- mean(covered, na.rm = TRUE)
    n_cov <- sum(!is.na(covered))
    data.frame(
      scenario = g$scenario[1], separation = g$separation[1], n = as.integer(g$n[1]),
      estimator = g$estimator[1], n_ok = sum(ok),
      bias = mean(err), bias_mcse = sd(e) / sqrt(sum(ok)),
      rmse = sqrt(mean(err^2)),
      coverage = cov_rate, coverage_mcse = sqrt(cov_rate * (1 - cov_rate) / n_cov),
      se_sd = mean(se, na.rm = TRUE) / sd(e)
    )
  })
  out <- do.call(rbind, rows)
  out$separation <- factor(out$separation, SEP_LEVELS)
  out$estimator <- factor(out$estimator,
                          c("two_step", "modal.bch", "prop.bch", "modal.ml", "prop.ml"))
  out <- out[order(out$scenario, out$separation, out$n, out$estimator), ]
  rownames(out) <- NULL
  out
}

# Tables in the layout of the manuscript: one row per separation and n, one
# column (or column pair) per estimator.
num <- function(x, digits = 3) sub("^(-?)0\\.", "\\1.", formatC(x, format = "f", digits = digits))
print_tables <- function(s) {
  labels <- c(two_step = "2-step", modal.bch = "BCH modal", prop.bch = "BCH prop.",
              modal.ml = "ML modal", prop.ml = "ML prop.")
  for (sc in unique(s$scenario)) {
    ss <- s[s$scenario == sc, ]
    ests <- levels(droplevels(ss$estimator))
    key <- unique(ss[c("separation", "n")])
    cell <- function(k, est, f) {
      r <- ss[ss$separation == key$separation[k] & ss$n == key$n[k] & ss$estimator == est, ]
      if (nrow(r) == 0L) "" else f(r)
    }
    bias_tab <- sapply(ests, function(est) vapply(seq_len(nrow(key)), cell, "", est = est,
      f = function(r) sprintf("%s (%s)", num(r$bias), num(r$bias_mcse))))
    cov_tab <- sapply(ests, function(est) vapply(seq_len(nrow(key)), cell, "", est = est,
      f = function(r) sprintf("%s (%s) %s", num(r$coverage, 2), num(r$coverage_mcse, 2),
                              num(r$se_sd, 2))))
    rn <- sprintf("%-4s %4d", key$separation, key$n)
    dimnames(bias_tab) <- dimnames(cov_tab) <- list(rn, labels[ests])
    cat(sprintf("\n== %s scenario: bias (MCSE); target %s = %g ==\n", sc, TARGET[[sc]], TRUTH[[sc]]))
    print(noquote(bias_tab))
    cat(sprintf("\n== %s scenario: coverage (MCSE) SE/SD ==\n", sc))
    print(noquote(cov_tab))
  }
}

summary_tab <- summarize_simulation(replicates, ALPHA)
utils::write.csv(summary_tab, file.path(OUT_DIR, sprintf("sim_summary%s.csv", if (QUICK) "_quick" else "")),
                 row.names = FALSE)
print_tables(summary_tab)
cat(sprintf("\nReplications with estimates (of %d): %s\n", N_REP,
            paste(range(summary_tab$n_ok), collapse = "-")))
