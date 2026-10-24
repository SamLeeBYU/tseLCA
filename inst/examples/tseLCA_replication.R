### R code from vignette source 'tseLCA.Rnw'

###################################################
### code chunk number 1: setup
###################################################
## Replication code for "tseLCA: An R Package for Three-step Estimation of
## Latent Class Models" (extracted from tseLCA.Rnw with Stangle).
##
## Run in a directory containing gss7677_tolerance.csv (otherwise it is
## downloaded) and, for a quick run, tseLCA_replication_cache.rds:
##   Rscript tseLCA.R                          all models (about 30 minutes)
##   TSELCA_REP_QUICK=TRUE Rscript tseLCA.R    stored multi-start fits (1 min)
## The tables are written to ./tables. The simulation study is replicated by
## tseLCA_sim.R (TSELCA_SIM_TABLES_ONLY=TRUE: tables from stored results).

library(tseLCA)

rm(list = ls())
gc()
r_opts <- options(
  prompt = "R> ",
  continue = "+  ",
  width = 77,
  digits = 4,
  useFancyQuotes = FALSE,
  warn = 1
)


###################################################
### code chunk number 2: sim-numbers
###################################################
## Numbers quoted in the text, computed from the saved simulation results
sim_f <- function(x, d = 2) sub("^(-?)0\\.", "\\1.", formatC(x, digits = d, format = "f"))
sim_and <- function(x, d = 2) sprintf("%s and %s", sim_f(min(x), d), sim_f(max(x), d))
sim_dash <- function(x, d = 2) sprintf("%s--%s", sim_f(min(x), d), sim_f(max(x), d))
sim_list <- function(x, d = 2) {
  x <- sim_f(x, d)
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}
sim <- read.csv(file.path("tseLCA_sim_output", "sim_summary.csv"))
sim_reps <- readRDS(file.path("tseLCA_sim_output", "sim_replicates.rds"))
sim_cov <- subset(sim, scenario == "covariate")
sim_dis <- subset(sim, scenario == "distal")
sim_pick <- function(d, est = NULL, sep = NULL, n = NULL) {
  d <- d[order(d$n), ]
  d[(is.null(est) | d$estimator %in% est) & (is.null(sep) | d$separation %in% sep) &
      (is.null(n) | d$n %in% n), ]
}
## number of replications, and the number of conditions of one scenario
sim_nrep <- max(sim$n_ok)
sim_ncond <- nrow(sim_pick(sim_cov, "two_step"))
## covariate effects: bias, coverage, and SE/SD
cb_low <- abs(sim_pick(sim_cov, sep = "low", n = 500)$bias)
cb_mid <- max(abs(sim_pick(sim_cov, sep = "mid")$bias))
cc_ml <- sim_pick(sim_cov, "prop.ml", "low")$coverage
cc_2s <- sim_pick(sim_cov, "two_step", "low")$coverage
cs_bch_mid <- sim_pick(sim_cov, c("modal.bch", "prop.bch"), "mid")$se_sd
cs_modal_low <- sim_pick(sim_cov, "modal.ml", "low")$se_sd
cs_prop_low <- sim_pick(sim_cov, "prop.ml", "low")$se_sd
## non-convergence of the BCH estimator of the covariate effects
cf <- function(est) sim_nrep - sim_pick(sim_cov, est)$n_ok
cf_prop <- sum(cf("prop.bch")); cf_modal <- sum(cf("modal.bch"))
cf_ml <- sum(cf(c("modal.ml", "prop.ml", "two_step")))
## the ML estimates under low separation, and the dependence of the estimates
## and their standard errors (ML with proportional assignment, n = 500)
lowml <- subset(sim_reps, scenario == "covariate" & separation == "low" &
  estimator %in% c("modal.ml", "prop.ml") & is.finite(estimate))
dev_max <- max(abs(lowml$estimate - 1))
pml <- subset(lowml, estimator == "prop.ml" & n == 500)
pml_cor <- cor(pml$se, pml$estimate)
pml_q <- cut(pml$se, quantile(pml$se, 0:4 / 4), include.lowest = TRUE)
pml_cover1 <- tapply(abs(pml$estimate - 1) <= qnorm(.975) * pml$se, pml_q, mean)[[1]]
## distal outcome
db_range <- function(sep) max(tapply(sim_pick(sim_dis, sep = sep)$bias,
  sim_pick(sim_dis, sep = sep)$n, function(x) max(x) - min(x)))
db_low <- db_range("low"); db_mid <- db_range("mid")
db_mh500 <- max(abs(sim_pick(sim_dis, sep = c("mid", "high"), n = 500)$bias))
dc_high <- sim_pick(sim_dis, sep = "high")$coverage
dc_mid <- sim_pick(sim_dis, sep = "mid")$coverage
dse_mid <- sim_pick(sim_dis, sep = "mid")$se_sd
dc_bch_low <- sim_pick(sim_dis, c("modal.bch", "prop.bch"), "low")$coverage
dc_ml_low <- sim_pick(sim_dis, c("modal.ml", "prop.ml"), "low")$coverage


###################################################
### code chunk number 3: sim-tables
###################################################
## The simulation tables are written by tseLCA_sim.R; rebuild them from the
## saved replication results (tseLCA_sim_output/sim_replicates.rds) if needed.
if (!file.exists(file.path("tseLCA_sim_output", "tables", "sim_covariate_bias.tex"))) {
  system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "tseLCA_sim.R"),
    env = c("TSELCA_SIM_TABLES_ONLY=TRUE", "TSELCA_SIM_DIR=tseLCA_sim_output"))
}


###################################################
### code chunk number 4: illus-setup
###################################################
library("poLCA")
N_STARTS <- 50L # random starts per measurement model
## Sample sizes quoted in the text of the two examples
data("election", package = "poLCA")
N_ELECTION <- nrow(election)
gss_raw <- read.csv("gss7677_tolerance.csv")
gss_items0 <- c("atheists", "communists", "militarists", "racists", "homosexuals")
G_N <- nrow(gss_raw)
G_N76 <- sum(gss_raw$year == 1976)
G_N77 <- sum(gss_raw$year == 1977)
G_N5 <- sum(complete.cases(gss_raw[gss_items0]))
G_NCC <- sum(complete.cases(gss_raw[c(gss_items0, "cohort", "education")]))
G_TOL <- colMeans(gss_raw[complete.cases(gss_raw[gss_items0]), gss_items0] == 1)
## Fits from many random starts are stored in tseLCA_replication_cache.rds.
## With TSELCA_REP_QUICK=TRUE they are read from it (fits missing from it are
## computed and added); otherwise all models are estimated anew.
QUICK <- isTRUE(as.logical(Sys.getenv("TSELCA_REP_QUICK", "FALSE")))
cache_file <- "tseLCA_replication_cache.rds"
cache <- if (QUICK && file.exists(cache_file)) readRDS(cache_file) else list()
cache_changed <- FALSE
fit_or_load <- function(name, expr) {
  if (QUICK && !is.null(cache[[name]])) return(cache[[name]])
  cache[[name]] <<- expr
  cache_changed <<- TRUE
  cache[[name]]
}
## LaTeX tables in ./tables, formatted as the other tables
dir.create("tables", showWarnings = FALSE)
num <- function(x, d = 2) {
  x[round(x, d) == 0] <- 0
  sub("^(-?)0\\.", "\\1.", formatC(x, digits = d, format = "f"))
}
est_se <- function(est, se, d = 2) paste0(num(est, d), " (", num(se, d), ")")
## All tables are set at the common width \tablewidth (the natural width of
## the widest table of the paper) and scaled to the text width, so that they
## have the same font size, with equally wide gaps between the columns and
## margins left and right (\eqtab in the preamble); tables/<label>_rows.tex
## holds the rows and tables/<label>_measure.tex the call that measures them.
write_table <- function(body, header, align, caption, label, footer = NULL) {
  rows <- c(if (is.matrix(body)) {
    paste0(apply(body, 1, paste, collapse = " & "), " \\\\")
  } else body, footer)
  inner <- c("\\toprule", header, "\\midrule", rows, "\\bottomrule")
  ## the rows, and the call that measures the table for the common width
  rows_file <- file.path("tables", paste0(label, "_rows"))
  n_col <- length(strsplit(gsub("[^a-zA-Z]", "", gsub("\\{[^}]*\\}", "", align)), "")[[1]])
  writeLines(inner, paste0(rows_file, ".tex"))
  writeLines(sprintf("\\measureeq{%d}{%s}{\\inputrows{%s}}", n_col, align, rows_file),
             file.path("tables", paste0(label, "_measure.tex")))
  writeLines(c("\\begin{table}[t!]", "\\centering", "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\eqtab{%d}{%s}{\\inputrows{%s}}}", n_col, align, rows_file),
               sprintf("\\caption{\\footnotesize %s}", caption),
               sprintf("\\label{tab:%s}", label), "\\end{table}"),
             file.path("tables", paste0(label, ".tex")))
}
tex_line <- function(...) paste(paste(c(...), collapse = " & "), "\\\\")
## a column header broken before its last word, e.g. "Intolerant of / right"
wrap_hdr <- function(x) vapply(x, function(s) {
  w <- strsplit(s, " ")[[1]]
  if (length(w) < 2L) s else sprintf("\\shortstack{%s\\\\%s}",
    paste(w[-length(w)], collapse = " "), w[length(w)])
}, "", USE.NAMES = FALSE)
## a row with a little extra space below it, for rows whose label spans two lines
tex_line_sp <- function(...) paste(paste(c(...), collapse = " & "), "\\\\[3pt]")
prob_se <- function(model, row = NULL) {
  x <- if (is.null(row)) class_sizes(model, se = TRUE) else item_probs(model, se = TRUE)
  if (is.null(row)) rbind(est = x$estimate, se = x$se)
  else rbind(est = x$estimate[row, ], se = x$se[row, ])
}


###################################################
### code chunk number 5: elec-data
###################################################
data("election", package = "poLCA")
items <- colnames(election)[1:12]
elec <- election[complete.cases(election[, c(items, "PARTY")]), ]
f_elec <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ 1
nrow(elec)


###################################################
### code chunk number 6: elec-select (eval = FALSE)
###################################################
## set.seed(20260928)
## sel <- tse_lca(f_elec, data = elec, nclass = 1:6,
##   control = tse_control(n_init = 50))
## sel


###################################################
### code chunk number 7: elec-select-run
###################################################
set.seed(20260928)
sel <- fit_or_load("sel_elec", tse_lca(f_elec, data = elec, nclass = 1:6,
  control = tse_control(n_init = N_STARTS)))
sel


###################################################
### code chunk number 8: elec-onestep
###################################################
elec_pol <- elec
elec_pol[items] <- lapply(elec_pol[items], as.integer)
f_onestep <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ PARTY
onestep <- function(n_cl) {
  ## each start separately: one failed start stops poLCA's own nrep loop
  best <- NULL
  for (s in seq_len(if (n_cl == 1) 1L else N_STARTS)) {
    fit <- tryCatch(poLCA(f_onestep, elec_pol, nclass = n_cl, maxiter = 5000,
      verbose = FALSE), error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$llik) &&
        (is.null(best) || fit$llik > best$llik)) best <- fit
  }
  best
}
set.seed(20260928)
one <- fit_or_load("onestep_elec", lapply(1:6, onestep))


###################################################
### code chunk number 9: elec-measurement
###################################################
m_elec <- sel[[3]]
m_elec


###################################################
### code chunk number 10: elec-labels
###################################################
G <- c("MORALG", "CARESG", "KNOWG", "LEADG", "INTELG")
B <- c("MORALB", "CARESB", "KNOWB", "LEADB", "INTELB")
label_elec <- function(post) {
  ## Gore minus Bush mean rating of the positive traits (lower = better)
  codes <- sapply(elec[c(G, B)], as.integer)
  means <- crossprod(post, codes) / colSums(post)
  score <- rowMeans(means[, G]) - rowMeans(means[, B])
  lab <- rep("Other", ncol(post))
  lab[which.min(score)] <- "Gore affinity"
  lab[which.max(score)] <- "Bush affinity"
  lab
}
lab_elec <- label_elec(posterior(m_elec))
ord_elec <- match(c("Gore affinity", "Bush affinity", "Other"), lab_elec)
traits <- c(MORAL = "Moral", CARES = "Caring", KNOW = "Knowledgeable",
  LEAD = "Good leader", DISHON = "Dishonest", INTEL = "Intelligent")
## one-step (poLCA) classes, labeled by the same rule
one3 <- one[[3]]
lab_one <- label_elec(one3$posterior)
ord_one <- match(c("Gore affinity", "Bush affinity", "Other"), lab_one)
cells <- function(three, one_est, one_se) {
  as.vector(rbind(est_se(three["est", ord_elec], three["se", ord_elec]),
    est_se(one_est[ord_one], one_se[ord_one])))
}
body <- rbind(c("", "Class size", cells(prob_se(m_elec), one3$P, one3$P.se)))
for (cand in c("Gore", "Bush")) for (tr in names(traits)) {
  item <- paste0(tr, substr(cand, 1, 1))
  body <- rbind(body, c(if (tr == "MORAL") cand else "", traits[[tr]],
    cells(prob_se(m_elec, sprintf("P(%s.0|C)", item)),
      one3$probs[[item]][, 1], one3$probs.se[[item]][, 1])))
}
cls <- c("Gore affinity", "Bush affinity", "Other")
## entropy R^2 from the posterior class probabilities, relative to the entropy
## of the class proportions (as reported by tseLCA)
entropy_r2 <- function(post) {
  if (ncol(post) < 2) return(NA_real_)
  h <- function(p) -sum(p[p > 1e-12] * log(p[p > 1e-12]))
  1 - mean(apply(post, 1, h)) / h(colMeans(post))
}
write_table(body,
  header = c(paste(" & &", paste(sprintf("\\multicolumn{2}{c}{%s}", cls),
      collapse = " & "), "\\\\"),
    "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
    tex_line("Candidate", "Trait", rep(c("Three-step", "One-step"), 3))),
  align = "ll cc cc cc",
  caption = paste("Three-class models for the ANES 2000 candidate evaluations:",
    "class sizes and probabilities that a trait describes the candidate",
    "``extremely well'', with standard errors in parentheses. Three-step: the",
    "measurement model; one-step: the latent class regression on \\code{PARTY}",
    "(\\pkg{poLCA}), whose class sizes are averages over respondents. The entropy",
    sprintf("$R^2$ is %s for the three-step and %s for the one-step model.",
      num(as.data.frame(sel)$entropy.R2[3]), num(entropy_r2(one3$posterior)))),
  label = "elec_measurement")


###################################################
### code chunk number 11: elec-classify
###################################################
cl_elec <- tse_classify(m_elec, assignment = "modal")
cl_elec


###################################################
### code chunk number 12: elec-covariate
###################################################
other <- which(lab_elec == "Other")
fc_elec <- tse_covariate(cl_elec, ~ PARTY, ref = other)
summary(fc_elec)


###################################################
### code chunk number 13: elec-estimators
###################################################
fits_elec <- list(
  "ML, corrected SE" = fc_elec,
  "ML, robust SE" = update(fc_elec, se = "robust"),
  "BCH" = update(fc_elec, method = "BCH"),
  "two-step" = tse_twostep(m_elec, ~ PARTY, ref = other, se = TRUE),
  "uncorrected" = update(fc_elec, method = "none"))


###################################################
### code chunk number 14: elec-party-table
###################################################
party <- grep("^PARTY:", names(coef(fc_elec)), value = TRUE)
est <- sapply(fits_elec, function(f) coef(f)[party])
se <- sapply(fits_elec, function(f) sqrt(diag(vcov(f)))[party])
rownames(est) <- rownames(se) <- paste("PARTY:", lab_elec[-other])
## poLCA's three-class model: coefficients against its first class,
## re-expressed against its "Other" class
one3 <- one[[3]]
lab_one <- label_elec(one3$posterior)
rebase_polca <- function(fit, ref) {
  n_cl <- ncol(fit$coeff) + 1L
  Q <- nrow(fit$coeff)
  M <- matrix(0, n_cl - 1L, n_cl - 1L)
  for (i in seq_len(n_cl - 1L)) {
    t <- seq_len(n_cl)[-ref][i]
    if (t != 1L) M[i, t - 1L] <- 1
    if (ref != 1L) M[i, ref - 1L] <- -1
  }
  A <- kronecker(M, diag(Q))
  list(est = matrix(A %*% as.vector(fit$coeff), Q, n_cl - 1L,
    dimnames = list(rownames(fit$coeff), seq_len(n_cl)[-ref])),
    vcov = A %*% fit$coeff.V %*% t(A))
}
one_rb <- rebase_polca(one3, which(lab_one == "Other"))
one_se <- matrix(sqrt(diag(one_rb$vcov)), nrow(one_rb$est))
agree <- mean(lab_one[max.col(one3$posterior)] == lab_elec[classes(cl_elec)])
gb <- paste("PARTY:", c("Gore affinity", "Bush affinity"))
labels_est <- c("ML, corrected SE" = "Three-step ML, corrected SE",
  "ML, robust SE" = "Three-step ML, robust SE", "BCH" = "Three-step BCH",
  "two-step" = "Two-step", "uncorrected" = "Three-step, uncorrected")
## estimate (SE), 95% Wald interval, and the joint Wald test of both PARTY
## coefficients (2 df)
ci <- function(b, se) sprintf("[%s, %s]", num(b - qnorm(.975) * se, 3),
  num(b + qnorm(.975) * se, 3))
wald_party <- function(b, V) {
  w <- drop(t(b) %*% solve(V, b))
  p <- pchisq(w, length(b), lower.tail = FALSE)
  c(sprintf("%.1f", w), if (p < .001) "$<$.001" else num(p, 3))
}
party_row <- function(label, b, se, V) {
  c(label, as.vector(rbind(est_se(b, se, d = 3), ci(b, se))), wald_party(b, V))
}
body <- t(sapply(names(labels_est), function(e) {
  f <- fits_elec[[e]]
  nm <- sprintf("PARTY:C%d", match(c("Gore affinity", "Bush affinity"), lab_elec))
  party_row(labels_est[[e]], coef(f)[nm], sqrt(diag(vcov(f)))[nm], vcov(f)[nm, nm])
}))
## poLCA: PARTY is the second row of each class's coefficients
one_cols <- match(c("Gore affinity", "Bush affinity"), lab_one[lab_one != "Other"])
idx_party <- 2 * one_cols
ratio_bush <- one_rb$est["PARTY", one_cols[2]] /
  coef(fits_elec[["ML, corrected SE"]])[[sprintf("PARTY:C%d", which(lab_elec == "Bush affinity"))]]
body <- rbind(body, party_row("One-step (\\pkg{poLCA})",
  one_rb$est["PARTY", one_cols], one_se[2, one_cols],
  one_rb$vcov[idx_party, idx_party, drop = FALSE]))
write_table(unname(body),
  header = c(paste(" & \\multicolumn{2}{c}{Gore affinity} & \\multicolumn{2}{c}{Bush affinity}",
      "& \\multicolumn{2}{c}{Wald test} \\\\"),
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
    tex_line("Estimator", rep(c("Estimate (SE)", "95\\% CI"), 2), "$\\chi^2_2$", "$p$")),
  align = "l cc cc cc",
  caption = paste("Effect of party identification on membership of the Gore- and",
    "Bush-affinity classes relative to the ``Other'' class: multinomial logit",
    "coefficients with standard errors, 95\\% Wald confidence intervals (all",
    "coefficients have $p < .001$), and the Wald test that party identification",
    "has no effect on class membership (2 df). Three-class models of the ANES",
    paste0("2000 data ($N = ", nrow(elec), "$); the three-step estimators use modal assignment.")),
  label = "elec_party")


###################################################
### code chunk number 15: elec-plot
###################################################
party_x <- seq(1, 7, length.out = 101)
p_three <- predict(fc_elec, newdata = data.frame(PARTY = party_x))
colnames(p_three) <- lab_elec
eta_one <- cbind(0, cbind(1, party_x) %*% one3$coeff)
p_one <- exp(eta_one) / rowSums(exp(eta_one))
colnames(p_one) <- lab_one
p_one <- p_one[, lab_elec]
par(mar = c(4, 4, 1, 1))
matplot(party_x, p_three, type = "l", lwd = 3, col = 1, lty = 1, ylim = c(0, 1),
  xlab = "Party ID: strong Democrat (1) to strong Republican (7)",
  ylab = "Probability of latent class membership")
matlines(party_x, p_one, lwd = 2, col = 1, lty = 2)
peak <- apply(p_three, 2, which.max)
for (k in seq_along(peak)) {
  text(party_x[peak[k]], max(p_three[, k]) + 0.05, lab_elec[k],
    adj = c((peak[k] - 1) / (length(party_x) - 1), 0.5))
}
legend("topright", c("tseLCA (three-step)", "poLCA (one-step)"),
  lty = 1:2, lwd = c(3, 2), col = 1, bty = "n")


###################################################
### code chunk number 16: elec-onecall
###################################################
fit_elec <- tseLCA(cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ PARTY, data = elec,
  nclass = 3, ref = other, start = item_probs(m_elec))
all.equal(coef(fit_elec), coef(fc_elec), tolerance = 1e-4)


###################################################
### code chunk number 17: gss-download
###################################################
## gss7677_tolerance.csv is an extract of GSScombined_1977original.tab from
## the replication archive of Bakk, Oberski & Vermunt (2014), Harvard
## Dataverse, doi:10.7910/DVN/24497 (CC0); rebuilt from it when missing
if (!file.exists("gss7677_tolerance.csv")) {
  sav <- tempfile(fileext = ".sav")
  utils::download.file(
    "https://dataverse.harvard.edu/api/access/datafile/2456935?format=original",
    sav, mode = "wb", quiet = TRUE)
  raw <- suppressWarnings(foreign::read.spss(sav, to.data.frame = TRUE,
    use.value.labels = FALSE))
  write.csv(data.frame(year = raw$year, id = raw$id, atheists = raw$atheists_,
    communists = raw$communists_, homosexuals = raw$homtest,
    militarists = raw$militarists_, racists = raw$racist_, cohort = raw$cohort_,
    education = raw$education, polviews = raw$polviews, partyid = raw$partyid,
    natrace = raw$natrace), "gss7677_tolerance.csv", row.names = FALSE, na = "")
}


###################################################
### code chunk number 18: gss-data
###################################################
gss <- read.csv("gss7677_tolerance.csv")
gss$cohort <- factor(gss$cohort, 0:3,
  c("after 1951", "1934-1951", "1915-1933", "1914 or before"))
gss$education <- factor(gss$education, 0:2,
  c("< 12 years", "12 years", "> 12 years"))
gss$natrace <- factor(gss$natrace, 1:3,
  c("too little", "about right", "too much"))
f_gss <- cbind(atheists, communists, militarists, racists,
  homosexuals) ~ 1


###################################################
### code chunk number 19: gss-select (eval = FALSE)
###################################################
## set.seed(20260929)
## sel_gss <- tse_lca(f_gss, data = gss, nclass = 1:5,
##   control = tse_control(n_init = 50))
## m_gss <- sel_gss[[4]]


###################################################
### code chunk number 20: gss-select-run
###################################################
set.seed(20260929)
sel_gss <- fit_or_load("sel_gss", tse_lca(f_gss, data = gss, nclass = 1:5,
  control = tse_control(n_init = N_STARTS)))
m_gss <- sel_gss[[4]]


###################################################
### code chunk number 21: gss-labels
###################################################
ip <- item_probs(m_gss)
rownames(ip) <- c("atheists", "communists", "militarists", "racists", "homosexuals")
avg <- colMeans(ip)
mid <- setdiff(1:4, c(which.max(avg), which.min(avg)))
lab_gss <- character(4)
lab_gss[which.max(avg)] <- "Tolerant"
lab_gss[which.min(avg)] <- "Intolerant"
lab_gss[mid] <- ifelse(colMeans(ip[c("atheists", "communists"), mid]) >
  colMeans(ip[c("militarists", "racists"), mid]),
  "Intolerant of right", "Intolerant of left")
gss_classes <- c("Intolerant", "Tolerant", "Intolerant of right", "Intolerant of left")
## one-step comparison: both approaches on the respondents with complete
## covariates, for one to five classes
items_gss <- c("atheists", "communists", "militarists", "racists", "homosexuals")
gss_cc <- gss[complete.cases(gss[, c(items_gss, "cohort", "education")]), ]
set.seed(20260930)
sel_gss_cc <- fit_or_load("sel_gss_cc", tse_lca(f_gss, data = gss_cc,
  nclass = 1:5, control = tse_control(n_init = N_STARTS)))
gss_pol <- gss_cc
gss_pol[items_gss] <- lapply(gss_pol[items_gss], function(x) as.integer(x) + 1L)
f_gss1 <- cbind(atheists, communists, militarists, racists, homosexuals) ~
  cohort + education
onestep_gss <- function(n_cl) {
  best <- NULL
  for (s in seq_len(if (n_cl == 1) 1L else N_STARTS)) {
    fit <- tryCatch(poLCA(f_gss1, gss_pol, nclass = n_cl, maxiter = 5000,
      verbose = FALSE), error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$llik) &&
        (is.null(best) || fit$llik > best$llik)) best <- fit
  }
  best
}
set.seed(20260930)
one_gss <- fit_or_load("onestep_gss", lapply(1:5, onestep_gss))
## BIC tables of the two examples (smallest value of each column in bold)
bic_col <- function(x) {
  out <- formatC(x, format = "f", digits = 1)
  out[which.min(x)] <- sprintf("\\textbf{%s}", out[which.min(x)])
  out
}
r2_col <- function(x) ifelse(is.na(x), "--", sub("^0", "", formatC(x, format = "f", digits = 2)))
bic_table <- function(sel_fits, one_fits, covariates, data_name, n, label) {
  three <- as.data.frame(sel_fits)
  write_table(cbind(seq_along(three$BIC),
      bic_col(three$AIC), bic_col(three$BIC), r2_col(three$entropy.R2),
      bic_col(sapply(one_fits, function(f) f$aic)),
      bic_col(sapply(one_fits, function(f) f$bic)),
      r2_col(sapply(one_fits, function(f) entropy_r2(f$posterior)))),
    header = c(" & \\multicolumn{3}{c}{Three-step} & \\multicolumn{3}{c}{One-step} \\\\",
      "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
      tex_line("Classes ($T$)", rep(c("AIC", "BIC", "Entropy $R^2$"), 2))),
    align = "c ccc ccc",
    caption = paste0("AIC, BIC, and entropy $R^2$ of the three-step measurement model and of the ",
      "one-step latent class regression (\\pkg{poLCA}) on ", covariates, ", ", data_name,
      " ($N = ", n, "$), for $T$ classes; the smallest AIC and BIC of each approach are in bold. ",
      "Both approaches use the respondents with complete covariates. The posterior class ",
      "probabilities of the one-step model also use the covariates, so its entropy $R^2$ is ",
      "similar to or higher than that of the three-step model, which uses the items only."),
    label = label)
}
bic_table(sel, one, "\\code{PARTY}", "ANES 2000", nrow(elec), "bic_elec")
bic_table(sel_gss_cc, one_gss, "cohort and education", "GSS 1976--77", nrow(gss_cc),
  "bic_gss")
## the four-class one-step model: classes labeled from their profiles
p1 <- one_gss[[4]]
tol1 <- sapply(items_gss, function(it) p1$probs[[it]][, 2]) # P(tolerant)
avg1 <- rowMeans(tol1)
lab_one_gss <- character(4)
lab_one_gss[which.max(avg1)] <- "Tolerant"
lab_one_gss[which.min(avg1)] <- "Intolerant"
mid1 <- setdiff(1:4, c(which.max(avg1), which.min(avg1)))
big1 <- mid1[which.max(p1$P[mid1])]
small1 <- setdiff(mid1, big1)
lab_one_gss[big1] <- "Partially tolerant"
lab_one_gss[small1] <- paste("Tolerant of", names(which.max(tol1[small1, ])))
one_classes <- c("Intolerant", "Tolerant", "Partially tolerant", lab_one_gss[small1])
ord1 <- match(one_classes, lab_one_gss)
## the three-class one-step model (the number of classes selected by its BIC) and
## the three-step model with three classes, labeled by the same rule
label3 <- function(avg) {
  lab <- character(3)
  lab[which.max(avg)] <- "Tolerant"
  lab[which.min(avg)] <- "Intolerant"
  lab[setdiff(1:3, c(which.max(avg), which.min(avg)))] <- "Partially tolerant"
  lab
}
p3 <- one_gss[[3]]
lab_p3 <- label3(rowMeans(sapply(items_gss, function(it) p3$probs[[it]][, 2])))
classes3 <- c("Intolerant", "Tolerant", "Partially tolerant")
ord3 <- match(classes3, lab_p3)
m3_gss <- sel_gss_cc[[3]]
lab_m3 <- label3(colMeans(item_probs(m3_gss)))
agree3 <- mean(lab_p3[max.col(p3$posterior)] ==
  lab_m3[classes(tse_classify(m3_gss, assignment = "modal"))])
ord_gss <- match(gss_classes, lab_gss)
bakk_est <- rbind(size = c(.56, .23, .11, .10),
  atheists = c(.03, .98, .41, .61), communists = c(.04, .95, .59, .27),
  militarists = c(.05, .92, .34, .38), racists = c(.08, .90, .02, .81),
  homosexuals = c(.13, .96, .72, .56))
bakk_se <- rbind(size = c(.02, .01, .03, .03),
  atheists = c(.01, .01, .06, .07), communists = c(.01, .02, .11, .07),
  militarists = c(.01, .02, .05, .06), racists = c(.01, .02, .06, .20),
  homosexuals = c(.01, .01, .07, .06))
row_label <- function(row) {
  if (row == "size") "Class size" else sprintf("\\shortstack[l]{Tolerant of\\\\%s}", row)
}
cm4 <- "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}"
panel_a <- sprintf("\\multicolumn{9}{l}{\\textit{Three-step: measurement model ($N = %d$)}} \\\\", G_N5)
bakk_diff <- bakk_est
for (row in rownames(bakk_est)) {
  s <- if (row == "size") prob_se(m_gss) else prob_se(m_gss, sprintf("P(%s|C)", row))
  bakk_diff[row, ] <- s["est", ord_gss] - bakk_est[row, ]
  panel_a <- c(panel_a, tex_line_sp(row_label(row),
    as.vector(rbind(est_se(s["est", ord_gss], s["se", ord_gss]),
      est_se(bakk_est[row, ], bakk_se[row, ])))))
}
## largest difference from the estimates of Bakk et al., except for the probability
## that the intolerant of right tolerate racists (on the boundary in our estimates)
bakk_diff["racists", which(gss_classes == "Intolerant of right")] <- NA
bakk_max <- max(abs(bakk_diff), na.rm = TRUE)
panel_b <- c("\\midrule",
  paste("\\multicolumn{9}{l}{\\textit{One-step, four classes, estimated with cohort and",
    sprintf("education (\\pkg{poLCA}, $N = %d$)}} \\\\[3pt]", G_NCC)),
  tex_line("", sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(one_classes))), cm4)
for (row in rownames(bakk_est)) {
  est <- if (row == "size") p1$P else p1$probs[[row]][, 2]
  se <- if (row == "size") p1$P.se else p1$probs.se[[row]][, 2]
  panel_b <- c(panel_b, tex_line_sp(row_label(row),
    sprintf("\\multicolumn{2}{c}{%s}", est_se(est[ord1], se[ord1]))))
}
panel_c <- c("\\midrule",
  paste("\\multicolumn{9}{l}{\\textit{One-step, three classes",
    sprintf("estimated with cohort and education (\\pkg{poLCA}, $N = %d$)}} \\\\[3pt]", G_NCC)),
  tex_line("", sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(classes3))),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}")
for (row in rownames(bakk_est)) {
  est <- if (row == "size") p3$P else p3$probs[[row]][, 2]
  se <- if (row == "size") p3$P.se else p3$probs.se[[row]][, 2]
  panel_c <- c(panel_c, tex_line_sp(row_label(row),
    sprintf("\\multicolumn{2}{c}{%s}", est_se(est[ord3], se[ord3]))))
}
write_table(c(panel_a, panel_b, panel_c),
  header = c(paste(" &", paste(sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(gss_classes)),
    collapse = " & "), "\\\\"), cm4,
    paste(" &", paste(rep("\\pkg{tseLCA} & Bakk et al.", 4), collapse = " & "), "\\\\")),
  align = "l cc cc cc cc",
  caption = paste("Latent class models of tolerance for nonconformity (GSS 1976--77),",
    "with class sizes and probabilities of a tolerant answer and standard errors in",
    "parentheses. Top, the four-class three-step measurement model, estimated by",
    "\\pkg{tseLCA} and reported by \\citet[Table~7]{Bakk2014}. Middle and bottom, the",
    "one-step latent class regressions on cohort and education (\\pkg{poLCA}) with four",
    "and with three classes, whose class sizes are averages over respondents.",
    sprintf("The entropy $R^2$ is %s for the three-step and %s and %s for the one-step models.",
      num(as.data.frame(sel_gss)$entropy.R2[4]), num(entropy_r2(p1$posterior)),
      num(entropy_r2(p3$posterior)))),
  label = "gss_measurement")


###################################################
### code chunk number 22: gss-covariate
###################################################
cl_gss <- tse_classify(m_gss, assignment = "proportional")
intolerant <- which(lab_gss == "Intolerant")
fc_gss <- tse_covariate(cl_gss, ~ cohort + education, ref = intolerant)
anova(fc_gss)


###################################################
### code chunk number 23: gss-covariate-table
###################################################
unc_gss <- tse_covariate(tse_classify(m_gss, assignment = "modal"),
  ~ cohort + education, method = "none", ref = intolerant)
terms <- c("cohort1934-1951" = "Born 1934--1951", "cohort1915-1933" = "Born 1915--1933",
  "cohort1914 or before" = "\\shortstack[l]{Born 1914\\\\or before}",
  "education12 years" = "\\shortstack[l]{12 years of\\\\education}",
  "education> 12 years" = "\\shortstack[l]{More than 12\\\\years of education}")
se_c <- sqrt(diag(vcov(fc_gss)))
se_r <- sqrt(diag(vcov(update(fc_gss, se = "robust"))))
se_u <- sqrt(diag(vcov(unc_gss)))
panel_a <- "\\multicolumn{7}{l}{\\textit{Three-step (classes of the measurement model)}} \\\\"
ratio_cu <- numeric(0) # corrected relative to uncorrected coefficients
for (tm in names(terms)) {
  cells <- character(0)
  for (k in gss_classes[-1]) {
    nm <- sprintf("%s:C%d", tm, which(lab_gss == k))
    ratio_cu <- c(ratio_cu, coef(fc_gss)[[nm]] / coef(unc_gss)[[nm]])
    cells <- c(cells, est_se(coef(fc_gss)[nm], se_c[nm]),
      est_se(coef(unc_gss)[nm], se_u[nm]))
  }
  panel_a <- c(panel_a, tex_line_sp(terms[[tm]], cells))
}
## one-step: poLCA's coefficients against its first class, re-expressed
## against its intolerant class
rb1 <- rebase_polca(p1, which(lab_one_gss == "Intolerant"))
se1 <- matrix(sqrt(diag(rb1$vcov)), nrow(rb1$est), dimnames = dimnames(rb1$est))
panel_b <- c("\\midrule",
  paste("\\multicolumn{7}{l}{\\textit{One-step: classes estimated with cohort and",
    "education (\\pkg{poLCA})}} \\\\[3pt]"),
  tex_line("", sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(one_classes[-1]))),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}")
for (tm in names(terms)) {
  cols <- as.character(match(one_classes[-1], lab_one_gss))
  panel_b <- c(panel_b, tex_line_sp(terms[[tm]],
    sprintf("\\multicolumn{2}{c}{%s}", est_se(rb1$est[tm, cols], se1[tm, cols]))))
}
## one-step with three classes, against its intolerant class
rb3 <- rebase_polca(p3, which(lab_p3 == "Intolerant"))
se3 <- matrix(sqrt(diag(rb3$vcov)), nrow(rb3$est), dimnames = dimnames(rb3$est))
panel_c <- c("\\midrule",
  paste("\\multicolumn{7}{l}{\\textit{One-step with three classes, cohort and",
    "education (\\pkg{poLCA})}} \\\\[3pt]"),
  tex_line("", sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(classes3[-1]))),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}")
for (tm in names(terms)) {
  cols <- as.character(match(classes3[-1], lab_p3))
  panel_c <- c(panel_c, tex_line_sp(terms[[tm]],
    sprintf("\\multicolumn{2}{c}{%s}", est_se(rb3$est[tm, cols], se3[tm, cols]))))
}
write_table(c(panel_a, panel_b, panel_c),
  header = c(paste(" &", paste(sprintf("\\multicolumn{2}{c}{%s}", wrap_hdr(gss_classes[-1])),
    collapse = " & "), "\\\\"),
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
    paste(" &", paste(rep("Corrected & Uncorrected", 3), collapse = " & "), "\\\\")),
  align = "l cc cc cc",
  caption = paste("Effects of birth cohort and education on class membership",
    "(multinomial logit, reference class ``Intolerant''; reference categories:",
    "born after 1951, less than 12 years of education), with standard errors in",
    "parentheses. Top, corrected: ML three-step estimator with proportional",
    "assignment and standard errors corrected for the Step-1 uncertainty;",
    "uncorrected: multinomial logit of the modal class assignments. Middle and",
    "bottom: the one-step latent class regression with four and with three",
    "classes, whose classes differ from those of the measurement model",
    "(Table~\\ref{tab:gss_measurement})."),
  label = "gss_covariate")
nm_old <- sprintf("cohort1914 or before:C%d", which(lab_gss == "Intolerant of right"))


###################################################
### code chunk number 24: gss-distal
###################################################
fd_pol <- tse_distal(cl_gss, polviews ~ 1, family = "gaussian")
omnibus_test(fd_pol)
fd_race <- tse_distal(cl_gss, natrace ~ 1, family = "multinomial")


###################################################
### code chunk number 25: gss-distal-bch
###################################################
## the BCH estimates, compared with the ML estimates in the text
fd_pol_bch <- tse_distal(cl_gss, polviews ~ 1, family = "gaussian", method = "BCH")
fd_race_bch <- tse_distal(cl_gss, natrace ~ 1, family = "multinomial", method = "BCH")
mu_left <- sprintf("mu_C%d", which(lab_gss == "Intolerant of left"))
d_pol <- abs(coef(fd_pol) - coef(fd_pol_bch)[names(coef(fd_pol))])
d_pol_other <- max(d_pol[setdiff(names(d_pol), mu_left)])
d_pol_left <- d_pol[[mu_left]]
d_race <- max(abs(coef(fd_race) - coef(fd_race_bch)[names(coef(fd_race))]))
se_ratio_left <- sqrt(diag(vcov(fd_pol_bch)))[[mu_left]] / sqrt(diag(vcov(fd_pol)))[[mu_left]]
## the same models with the covariates
fc_gss_bch <- tse_covariate(cl_gss, ~ cohort + education, method = "BCH", ref = intolerant)
fb_pol_bch <- tse_distal(fc_gss_bch, polviews ~ 1, family = "gaussian")
fb_race_bch <- tse_distal(fc_gss_bch, natrace ~ 1, family = "multinomial")


###################################################
### code chunk number 26: gss-combined
###################################################
fb_pol <- tse_distal(fc_gss, polviews ~ 1, family = "gaussian")
fb_race <- tse_distal(fc_gss, natrace ~ 1, family = "multinomial")
summary(distal(fb_pol))


###################################################
### code chunk number 27: gss-distal-tables
###################################################
## two panels (distal outcome alone; with covariates), each with the ML and BCH estimators
panel_title <- function(txt, ncol) sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncol, txt)
pol_panel <- function(fits, title) {
  rows <- sapply(gss_classes, function(k) {
    nm <- sprintf("mu_C%d", which(lab_gss == k))
    tex_line(k, sapply(fits, function(f) est_se(coef(f)[nm], sqrt(diag(vcov(f)))[nm])))
  })
  c(panel_title(title, 3), unname(rows))
}
## omnibus Wald statistics of the four fits, for the captions
wald_caption <- function(fits, what) {
  w <- sapply(fits, function(f) sprintf("%.1f", omnibus_test(f)$statistic))
  df <- as.integer(omnibus_test(fits[[1]])$parameter)
  sprintf(paste("The omnibus Wald tests of equal %s across classes (%d df) are %s (ML) and",
    "%s (BCH) for the distal outcome alone and %s (ML) and %s (BCH) with covariates."),
    what, df, w[1], w[2], w[3], w[4])
}
write_table(c(pol_panel(list(fd_pol, fd_pol_bch), "Distal outcome only"), "\\midrule",
    pol_panel(list(fb_pol, fb_pol_bch), "With covariates")),
  header = "Class & ML & BCH \\\\",
  align = "l cc",
  caption = paste("Class means of liberal (1) to conservative (7) self-placement",
    "(\\code{polviews}), with standard errors in parentheses. Proportional",
    "assignment; ``with covariates'': class membership depends on cohort and",
    "education.", wald_caption(list(fd_pol, fd_pol_bch, fb_pol, fb_pol_bch), "means")),
  label = "gss_polviews")
cats <- levels(gss$natrace)
race_panel <- function(fits, title) {
  rows <- sapply(gss_classes, function(k) {
    t <- which(lab_gss == k)
    nm <- sprintf("C%d:%s", t, cats)
    tex_line(k, unlist(lapply(fits, function(f)
      est_se(coef(f)[nm], sqrt(diag(vcov(f)))[nm]))))
  })
  c(panel_title(title, 7), unname(rows))
}
cat_head <- paste(tools::toTitleCase(cats), collapse = " & ")
write_table(c(race_panel(list(fd_race, fd_race_bch), "Distal outcome only"), "\\midrule",
    race_panel(list(fb_race, fb_race_bch), "With covariates")),
  header = c(" & \\multicolumn{3}{c}{ML} & \\multicolumn{3}{c}{BCH} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
    paste("Class &", paste(rep(cat_head, 2), collapse = " & "), "\\\\")),
  align = "l ccc ccc",
  caption = paste("Class distributions of views on spending to improve the",
    "conditions of Black Americans (\\code{natrace}), with standard errors in",
    "parentheses. Proportional assignment; ``with covariates'': class membership",
    "depends on cohort and education.",
    wald_caption(list(fd_race, fd_race_bch, fb_race, fb_race_bch), "distributions")),
  label = "gss_natrace")


###################################################
### code chunk number 28: gss-onecall (eval = FALSE)
###################################################
## tseLCA(cbind(atheists, communists, militarists, racists, homosexuals) ~
##   cohort + education | polviews, data = gss, nclass = 4,
##   assignment = "proportional", ref = intolerant,
##   start = item_probs(m_gss))


###################################################
### code chunk number 29: further-fiml
###################################################
m_fiml <- tse_lca(f_elec, data = election, nclass = 3, missing = "fiml")
nobs(m_fiml)
round(class_sizes(m_fiml), 3)


###################################################
### code chunk number 30: further-newdata (eval = FALSE)
###################################################
## set.seed(20260930)
## m_76 <- tse_lca(f_gss, data = subset(gss, year == 1976), nclass = 4,
##   control = tse_control(n_init = 50))
## cl_77 <- tse_classify(m_76, newdata = subset(gss, year == 1977),
##   assignment = "proportional")


###################################################
### code chunk number 31: further-newdata-run
###################################################
set.seed(20260930)
m_76 <- fit_or_load("m_76", tse_lca(f_gss, data = subset(gss, year == 1976),
  nclass = 4, control = tse_control(n_init = N_STARTS)))
cl_77 <- tse_classify(m_76, newdata = subset(gss, year == 1977),
  assignment = "proportional")
## the same model estimated on the 1977 respondents alone, with its classes
## matched to the 1976 classes by the smallest total difference in item probabilities
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
set.seed(20260930)
m_77 <- fit_or_load("m_77", tse_lca(f_gss, data = subset(gss, year == 1977),
  nclass = 4, control = tse_control(n_init = N_STARTS)))
ip76 <- item_probs(m_76)
ip77 <- item_probs(m_77)
pm <- perms(1:4)
match77 <- pm[[which.min(sapply(pm, function(p) sum(abs(ip77[, p] - ip76))))]]
diff77 <- max(abs(ip77[, match77] - ip76))


###################################################
### code chunk number 32: save-cache
###################################################
if (cache_changed) saveRDS(cache, cache_file)
## the common table width: measure every generated table
measure <- c(list.files("tables", "_measure\\.tex$", full.names = TRUE),
  list.files(file.path("tseLCA_sim_output", "tables"), "_measure\\.tex$",
    full.names = TRUE))
writeLines(sprintf("\\input{%s}", sub("\\.tex$", "", measure)),
  file.path("tables", "table_width.tex"))


