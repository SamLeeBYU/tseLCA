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
### code chunk number 2: sim-tables
###################################################
## The simulation tables are written by tseLCA_sim.R; rebuild them from the
## saved replication results (tseLCA_sim_output/sim_replicates.rds) if needed.
if (!file.exists(file.path("tseLCA_sim_output", "tables", "sim_box.tex"))) {
  system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "tseLCA_sim.R"),
    env = c("TSELCA_SIM_TABLES_ONLY=TRUE", "TSELCA_SIM_DIR=tseLCA_sim_output"))
}


###################################################
### code chunk number 3: illus-setup
###################################################
library("poLCA")
N_STARTS <- 50L # random starts per measurement model
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
write_table <- function(body, header, align, caption, label, footer = NULL,
                        resize = FALSE) {
  rows <- c(paste0(apply(body, 1, paste, collapse = " & "), " \\\\"), footer)
  writeLines(c("\\begin{table}[t!]", "\\centering",
               if (resize) "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\begin{tabular}{%s}", align), "\\toprule", header,
               "\\midrule", rows, "\\bottomrule",
               if (resize) c("\\end{tabular}%", "}") else "\\end{tabular}",
               sprintf("\\caption{\\footnotesize %s}", caption),
               sprintf("\\label{tab:%s}", label), "\\end{table}"),
             file.path("tables", paste0(label, ".tex")))
}
prob_se <- function(model, row = NULL) {
  x <- if (is.null(row)) class_sizes(model, se = TRUE) else item_probs(model, se = TRUE)
  if (is.null(row)) rbind(est = x$estimate, se = x$se)
  else rbind(est = x$estimate[row, ], se = x$se[row, ])
}


###################################################
### code chunk number 4: elec-data
###################################################
data("election", package = "poLCA")
items <- colnames(election)[1:12]
elec <- election[complete.cases(election[, c(items, "PARTY")]), ]
f_elec <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ 1
nrow(elec)


###################################################
### code chunk number 5: elec-select (eval = FALSE)
###################################################
## set.seed(20260928)
## sel <- tse_lca(f_elec, data = elec, nclass = 1:6,
##   control = tse_control(n_init = 50))
## sel


###################################################
### code chunk number 6: elec-select-run
###################################################
set.seed(20260928)
sel <- fit_or_load("sel_elec", tse_lca(f_elec, data = elec, nclass = 1:6,
  control = tse_control(n_init = N_STARTS)))
sel


###################################################
### code chunk number 7: elec-onestep
###################################################
elec_pol <- elec
elec_pol[items] <- lapply(elec_pol[items], as.integer)
f_onestep <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ PARTY
onestep <- function(K) {
  ## each start separately: one failed start stops poLCA's own nrep loop
  best <- NULL
  for (s in seq_len(if (K == 1) 1L else N_STARTS)) {
    fit <- tryCatch(poLCA(f_onestep, elec_pol, nclass = K, maxiter = 5000,
      verbose = FALSE), error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$llik) &&
        (is.null(best) || fit$llik > best$llik)) best <- fit
  }
  best
}
set.seed(20260928)
one <- fit_or_load("onestep_elec", lapply(1:6, onestep))
ic <- as.data.frame(sel)
body <- cbind(ic$nclass, formatC(ic$BIC, format = "f", digits = 1),
  formatC(sapply(one, function(f) f$bic), format = "f", digits = 1))
write_table(body,
  header = "Classes & Three-step (measurement model) & One-step (\\pkg{poLCA}) \\\\",
  align = "c cc",
  caption = paste("BIC of the three-step measurement model and of the one-step",
    "latent class regression on \\code{PARTY} for one to six classes (ANES 2000,",
    "$N = 1300$)."),
  label = "elec_bic")


###################################################
### code chunk number 8: elec-measurement
###################################################
m_elec <- sel[[3]]
m_elec


###################################################
### code chunk number 9: elec-labels
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
s <- prob_se(m_elec)
body <- rbind(c("", "Class size", est_se(s["est", ord_elec], s["se", ord_elec])))
for (cand in c("Gore", "Bush")) for (tr in names(traits)) {
  s <- prob_se(m_elec, sprintf("P(%s%s.0|C)", tr, substr(cand, 1, 1)))
  body <- rbind(body, c(if (tr == "MORAL") cand else "", traits[[tr]],
    est_se(s["est", ord_elec], s["se", ord_elec])))
}
write_table(body,
  header = "Candidate & Trait & Gore affinity & Bush affinity & Other \\\\",
  align = "ll ccc",
  caption = paste("Three-class measurement model for the ANES 2000 candidate",
    "evaluations: class sizes and probabilities that a trait describes the",
    "candidate ``extremely well'', with standard errors in parentheses."),
  label = "elec_measurement")


###################################################
### code chunk number 10: elec-classify
###################################################
cl_elec <- tse_classify(m_elec, assignment = "modal")
cl_elec


###################################################
### code chunk number 11: elec-covariate
###################################################
other <- which(lab_elec == "Other")
fc_elec <- tse_covariate(cl_elec, ~ PARTY, ref = other)
summary(fc_elec)


###################################################
### code chunk number 12: elec-estimators
###################################################
fits_elec <- list(
  "ML, corrected SE" = fc_elec,
  "ML, robust SE" = update(fc_elec, se = "robust"),
  "BCH" = update(fc_elec, method = "BCH"),
  "two-step" = tse_twostep(m_elec, ~ PARTY, ref = other, se = TRUE),
  "uncorrected" = update(fc_elec, method = "none"))


###################################################
### code chunk number 13: elec-party-table
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
  K <- ncol(fit$coeff) + 1L
  Q <- nrow(fit$coeff)
  M <- matrix(0, K - 1L, K - 1L)
  for (i in seq_len(K - 1L)) {
    t <- seq_len(K)[-ref][i]
    if (t != 1L) M[i, t - 1L] <- 1
    if (ref != 1L) M[i, ref - 1L] <- -1
  }
  A <- kronecker(M, diag(Q))
  list(est = matrix(A %*% as.vector(fit$coeff), Q, K - 1L,
    dimnames = list(rownames(fit$coeff), seq_len(K)[-ref])),
    vcov = A %*% fit$coeff.V %*% t(A))
}
one_rb <- rebase_polca(one3, which(lab_one == "Other"))
one_se <- matrix(sqrt(diag(one_rb$vcov)), nrow(one_rb$est))
agree <- mean(lab_one[max.col(one3$posterior)] == lab_elec[classes(cl_elec)])
gb <- paste("PARTY:", c("Gore affinity", "Bush affinity"))
labels_est <- c("ML, corrected SE" = "Three-step ML, corrected SE",
  "ML, robust SE" = "Three-step ML, robust SE", "BCH" = "Three-step BCH",
  "two-step" = "Two-step", "uncorrected" = "Three-step, uncorrected")
body <- t(sapply(names(labels_est), function(e)
  c(labels_est[[e]], est_se(est[gb, e], se[gb, e], d = 3))))
one_cols <- match(c("Gore affinity", "Bush affinity"), lab_one[lab_one != "Other"])
body <- rbind(body, c("One-step (\\pkg{poLCA})",
  est_se(one_rb$est["PARTY", one_cols], one_se[2, one_cols], d = 3)))
write_table(unname(body),
  header = c(" & \\multicolumn{2}{c}{Effect of \\code{PARTY} (reference: Other)} \\\\",
    "\\cmidrule(lr){2-3}", "Estimator & Gore affinity & Bush affinity \\\\"),
  align = "l cc",
  caption = paste("Effect of party identification on membership of the Gore- and",
    "Bush-affinity classes relative to the ``Other'' class (multinomial logit",
    "coefficients, standard errors in parentheses), three-class models of the",
    "ANES 2000 data ($N = 1300$). The three-step estimators use modal assignment."),
  label = "elec_party")


###################################################
### code chunk number 14: elec-plot
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
### code chunk number 15: elec-onecall
###################################################
fit_elec <- tseLCA(cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
  MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ PARTY, data = elec,
  nclass = 3, ref = other, start = item_probs(m_elec))
all.equal(coef(fit_elec), coef(fc_elec), tolerance = 1e-4)


###################################################
### code chunk number 16: gss-download
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
### code chunk number 17: gss-data
###################################################
gss <- read.csv("gss7677_tolerance.csv")
gss$cohort <- factor(gss$cohort, 0:3,
  c("after 1951", "1934-1951", "1915-1933", "1914 or before"))
gss$education <- factor(gss$education, 0:2,
  c("< 12 years", "12 years", "> 12 years"))
gss$natrace <- factor(gss$natrace, 1:3,
  c("too little", "about right", "too much"))
f_gss <- cbind(atheists, communists, militarists, racists, homosexuals) ~ 1


###################################################
### code chunk number 18: gss-select (eval = FALSE)
###################################################
## set.seed(20260929)
## sel_gss <- tse_lca(f_gss, data = gss, nclass = 1:5,
##   control = tse_control(n_init = 50))
## m_gss <- sel_gss[[4]]


###################################################
### code chunk number 19: gss-select-run
###################################################
set.seed(20260929)
sel_gss <- fit_or_load("sel_gss", tse_lca(f_gss, data = gss, nclass = 1:5,
  control = tse_control(n_init = N_STARTS)))
m_gss <- sel_gss[[4]]


###################################################
### code chunk number 20: gss-labels
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
ord_gss <- match(gss_classes, lab_gss)
bakk_est <- rbind(size = c(.56, .23, .11, .10),
  atheists = c(.03, .98, .41, .61), communists = c(.04, .95, .59, .27),
  militarists = c(.05, .92, .34, .38), racists = c(.08, .90, .02, .81),
  homosexuals = c(.13, .96, .72, .56))
bakk_se <- rbind(size = c(.02, .01, .03, .03),
  atheists = c(.01, .01, .06, .07), communists = c(.01, .02, .11, .07),
  militarists = c(.01, .02, .05, .06), racists = c(.01, .02, .06, .20),
  homosexuals = c(.01, .01, .07, .06))
body <- NULL
for (row in rownames(bakk_est)) {
  s <- if (row == "size") prob_se(m_gss) else prob_se(m_gss, sprintf("P(%s|C)", row))
  cells <- as.vector(rbind(est_se(s["est", ord_gss], s["se", ord_gss]),
    est_se(bakk_est[row, ], bakk_se[row, ])))
  body <- rbind(body, c(if (row == "size") "Class size" else
    paste0("Tolerant of ", row), cells))
}
write_table(body,
  header = c(paste(" &", paste(sprintf("\\multicolumn{2}{c}{%s}", gss_classes),
    collapse = " & "), "\\\\"),
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
    paste(" &", paste(rep("\\pkg{tseLCA} & Bakk et al.", 4), collapse = " & "), "\\\\")),
  align = "l cc cc cc cc",
  caption = paste("Four-class model of tolerance for nonconformity (GSS 1976--77,",
    "$N = 2689$): class sizes and probabilities of a tolerant answer, with",
    "standard errors in parentheses, estimated by \\pkg{tseLCA} and reported by",
    "\\citet[Table~7]{Bakk2014}."),
  label = "gss_measurement", resize = TRUE)


###################################################
### code chunk number 21: gss-covariate
###################################################
cl_gss <- tse_classify(m_gss, assignment = "proportional")
intolerant <- which(lab_gss == "Intolerant")
fc_gss <- tse_covariate(cl_gss, ~ cohort + education, ref = intolerant)
anova(fc_gss)


###################################################
### code chunk number 22: gss-covariate-table
###################################################
unc_gss <- tse_covariate(tse_classify(m_gss, assignment = "modal"),
  ~ cohort + education, method = "none", ref = intolerant)
terms <- c("cohort1934-1951" = "Born 1934--1951", "cohort1915-1933" = "Born 1915--1933",
  "cohort1914 or before" = "Born 1914 or before",
  "education12 years" = "12 years of education",
  "education> 12 years" = "More than 12 years of education")
se_c <- sqrt(diag(vcov(fc_gss)))
se_r <- sqrt(diag(vcov(update(fc_gss, se = "robust"))))
se_u <- sqrt(diag(vcov(unc_gss)))
body <- NULL
for (tm in names(terms)) {
  cells <- character(0)
  for (k in gss_classes[-1]) {
    nm <- sprintf("%s:C%d", tm, which(lab_gss == k))
    cells <- c(cells, est_se(coef(fc_gss)[nm], se_c[nm]),
      est_se(coef(unc_gss)[nm], se_u[nm]))
  }
  body <- rbind(body, c(terms[[tm]], cells))
}
write_table(body,
  header = c(paste(" &", paste(sprintf("\\multicolumn{2}{c}{%s}", gss_classes[-1]),
    collapse = " & "), "\\\\"),
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
    paste(" &", paste(rep("Corrected & Uncorrected", 3), collapse = " & "), "\\\\")),
  align = "l cc cc cc",
  caption = paste("Effects of birth cohort and education on class membership",
    "(multinomial logit, reference class ``Intolerant''; reference categories:",
    "born after 1951, less than 12 years of education). Corrected: ML three-step",
    "estimator with proportional assignment and standard errors corrected for the",
    "Step-1 uncertainty; uncorrected: multinomial logit of the modal class",
    "assignments. Standard errors in parentheses."),
  label = "gss_covariate", resize = TRUE)
nm_old <- sprintf("cohort1914 or before:C%d", which(lab_gss == "Intolerant of right"))


###################################################
### code chunk number 23: gss-distal
###################################################
fd_pol <- tse_distal(cl_gss, polviews ~ 1, family = "gaussian")
omnibus_test(fd_pol)
fd_race <- tse_distal(cl_gss, natrace ~ 1, family = "multinomial")


###################################################
### code chunk number 24: gss-combined
###################################################
fb_pol <- tse_distal(fc_gss, polviews ~ 1, family = "gaussian")
fb_race <- tse_distal(fc_gss, natrace ~ 1, family = "multinomial")
summary(distal(fb_pol))


###################################################
### code chunk number 25: gss-distal-tables
###################################################
fd_pol_bch <- tse_distal(cl_gss, polviews ~ 1, family = "gaussian", method = "BCH")
pol_fits <- list(fd_pol, fd_pol_bch, fb_pol)
body <- NULL
for (k in gss_classes) {
  nm <- sprintf("mu_C%d", which(lab_gss == k))
  body <- rbind(body, c(k, sapply(pol_fits, function(f)
    est_se(coef(f)[nm], sqrt(diag(vcov(f)))[nm]))))
}
wald <- sapply(pol_fits, function(f) {
  h <- omnibus_test(f)
  sprintf("%.1f (%d)", h$statistic, as.integer(h$parameter))
})
write_table(body,
  footer = c("\\midrule", paste("Wald test (df) &", paste(wald, collapse = " & "), "\\\\")),
  header = c(" & \\multicolumn{2}{c}{Distal outcome only} & With covariates \\\\",
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-4}", "Class & ML & BCH & ML \\\\"),
  align = "l cc c",
  caption = paste("Class means of liberal (1) to conservative (7) self-placement",
    "(\\code{polviews}), with standard errors in parentheses, and Wald tests of",
    "equal means across classes. Proportional assignment; ``with covariates'':",
    "class membership depends on cohort and education."),
  label = "gss_polviews")
cats <- levels(gss$natrace)
race_fits <- list(fd_race, fb_race)
body <- NULL
for (k in gss_classes) {
  t <- which(lab_gss == k)
  cells <- character(0)
  for (f in race_fits) {
    nm <- sprintf("C%d:%s", t, cats)
    cells <- c(cells, est_se(coef(f)[nm], sqrt(diag(vcov(f)))[nm]))
  }
  body <- rbind(body, c(k, cells))
}
wald <- sapply(race_fits, function(f) {
  h <- omnibus_test(f)
  sprintf("\\multicolumn{3}{c}{%.1f (%d)}", h$statistic, as.integer(h$parameter))
})
cat_head <- paste(tools::toTitleCase(cats), collapse = " & ")
write_table(body,
  footer = c("\\midrule", paste("Wald test (df) &", wald[1], "&", wald[2], "\\\\")),
  header = c(paste(" & \\multicolumn{3}{c}{Distal outcome only} &",
    "\\multicolumn{3}{c}{With covariates} \\\\"),
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
    paste("Class &", cat_head, "&", cat_head, "\\\\")),
  align = "l ccc ccc",
  caption = paste("Class distributions of views on spending to improve the",
    "conditions of Black Americans (\\code{natrace}), with standard errors in",
    "parentheses, and Wald tests of equal distributions across classes. ML",
    "estimator with proportional assignment; ``with covariates'': class",
    "membership depends on cohort and education."),
  label = "gss_natrace", resize = TRUE)


###################################################
### code chunk number 26: gss-onecall (eval = FALSE)
###################################################
## tseLCA(cbind(atheists, communists, militarists, racists, homosexuals) ~
##   cohort + education | polviews, data = gss, nclass = 4,
##   assignment = "proportional", ref = intolerant, start = item_probs(m_gss))


###################################################
### code chunk number 27: further-fiml
###################################################
m_fiml <- tse_lca(f_elec, data = election, nclass = 3, missing = "fiml")
nobs(m_fiml)
round(class_sizes(m_fiml), 3)


###################################################
### code chunk number 28: further-newdata (eval = FALSE)
###################################################
## set.seed(20260930)
## m_76 <- tse_lca(f_gss, data = subset(gss, year == 1976), nclass = 4,
##   control = tse_control(n_init = 50))
## cl_77 <- tse_classify(m_76, newdata = subset(gss, year == 1977),
##   assignment = "proportional")


###################################################
### code chunk number 29: further-newdata-run
###################################################
set.seed(20260930)
m_76 <- fit_or_load("m_76", tse_lca(f_gss, data = subset(gss, year == 1976),
  nclass = 4, control = tse_control(n_init = N_STARTS)))
cl_77 <- tse_classify(m_76, newdata = subset(gss, year == 1977),
  assignment = "proportional")


###################################################
### code chunk number 30: save-cache
###################################################
if (cache_changed) saveRDS(cache, cache_file)


