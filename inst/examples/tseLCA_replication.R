#####################################################
### Replication script for the tseLCA package:
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee
### E-Mail: samlee@arizona.edu
#####################################################
###
### Sections
###   1. ANES 2000 (poLCA `election` data): the three-step workflow with party
###      identification as a covariate, for the three-class model of Linzer
###      and Lewis (2011, Section 5.2)
###   2. GSS 1976-77 tolerance for nonconformity (McCutcheon 1985; Bakk,
###      Oberski & Vermunt 2014): covariates, distal outcomes, and combined
###      models
###   3. One-step (poLCA) versus three-step estimation on `election`: class
###      enumeration, and the three-class models compared
###
### Only the exported tseLCA interface is used. The simulation study is
### replicated by the separate script tseLCA_sim.R.
###
### Data: `election` ships with poLCA. The GSS extract gss7677_tolerance.csv
### (public domain, CC0) comes from the replication archive of Bakk, Oberski
### & Vermunt (2014), Harvard Dataverse, doi:10.7910/DVN/24497; if it is not
### next to this script, it is rebuilt from that archive (internet needed).
###
### Run time. The measurement models are fitted from 50 random starts per
### number of classes, and the one-step poLCA models from 50 starts each;
### these fits take about 30 minutes. A full run saves them to
### tseLCA_replication_cache.rds next to this script. With the environment
### variable TSELCA_REP_QUICK=TRUE, the script instead loads them from that
### file and recomputes everything else (all Step-2 and Step-3 models,
### tables, and tests) in about a minute:
###   TSELCA_REP_QUICK=TRUE Rscript tseLCA_replication.R
#####################################################

###################################################
### preliminaries
###################################################
library("tseLCA")
library("poLCA")
options(digits = 4, width = 80)
N_STARTS <- 50L # random starts per measurement model

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(normalizePath(f)) else getwd()
  }
)

## The time-consuming fits: computed and cached in a full run, loaded from
## the cache in a quick run.
QUICK <- isTRUE(as.logical(Sys.getenv("TSELCA_REP_QUICK", "FALSE")))
cache_file <- file.path(script_dir, "tseLCA_replication_cache.rds")
if (QUICK) {
  if (!file.exists(cache_file)) {
    stop("TSELCA_REP_QUICK=TRUE needs ", cache_file, "; run the script once without it.")
  }
  cache <- readRDS(cache_file)
  message("Quick run: the measurement-model and poLCA fits are loaded from ", cache_file)
} else {
  cache <- list()
}
fit_or_load <- function(name, expr) {
  if (QUICK) return(cache[[name]])
  cache[[name]] <<- expr
  cache[[name]]
}

###################################################
### 1. ANES 2000: data
###################################################
## Twelve evaluations of the candidates Gore (G) and Bush (B): how well is
## each described as moral, caring, knowledgeable, a good leader, dishonest,
## and intelligent (1 extremely well ... 4 not well at all), and party
## identification (PARTY: 1 strong Democrat ... 7 strong Republican).
## Complete cases on the items and PARTY, so that Sections 1 and 3 use the
## same sample.
data("election", package = "poLCA")
items <- colnames(election)[1:12]
elec <- election[complete.cases(election[, c(items, "PARTY")]), c(items, "PARTY")]
nrow(elec)
f_elec <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
                MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ 1

## Class labels following Linzer and Lewis (2011): the class that rates Gore
## best relative to Bush ("Gore affinity"), the reverse ("Bush affinity"),
## and the rest ("Other"); from posterior-weighted mean ratings of the five
## positive traits (lower = better).
G <- c("MORALG", "CARESG", "KNOWG", "LEADG", "INTELG")
B <- c("MORALB", "CARESB", "KNOWB", "LEADB", "INTELB")
label_elec <- function(post) {
  codes <- sapply(elec[c(G, B)], as.integer)
  means <- crossprod(post, codes) / colSums(post)
  score <- rowMeans(means[, G]) - rowMeans(means[, B])
  lab <- rep("Other", ncol(post))
  lab[which.min(score)] <- "Gore affinity"
  lab[which.max(score)] <- "Bush affinity"
  lab
}

###################################################
### 1a. Step 1: the measurement model
###################################################
## The measurement model alone defines the classes: covariates and distal
## outcomes play no part in choosing the number of classes.
set.seed(20260928)
sel <- fit_or_load("sel_elec", tse_lca(f_elec, data = elec, nclass = 1:6,
                                       control = tse_control(n_init = N_STARTS)))
sel
plot(sel)

## The information criteria keep decreasing up to six classes. As Linzer and
## Lewis (2011), we analyze the three-class model of Gore supporters, Bush
## supporters, and a neutral group.
m_elec <- sel[[3]]
lab_elec <- label_elec(posterior(m_elec))
lab_elec
m_elec
round(item_probs(m_elec), 3)
plot(m_elec)
logLik(m_elec)
BIC(m_elec)

###################################################
### 1b. Step 2: class assignment and classification error
###################################################
cl_elec <- tse_classify(m_elec, assignment = "modal")
cl_elec # classification-error matrix P(W = s | X = t) and entropy R^2
head(posterior(cl_elec))
table(classes(cl_elec))

###################################################
### 1c. Step 3: party identification as a covariate
###################################################
other <- which(lab_elec == "Other") # reference class
fc_elec <- tse_covariate(cl_elec, ~ PARTY, method = "ML", se = "corrected", ref = other)
summary(fc_elec)
confint(fc_elec)
anova(fc_elec)
logLik(fc_elec)
AIC(fc_elec)
printCoefmat(coef(summary(fc_elec)))

## class membership probabilities along the party scale (cf. Linzer and
## Lewis 2011, Figure 2)
p_party <- predict(fc_elec, newdata = data.frame(PARTY = 1:7))
colnames(p_party) <- lab_elec
round(p_party, 3)
matplot(1:7, p_party, type = "l", lwd = 3, lty = 1:3, col = 1, ylim = c(0, 1),
        xlab = "Party ID: strong Democrat (1) to strong Republican (7)",
        ylab = "Probability of latent class membership")
legend("right", lab_elec, lty = 1:3, lwd = 3, bty = "n")

## another reference class
summary(relevel(fc_elec, ref = which(lab_elec == "Gore affinity")))

## the same model from other estimators
fits_elec <- list(
  "ML, corrected SE" = fc_elec,
  "ML, robust SE" = tse_covariate(cl_elec, ~ PARTY, se = "robust", ref = other),
  "BCH" = tse_covariate(cl_elec, ~ PARTY, method = "BCH", ref = other),
  "two-step" = tse_twostep(m_elec, ~ PARTY, se = TRUE, ref = other),
  "uncorrected" = tse_covariate(cl_elec, ~ PARTY, method = "none", ref = other)
)
party <- grep("^PARTY:", names(coef(fc_elec)), value = TRUE)
est <- sapply(fits_elec, function(f) coef(f)[party])
se <- sapply(fits_elec, function(f) sqrt(diag(vcov(f)))[party])
rownames(est) <- rownames(se) <- paste("PARTY:", lab_elec[-other])
round(est, 3)
round(se, 3)

## all three steps in one call
## (starting from the chosen measurement model, so the classes are the same)
fit_elec <- tseLCA(update(f_elec, . ~ PARTY), data = elec, nclass = 3,
                   ref = other, start = item_probs(m_elec))
fit_elec
all.equal(coef(fit_elec), coef(fc_elec), tolerance = 1e-4)

## the methods available for fitted models
methods(class = "tseLCA")
methods(class = "tseLCA_structural")

###################################################
### 2. GSS 1976-77: data
###################################################
## Five tolerance indicators (1 = would allow a communist, atheist,
## homosexual, militarist, or racist to speak, to teach, and would not remove
## their book from the library), birth cohort, education, and two outcomes.
gss_file <- file.path(script_dir, "gss7677_tolerance.csv")
if (!file.exists(gss_file)) {
  sav <- tempfile(fileext = ".sav")
  utils::download.file(
    "https://dataverse.harvard.edu/api/access/datafile/2456935?format=original",
    sav, mode = "wb", quiet = TRUE
  )
  raw <- suppressWarnings(foreign::read.spss(sav, to.data.frame = TRUE,
                                             use.value.labels = FALSE))
  write.csv(data.frame(
    year = raw$year, id = raw$id, atheists = raw$atheists_,
    communists = raw$communists_, homosexuals = raw$homtest,
    militarists = raw$militarists_, racists = raw$racist_,
    cohort = raw$cohort_, education = raw$education, polviews = raw$polviews,
    partyid = raw$partyid, natrace = raw$natrace
  ), gss_file, row.names = FALSE, na = "")
}
gss <- read.csv(gss_file)
gss$cohort <- factor(gss$cohort, 0:3,
                     c("after 1951", "1934-1951", "1915-1933", "1914 or before"))
gss$education <- factor(gss$education, 0:2, c("< 12 years", "12 years", "> 12 years"))
gss$natrace <- factor(gss$natrace, 1:3, c("too little", "about right", "too much"))
f_gss <- cbind(atheists, communists, militarists, racists, homosexuals) ~ 1

###################################################
### 2a. Step 1: enumeration and the four-class model
###################################################
## With five binary items, at most five classes are identified.
set.seed(20260929)
sel_gss <- fit_or_load("sel_gss", tse_lca(f_gss, data = gss, nclass = 1:5,
                                          control = tse_control(n_init = N_STARTS)))
sel_gss
m_gss <- sel_gss[[4]] # as McCutcheon (1985) and Bakk et al. (2014)
ip <- item_probs(m_gss)
rownames(ip) <- c("atheists", "communists", "militarists", "racists", "homosexuals")

## class labels from the profiles
avg <- colMeans(ip)
mid <- setdiff(1:4, c(which.max(avg), which.min(avg)))
lab_gss <- character(4)
lab_gss[which.max(avg)] <- "Tolerant"
lab_gss[which.min(avg)] <- "Intolerant"
lab_gss[mid] <- ifelse(colMeans(ip[c("atheists", "communists"), mid]) >
                         colMeans(ip[c("militarists", "racists"), mid]),
                       "Intolerant of right", "Intolerant of left")
tab_gss <- rbind(size = class_sizes(m_gss), ip)
colnames(tab_gss) <- lab_gss
round(tab_gss, 3)

###################################################
### 2b. Steps 2 and 3: cohort and education
###################################################
cl_gss <- tse_classify(m_gss, assignment = "proportional")
cl_gss
intolerant <- which(lab_gss == "Intolerant")
fc_gss <- tse_covariate(cl_gss, ~ cohort + education, ref = intolerant)
summary(fc_gss)
anova(fc_gss)

## corrected vs uncorrected estimates (uncorrected: modal assignment, as in
## McCutcheon 1985)
unc_gss <- tse_covariate(tse_classify(m_gss, assignment = "modal"),
                         ~ cohort + education, method = "none", ref = intolerant)
slopes <- grep("^(cohort|education)", names(coef(fc_gss)), value = TRUE)
round(cbind(corrected = coef(fc_gss)[slopes], uncorrected = coef(unc_gss)[slopes],
            ratio = coef(fc_gss)[slopes] / coef(unc_gss)[slopes]), 3)

## cohort x education interaction: likelihood-ratio test. (The Wald test is
## unreliable here: one class probability is on the boundary.)
fc_int <- tse_covariate(cl_gss, ~ cohort * education, ref = intolerant)
lr <- 2 * (as.numeric(logLik(fc_int)) - as.numeric(logLik(fc_gss)))
lr_df <- attr(logLik(fc_int), "df") - attr(logLik(fc_gss), "df")
c(LR = lr, df = lr_df, p = pchisq(lr, lr_df, lower.tail = FALSE))

###################################################
### 2c. Distal outcomes
###################################################
## liberal (1) - conservative (7) self-placement
fd_pol <- tse_distal(cl_gss, polviews ~ 1, family = "gaussian")
summary(fd_pol)
omnibus_test(fd_pol)
summary(tse_distal(cl_gss, polviews ~ 1, family = "gaussian", method = "BCH"))

## spending on improving the conditions of Blacks (too little, about right,
## too much)
fd_race <- tse_distal(cl_gss, natrace ~ 1, family = "multinomial")
summary(fd_race)
omnibus_test(fd_race)

###################################################
### 2d. Distal outcomes combined with the covariates
###################################################
fb_pol <- tse_distal(fc_gss, polviews ~ 1, family = "gaussian")
summary(fb_pol)
omnibus_test(fb_pol)
fb_race <- tse_distal(fc_gss, natrace ~ 1, family = "multinomial")
omnibus_test(fb_race)

## in one call
fit_gss <- tseLCA(cbind(atheists, communists, militarists, racists, homosexuals) ~
                    cohort + education | polviews,
                  data = gss, nclass = 4, family = "gaussian",
                  assignment = "proportional", ref = intolerant,
                  start = item_probs(m_gss))
summary(distal(fit_gss))
all.equal(coef(distal(fit_gss)), coef(distal(fb_pol)), tolerance = 1e-4)

###################################################
### 3. One-step vs three-step estimation (election)
###################################################
## poLCA estimates the classes jointly with the covariate (one-step); its
## BIC refers to that joint model. Each start is run separately because a
## single failed start stops poLCA's own nrep loop.
elec_pol <- elec
elec_pol[items] <- lapply(elec_pol[items], as.integer)
f_onestep <- update(f_elec, . ~ PARTY)
onestep <- function(K) {
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

## 3a. class enumeration: BIC of the measurement model (three-step) and of
## the joint model (one-step)
bic <- data.frame(K = 1:6, three_step = as.data.frame(sel)$BIC,
                  one_step = sapply(one, function(f) f$bic))
bic
matplot(bic$K, bic[, -1], type = "b", pch = c(16, 1), lty = 1:2, col = 1,
        xlab = "Number of classes", ylab = "BIC")
legend("topright", c("three-step (measurement model)", "one-step (poLCA)"),
       pch = c(16, 1), lty = 1:2, bty = "n")

## 3b. the three-class models: modal assignments and PARTY effects
one3 <- one[[3]]
lab_one <- label_elec(one3$posterior)
table(one_step = lab_one[max.col(one3$posterior)],
      three_step = lab_elec[classes(cl_elec)])
mean(lab_one[max.col(one3$posterior)] == lab_elec[classes(cl_elec)])

## poLCA's coefficients are against its first class; re-express them against
## the "Other" class, with their variance, to compare with tseLCA.
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
nonref <- paste("PARTY:", lab_one[lab_one != "Other"])
cmp <- rbind(
  one_step = c(one_rb$est["PARTY", ], one_se[2, ]),
  three_step_ML = c(est[nonref, "ML, corrected SE"], se[nonref, "ML, corrected SE"])
)
colnames(cmp) <- c(paste("estimate", nonref), paste("se", nonref))
round(cmp, 3)

## predicted class membership along the party scale: three-step (tseLCA)
## and one-step (poLCA)
party_x <- seq(1, 7, length.out = 101)
p_three <- predict(fc_elec, newdata = data.frame(PARTY = party_x))
colnames(p_three) <- lab_elec
eta_one <- cbind(0, cbind(1, party_x) %*% one3$coeff)
p_one <- exp(eta_one) / rowSums(exp(eta_one))
colnames(p_one) <- lab_one
p_one <- p_one[, lab_elec]
matplot(party_x, p_three, type = "l", lwd = 3, col = 1, lty = 1, ylim = c(0, 1),
        xlab = "Party ID: strong Democrat (1) to strong Republican (7)",
        ylab = "Probability of latent class membership",
        main = "Party ID as a predictor of candidate affinity class")
matlines(party_x, p_one, lwd = 2, col = 1, lty = 2)
peak <- apply(p_three, 2, which.max)
for (k in seq_along(peak)) { # labels at the peaks, kept inside the plot
  text(party_x[peak[k]], max(p_three[, k]) + 0.05, lab_elec[k],
       adj = c((peak[k] - 1) / (length(party_x) - 1), 0.5))
}
legend("topright", c("tseLCA (three-step)", "poLCA (one-step)"),
       lty = 1:2, lwd = c(3, 2), col = 1, bty = "n")

if (!QUICK) {
  saveRDS(cache, cache_file)
  message("Saved the measurement-model and poLCA fits to ", cache_file)
}

sessionInfo()
