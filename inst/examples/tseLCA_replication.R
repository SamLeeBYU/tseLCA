#####################################################
### Replication script for the tseLCA package:
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee
### E-Mail: samlee@arizona.edu
#####################################################
###
### Sections
###   1. ANES 2000 (poLCA `election` data): choosing the number of classes
###      from the measurement model, and the three-step workflow with a
###      covariate (PARTY)
###   2. GSS 1976-77 tolerance for nonconformity (McCutcheon 1985; Bakk,
###      Oberski & Vermunt 2014): covariates, distal outcomes, and combined
###      models
###   3. One-step (poLCA) versus three-step class enumeration on `election`
###
### Only the exported tseLCA interface is used. The simulation study is
### replicated by the separate script tseLCA_sim.R.
###
### Data: `election` ships with poLCA. The GSS extract gss7677_tolerance.csv
### (public domain, CC0) comes from the replication archive of Bakk, Oberski
### & Vermunt (2014), Harvard Dataverse, doi:10.7910/DVN/24497; if it is not
### next to this script, it is rebuilt from that archive (internet needed).
###
### Run time: about 15 minutes, most of it in Section 3 (poLCA with random
### starts). Set N_STARTS lower for a quicker, less thorough run.
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

###################################################
### 1. ANES 2000: data
###################################################
## Twelve evaluations of the candidates (Gore and Bush) on four-point scales
## and party identification (PARTY: 1 strong Democrat ... 7 strong Republican).
## Complete cases on the items and PARTY, so that Sections 1 and 3 use the
## same sample.
data("election", package = "poLCA")
items <- colnames(election)[1:12]
elec <- election[complete.cases(election[, c(items, "PARTY")]), c(items, "PARTY")]
nrow(elec)
f_elec <- cbind(MORALG, CARESG, KNOWG, LEADG, DISHONG, INTELG,
                MORALB, CARESB, KNOWB, LEADB, DISHONB, INTELB) ~ 1

###################################################
### 1a. Step 1: number of classes from the measurement model
###################################################
## The measurement model alone determines the classes; covariates and distal
## outcomes play no part in choosing K.
set.seed(20260928)
sel <- tse_lca(f_elec, data = elec, nclass = 1:6,
               control = tse_control(n_init = N_STARTS))
sel
plot(sel)
ic <- as.data.frame(sel)
K_elec <- ic$nclass[which.min(ic$BIC)]
m_elec <- best_model(sel, criterion = "BIC")

###################################################
### 1b. Inspecting the measurement model
###################################################
m_elec
class_sizes(m_elec)
round(item_probs(m_elec), 3)
plot(m_elec)
logLik(m_elec)
BIC(m_elec)

###################################################
### 1c. Step 2: class assignment and classification error
###################################################
cl_elec <- tse_classify(m_elec, assignment = "modal")
cl_elec # classification-error matrix P(W = s | X = t) and entropy R^2
head(posterior(cl_elec))
table(classes(cl_elec))

###################################################
### 1d. Step 3: party identification as a covariate
###################################################
fc_elec <- tse_covariate(cl_elec, ~ PARTY, method = "ML", se = "corrected")
summary(fc_elec)
confint(fc_elec)
anova(fc_elec)
logLik(fc_elec)
AIC(fc_elec)

## class membership probabilities along the party scale
round(predict(fc_elec, newdata = data.frame(PARTY = 1:7)), 3)

## another reference class
summary(relevel(fc_elec, ref = 2))

## the coefficients as a matrix of the familiar printCoefmat() form
printCoefmat(coef(summary(fc_elec)))

## the same model from other estimators
fits_elec <- list(
  "ML, corrected SE" = fc_elec,
  "ML, robust SE" = tse_covariate(cl_elec, ~ PARTY, method = "ML", se = "robust"),
  "BCH" = tse_covariate(cl_elec, ~ PARTY, method = "BCH"),
  "two-step" = tse_twostep(m_elec, ~ PARTY, se = TRUE),
  "uncorrected" = tse_covariate(cl_elec, ~ PARTY, method = "none")
)
party <- grep("^PARTY:", names(coef(fc_elec)), value = TRUE)
round(sapply(fits_elec, function(f) coef(f)[party]), 3)
round(sapply(fits_elec, function(f) sqrt(diag(vcov(f)))[party]), 3)

## all three steps in one call
## (starting from the chosen measurement model, so the classes are the same)
fit_elec <- tseLCA(update(f_elec, . ~ PARTY), data = elec, nclass = K_elec,
                   start = item_probs(m_elec))
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
## their book from the library), birth cohort, education, and three outcomes.
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
sel_gss <- tse_lca(f_gss, data = gss, nclass = 1:5,
                   control = tse_control(n_init = N_STARTS))
sel_gss
m_gss <- sel_gss[[4]] # as McCutcheon (1985) and Bakk et al. (2014)
ip <- item_probs(m_gss)
rownames(ip) <- c("atheists", "communists", "militarists", "racists", "homosexuals")

## class labels from the profiles
avg <- colMeans(ip)
mid <- setdiff(1:4, c(which.max(avg), which.min(avg)))
lab <- character(4)
lab[which.max(avg)] <- "Tolerant"
lab[which.min(avg)] <- "Intolerant"
lab[mid] <- ifelse(colMeans(ip[c("atheists", "communists"), mid]) >
                     colMeans(ip[c("militarists", "racists"), mid]),
                   "Intolerant of right", "Intolerant of left")
round(rbind(size = class_sizes(m_gss), ip), 3)
lab

###################################################
### 2b. Steps 2 and 3: cohort and education
###################################################
cl_gss <- tse_classify(m_gss, assignment = "proportional")
cl_gss
intolerant <- which(lab == "Intolerant")
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
### 3. One-step vs three-step class enumeration (election)
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
one <- lapply(1:6, onestep)
bic <- data.frame(K = 1:6, three_step = ic$BIC,
                  one_step = sapply(one, function(f) f$bic))
bic
matplot(bic$K, bic[, -1], type = "b", pch = c(16, 1), lty = 1:2, col = 1,
        xlab = "Number of classes", ylab = "BIC")
legend("topright", c("three-step (measurement model)", "one-step (poLCA)"),
       pch = c(16, 1), lty = 1:2, bty = "n")

sessionInfo()
