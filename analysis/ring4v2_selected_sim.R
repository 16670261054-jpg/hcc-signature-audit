# Simulation with selected candidate signatures
# Candidate signatures and gate-1 references are produced by the template pipeline on the simulated TCGA outcome; gates pooling the two external cohorts and a "missing proxy" scenario are added
# R4_SMOKE=1 -> 3 runs per scenario, 10 references
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("R4_SMOKE", "0") == "1"
OUTD <- file.path(B, if (SMOKE) "N11_out/ring4v2_smoke" else "N11_out/ring4v2"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
NSIM <- if (SMOKE) 3L else 200L; NNULL <- if (SMOKE) 10L else 100L
BETAS <- c(0, 0.2, 0.3, 0.4, 0.6); K <- 5L; NPOOL <- 100L; GAMMA_OX <- 0.3
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %5.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
O <- OUTD   # load_cohorts.R writes cohort_meta to O
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
O <- OUTD

## ---- Covariates (as in ring2_proxy.R) ----
ids <- lapply(cohorts, function(c) colnames(c$z))
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
cohorts$TCGA$age <- as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids$TCGA, ph[[1]])])
cohorts$TCGA$sex <- toupper(substr(ph$gender[match(ids$TCGA, ph[[1]])], 1, 1))
su <- fread(text = readLines(gzfile(file.path(GE, "GSE14520/GSE14520_Extra_Supplement.txt.gz"))), data.table = FALSE)
su <- su[su[["Tissue Type"]] == "Tumor", ]; su <- su[!duplicated(su$Affy_GSM), ]
cohorts$GSE14520$age <- as.numeric(su$Age[match(ids$GSE14520, su$Affy_GSM)])
cohorts$GSE14520$sex <- toupper(substr(su$Gender[match(ids$GSE14520, su$Affy_GSM)], 1, 1))
ls2 <- rd(file.path(GE, "HCCDB18_LIRI-JP/HCCDB18.sample.txt")); lp2 <- rd(file.path(GE, "HCCDB18_LIRI-JP/HCCDB18.patient.txt"))
pid <- ls2$PATIENT_ID[match(ids$LIRI, ls2$SAMPLE_ID)]
cohorts$LIRI$age <- as.numeric(lp2$AGE[match(pid, lp2$PATIENT_ID)])
cohorts$LIRI$sex <- toupper(substr(lp2$GENDER[match(pid, lp2$PATIENT_ID)], 1, 1))
suppressPackageStartupMessages(library(estimate)); utils::data("common_genes", package = "estimate", envir = environment())
plat <- c(TCGA = "illumina", GSE14520 = "affymetrix", LIRI = "illumina")
for (k in names(cohorts)) {
  X <- cohorts[[k]]$x[, ids[[k]]]; X <- X[rowSums(X > 0) > 0, ]
  tf <- tempfile(fileext = ".txt"); tg <- tempfile(fileext = ".gct"); te <- tempfile(fileext = ".gct")
  write.table(data.frame(NAME = rownames(X), X, check.names = FALSE), tf, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(tf, tg, id = "GeneSymbol"); estimateScore(tg, te, platform = plat[[k]])
  es <- read.table(te, skip = 2, header = TRUE, sep = "\t", check.names = FALSE); em <- as.matrix(es[, -(1:2)]); rownames(em) <- es$NAME
  stopifnot(ncol(em) == length(ids[[k]]))
  cohorts[[k]]$stromal <- as.numeric(scale(em["StromalScore", ])); cohorts[[k]]$immune <- as.numeric(scale(em["ImmuneScore", ]))
}
BASE <- c("stage", "prolif", "stromal", "immune", "age", "sex")
D <- lapply(cohorts, function(c) { sx <- ifelse(c$sex %in% c("M", "F"), c$sex, NA)
  data.frame(stage = factor(c$stage, levels = c("I-II", "III-IV")), prolif = as.numeric(scale(c$prolif)), oxphos = as.numeric(scale(c$oxphos)),
             stromal = c$stromal, immune = c$immune, age = c$age, sex = factor(sx, levels = c("F", "M"))) })
CC <- lapply(D, function(d) complete.cases(d[, BASE]))     # simulation and tests use only patients with complete base covariates
tm("covariates ready")

## ---- Baseline Weibull PH parameters (observed outcome) ----
fitw <- function(k) { d <- D[[k]][CC[[k]], ]; y <- cohorts[[k]]$y[CC[[k]]]
  f <- survreg(as.formula(paste("y ~", paste(BASE, collapse = "+"))), data = d, dist = "weibull")
  X <- model.matrix(as.formula(paste("~", paste(BASE, collapse = "+"))), d)
  lp0 <- drop(X[, -1, drop = FALSE] %*% (-coef(f)[-1] / f$scale))
  fc <- survreg(Surv(y[, 1], 1 - y[, 2]) ~ 1, dist = "weibull")
  list(gamma = 1 / f$scale, lambda = exp(-coef(f)[1] / f$scale), lp0 = lp0,
       cg = 1 / fc$scale, cl = exp(-coef(fc)[1] / fc$scale), tmax = max(y[, 1]), n = nrow(d), obs_event_rate = mean(y[, 2])) }
W <- lapply(setNames(names(cohorts), names(cohorts)), fitw)
fwrite(rbindlist(lapply(names(W), function(k) data.table(cohort = k, n = W[[k]]$n, gamma = W[[k]]$gamma, lambda = W[[k]]$lambda,
  cens_gamma = W[[k]]$cg, cens_lambda = W[[k]]$cl, tmax = W[[k]]$tmax, obs_event_rate = W[[k]]$obs_event_rate))), file.path(O, "baseline_weibull.tsv"), sep = "\t")
Zc <- lapply(names(cohorts), function(k) cohorts[[k]]$z[, CC[[k]]]); names(Zc) <- names(cohorts)

simy <- function(k, S, beta) { w <- W[[k]]; n <- length(S)
  T <- (-log(runif(n)) / (w$lambda * exp(w$lp0 + beta * S)))^(1 / w$gamma)
  C <- pmin((-log(runif(n)) / w$cl)^(1 / w$cg), w$tmax)
  Surv(pmin(T, C), as.integer(T <= C)) }
eqscore <- function(g, k, sgn) { s <- colMeans(Zc[[k]][g, , drop = FALSE]) * sgn; as.numeric(scale(s)) }
cidx <- function(s, y) unname(concordance(y ~ s, reverse = TRUE)$concordance)
lrtp <- function(s, y, k) { d <- cbind(D[[k]][CC[[k]], ], s = s)
  m0 <- coxph(as.formula(paste("y ~", paste(BASE, collapse = "+"))), data = d)
  m1 <- coxph(as.formula(paste("y ~", paste(c(BASE, "s"), collapse = "+"))), data = d); anova(m0, m1)[["Pr(>|Chi|)"]][2] }
uni_p <- function(s, y) summary(coxph(y ~ s))$coefficients[1, 5]
EXT <- c("GSE14520", "LIRI")


fitsig_cap <- function(g, Z, y, k) {   # as in ring3_armA_null.R
  X <- t(Z[g, , drop = FALSE])
  p <- apply(X, 2, function(v) { s <- tryCatch(summary(coxph(y ~ v))$coefficients[5], error = function(e) 1); s })
  kp <- names(p)[p < 0.05]
  if (length(kp) >= 2) {
    f <- tryCatch(suppressWarnings(cv.glmnet(X[, kp], y, family = "cox", nfolds = 5)), error = function(e) NULL)
    if (!is.null(f)) { b <- as.matrix(coef(f, s = "lambda.min"))[, 1]; b <- b[b != 0]
      if (length(b) > k) { gf <- f$glmnet.fit; lam <- gf$lambda[max(which(gf$df <= k))]; b <- as.matrix(coef(gf, s = lam))[, 1]; b <- b[b != 0] }
      if (length(b)) return(b) }
  }
  if (length(kp) == 1) return(setNames(sign(coef(coxph(y ~ X[, kp]))), kp))
  sg <- sign(coef(coxph(y ~ rowMeans(X)))); setNames(rep(sg / length(g), length(g)), g)
}
wscore <- function(b, k) as.numeric(scale(drop(t(Zc[[k]][names(b), , drop = FALSE]) %*% b)))
PD <- do.call(rbind, lapply(EXT, function(k) cbind(D[[k]][CC[[k]], ], cohort = k)))
pooled_C <- function(sl, Y) { s <- unlist(sl[EXT]); y <- do.call(rbind, lapply(EXT, function(k) as.matrix(Y[[k]]))); y <- Surv(y[, 1], y[, 2])
  unname(concordance(y ~ s + strata(PD$cohort), reverse = TRUE)$concordance) }
pooled_lrt <- function(sl, Y) { d <- cbind(PD, s = unlist(lapply(EXT, function(k) sl[[k]])))
  y <- do.call(rbind, lapply(EXT, function(k) as.matrix(Y[[k]]))); y <- Surv(y[, 1], y[, 2])
  f0 <- as.formula(paste("y ~", paste(BASE, collapse = "+"), "+ strata(cohort)")); f1 <- update(f0, . ~ . + s)
  anova(coxph(f0, data = d), coxph(f1, data = d))[["Pr(>|Chi|)"]][2] }
judge2 <- function(b, Y, refs) {
  sl <- lapply(setNames(EXT, EXT), function(k) wscore(b, k))
  rl <- lapply(refs, function(r) lapply(setNames(EXT, EXT), function(k) wscore(r, k)))
  per <- lapply(EXT, function(k) { c0 <- cidx(sl[[k]], Y[[k]]); cn <- vapply(rl, function(r) cidx(r[[k]], Y[[k]]), 1)
    c(uni = uni_p(sl[[k]], Y[[k]]), C = c0, empP = (sum(cn >= c0) + 1) / (length(cn) + 1), lrt = lrtp(sl[[k]], Y[[k]], k)) }); names(per) <- EXT
  cp <- pooled_C(sl, Y); cpn <- vapply(rl, function(r) pooled_C(r, Y), 1)
  empPp <- (sum(cpn >= cp) + 1) / (length(cpn) + 1); lp <- pooled_lrt(sl, Y)
  data.table(k_sig = length(b), T0 = all(sapply(per, `[[`, "uni") < 0.05), T1 = all(sapply(per, `[[`, "empP") < 0.05), T2 = all(sapply(per, `[[`, "lrt") < 0.05),
             T1p = empPp < 0.05, T2p = lp < 0.05, C_GSE = per$GSE14520[["C"]], C_LIRI = per$LIRI[["C"]], empP_GSE = per$GSE14520[["empP"]], empP_LIRI = per$LIRI[["empP"]],
             lrt_GSE = per$GSE14520[["lrt"]], lrt_LIRI = per$LIRI[["lrt"]], C_pooled = cp, empP_pooled = empPp, lrt_pooled = lp) }

## ---- Simulation ----
jobs <- rbind(CJ(scen = "selected", b = seq_along(BETAS), sim = seq_len(NSIM)), CJ(scen = "misspec_oxphos", b = 1L, sim = seq_len(NSIM)))
oxz <- lapply(names(cohorts), function(k) D[[k]]$oxphos[CC[[k]]]); names(oxz) <- names(cohorts)
one <- function(i) { j <- jobs[i]; beta <- BETAS[j$b]; seed <- 21e6 + (j$scen == "misspec_oxphos") * 1e6 + 1e3 * j$b + j$sim; set.seed(seed)
  g <- sample(U, K)
  Y <- lapply(setNames(names(cohorts), names(cohorts)), function(k) { S <- eqscore(g, k, 1)
    if (j$scen == "misspec_oxphos") S <- 0 * S
    w <- W[[k]]; n <- length(S); extra <- if (j$scen == "misspec_oxphos") GAMMA_OX * oxz[[k]] else 0
    T <- (-log(runif(n)) / (w$lambda * exp(w$lp0 + beta * S + extra)))^(1 / w$gamma)
    C <- pmin((-log(runif(n)) / w$cl)^(1 / w$cg), w$tmax); Surv(pmin(T, C), as.integer(T <= C)) })
  pool <- if (beta > 0) c(g, sample(setdiff(U, g), NPOOL - K)) else sample(U, NPOOL)
  r <- tryCatch({ b <- fitsig(pool, Zc$TCGA, Y$TCGA); kk <- length(b)
    set.seed(seed + 5e5); refs <- lapply(seq_len(NNULL), function(t) fitsig_cap(sample(U, NPOOL), Zc$TCGA, Y$TCGA, kk))
    cbind(judge2(b, Y, refs), n_true_in_sig = sum(names(b) %in% g)) }, error = function(e) data.table(err = conditionMessage(e)))
  data.table(scen = j$scen, beta = beta, sim = j$sim, ev_GSE = sum(Y$GSE14520[, 2]), ev_LIRI = sum(Y$LIRI[, 2]), r) }
cl <- makeCluster(as.integer(Sys.getenv("NCORE", "16")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table); library(glmnet) })))
clusterExport(cl, c("jobs", "BETAS", "K", "NPOOL", "GAMMA_OX", "U", "cohorts", "W", "Zc", "D", "CC", "BASE", "EXT", "NNULL", "PD", "oxz",
                    "eqscore", "cidx", "lrtp", "uni_p", "fitsig", "fitsig_cap", "wscore", "pooled_C", "pooled_lrt", "judge2", "one"))
R <- rbindlist(parLapply(cl, seq_len(nrow(jobs)), one), fill = TRUE); stopCluster(cl)
fwrite(R, file.path(O, "sim_rows.tsv"), sep = "\t"); nerr <- if ("err" %in% names(R)) sum(!is.na(R$err)) else 0
tm(sprintf("simulation rows %d, errors %d", nrow(R), nerr))
cpci <- function(x, n) binom.test(x, n)$conf.int
S <- rbindlist(lapply(split(R[if ("err" %in% names(R)) is.na(err) else TRUE], by = c("scen", "beta")), function(d) { n <- nrow(d)
  rbindlist(lapply(c("T0", "T1", "T2", "FW", "T1p", "T2p", "FWp"), function(cr) {
    x <- switch(cr, FW = sum(d$T1 & d$T2), FWp = sum(d$T1p & d$T2p), sum(d[[cr]])); ci <- cpci(x, n)
    data.table(scen = d$scen[1], beta = d$beta[1], HR_per_SD = exp(d$beta[1]), criterion = cr, n = n, positive = x, rate = x / n, CP_lo = ci[1], CP_hi = ci[2]) })) }))
S[, over_bradley := CP_lo > 0.075]
fwrite(S, file.path(O, "sim_power.tsv"), sep = "\t")
mde <- S[scen == "selected" & rate >= 0.8, .(MDE_beta = min(beta)), by = criterion]
fwrite(mde, file.path(O, "sim_mde.tsv"), sep = "\t")
fwrite(R[, .(median_k_sig = median(k_sig, na.rm = TRUE), mean_true_in_sig = mean(n_true_in_sig, na.rm = TRUE)), by = .(scen, beta)], file.path(O, "sim_signature_size.tsv"), sep = "\t")
print(S); print(mde); tm("DONE")
