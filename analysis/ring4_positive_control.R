# Positive-control simulation, earlier design (analysis plan fixed before running)
# Part 1, plasmode simulation: observed expression + Weibull PH simulated survival (LP0 from the observed-outcome fit + beta*S); criteria T0 template / T1 random null / T2 proxy-adjusted LRT
# Part 2, Nault 5-gene score as a real example
# R4_SMOKE=1 -> 4 simulations per beta and 10 random sets, to check that results are produced
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("R4_SMOKE", "0") == "1"
OUTD <- file.path(B, if (SMOKE) "N11_out/ring4_smoke" else "N11_out/ring4"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
NSIM <- if (SMOKE) 4L else 200L; NNULL <- if (SMOKE) 10L else 100L; NNULL_NAULT <- if (SMOKE) 20L else 1000L
BETAS <- c(0, 0.1, 0.2, 0.3, 0.4, 0.6); K <- 5L
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

# Three criteria for one signature (given an outcome)
# Direction: pergene=FALSE -> one sign for the whole signature (simulation: true genes share the sign); pergene=TRUE -> each gene takes its TCGA univariable Cox sign (Nault: includes the protective gene RAMP3)
gsign <- function(h, Y, pergene) { f <- function(s) { v <- sign(coef(coxph(Y$TCGA ~ s)))[[1]]; if (!is.finite(v) || v == 0) 1 else v }
  if (pergene) vapply(h, function(x) f(as.numeric(scale(Zc$TCGA[x, ]))), 1) else rep(f(eqscore(h, "TCGA", 1)), length(h)) }
sscore <- function(h, k, w) as.numeric(scale(colMeans(Zc[[k]][h, , drop = FALSE] * w)))
judge <- function(g, Y, nnull, seed, pergene = FALSE) {
  w0 <- gsign(g, Y, pergene); sgn <- if (pergene) NA_real_ else w0[1]
  set.seed(seed); nulls <- lapply(seq_len(nnull), function(i) sample(U, length(g)))
  nw <- lapply(nulls, function(h) gsign(h, Y, pergene))
  r <- lapply(EXT, function(k) { s <- sscore(g, k, w0); c0 <- cidx(s, Y[[k]])
    cn <- vapply(seq_len(nnull), function(i) cidx(sscore(nulls[[i]], k, nw[[i]]), Y[[k]]), 1)
    c(uni = uni_p(s, Y[[k]]), C = c0, empP = (sum(cn >= c0) + 1) / (nnull + 1), lrt = lrtp(s, Y[[k]], k)) })
  names(r) <- EXT
  data.table(T0 = all(sapply(r, `[[`, "uni") < 0.05), T1 = all(sapply(r, `[[`, "empP") < 0.05), T2 = all(sapply(r, `[[`, "lrt") < 0.05),
             C_GSE = r$GSE14520[["C"]], C_LIRI = r$LIRI[["C"]], empP_GSE = r$GSE14520[["empP"]], empP_LIRI = r$LIRI[["empP"]],
             lrt_GSE = r$GSE14520[["lrt"]], lrt_LIRI = r$LIRI[["lrt"]], sign_TCGA = sgn) }

## ---- Part 1: simulation ----
jobs <- CJ(b = seq_along(BETAS), sim = seq_len(NSIM))
one <- function(i) { j <- jobs[i]; beta <- BETAS[j$b]; seed <- 11e6 + 1e3 * j$b + j$sim; set.seed(seed)
  g <- sample(U, K)
  Y <- lapply(setNames(names(cohorts), names(cohorts)), function(k) simy(k, eqscore(g, k, 1), beta))
  r <- tryCatch(judge(g, Y, NNULL, seed + 5e5), error = function(e) data.table(err = conditionMessage(e)))
  data.table(beta = beta, sim = j$sim, genes = paste(g, collapse = ","), ev_GSE = sum(Y$GSE14520[, 2]), ev_LIRI = sum(Y$LIRI[, 2]), r) }
cl <- makeCluster(min(12, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table) })))
clusterExport(cl, c("jobs", "BETAS", "K", "U", "cohorts", "W", "Zc", "D", "CC", "BASE", "EXT", "NNULL", "simy", "eqscore", "gsign", "sscore", "cidx", "lrtp", "uni_p", "judge", "one"))
R <- rbindlist(parLapply(cl, seq_len(nrow(jobs)), one), fill = TRUE); stopCluster(cl)
fwrite(R, file.path(O, "sim_rows.tsv"), sep = "\t"); tm(sprintf("simulation rows %d, errors %d", nrow(R), if ("err" %in% names(R)) sum(!is.na(R$err)) else 0))
cpci <- function(x, n) { b <- binom.test(x, n)$conf.int; c(b[1], b[2]) }
S <- rbindlist(lapply(BETAS, function(bt) { d <- R[beta == bt & (if ("err" %in% names(R)) is.na(err) else TRUE)]; n <- nrow(d)
  rbindlist(lapply(c("T0", "T1", "T2", "FW"), function(cr) { x <- if (cr == "FW") sum(d$T1 & d$T2) else sum(d[[cr]]); ci <- cpci(x, n)
    data.table(beta = bt, HR_per_SD = exp(bt), criterion = cr, n = n, positive = x, rate = x / n, CP_lo = ci[1], CP_hi = ci[2]) })) }))
S[, criterion := factor(criterion, levels = c("T0", "T1", "T2", "FW"), labels = c("T0_template_uni", "T1_random_null", "T2_proxy_LRT", "framework_T1andT2"))]
fwrite(S, file.path(O, "sim_power.tsv"), sep = "\t")
fw <- S[criterion == "framework_T1andT2"]; mde <- fw[rate >= 0.8, min(beta)]
a2 <- fread(file.path(B, "N11_out/ring2/ring2_HR_attenuation.tsv"))[comparison == "A_vs_fitsig_k" & group == "real"]
fwrite(data.table(typeI_framework_beta0 = fw[beta == 0]$rate, typeI_CP_lo = fw[beta == 0]$CP_lo, typeI_ok_bradley = fw[beta == 0]$CP_lo <= 0.075,
                  typeI_template_beta0 = S[criterion == "T0_template_uni" & beta == 0]$rate,
                  MDE_beta_power80 = if (is.finite(mde)) mde else NA_real_, MDE_HR = if (is.finite(mde)) exp(mde) else NA_real_,
                  published_adjHR_median_GSE = exp(a2[cohort == "GSE14520"]$median_abslogHR_adj), published_adjHR_median_LIRI = exp(a2[cohort == "LIRI"]$median_abslogHR_adj)),
       file.path(O, "sim_verdict.tsv"), sep = "\t")
tm("simulation summary written")

## ---- Part 2: Nault 5-gene ----
nault <- list(main_TAF9 = c("JPT1", "RAN", "RAMP3", "KRT19", "TAF9"), sens_TAF9B = c("JPT1", "RAN", "RAMP3", "KRT19", "TAF9B"))
Yreal <- lapply(setNames(names(cohorts), names(cohorts)), function(k) cohorts[[k]]$y[CC[[k]]])
NR <- rbindlist(lapply(names(nault), function(nm) { g <- nault[[nm]]; stopifnot(all(g %in% rownames(Zc$TCGA)))
  r <- judge(g, Yreal, NNULL_NAULT, 12e6 + match(nm, names(nault)), pergene = TRUE)
  sgn <- sapply(g, function(x) sign(coef(coxph(Yreal$TCGA ~ as.numeric(scale(Zc$TCGA[x, ]))))))
  r2 <- sapply(names(cohorts), function(k) { s <- sscore(g, k, gsign(g, Yreal, TRUE)); d <- D[[k]][CC[[k]], ]
    summary(lm(s ~ prolif + stromal + immune + oxphos, data = d))$r.squared })   # R² defined as in the proxy analysis
  data.table(set = nm, genes = paste(g, collapse = ","), per_gene_TCGA_sign = paste(sprintf("%s:%+d", g, sgn), collapse = ","), r,
             R2_proxy_TCGA = r2[["TCGA"]], R2_proxy_GSE = r2[["GSE14520"]], R2_proxy_LIRI = r2[["LIRI"]]) }))
fwrite(NR, file.path(O, "nault_5gene.tsv"), sep = "\t")
print(S); print(fread(file.path(O, "sim_verdict.tsv"))); print(NR); tm("DONE")
