# The 148 published signatures through the two gates
# Candidate signatures and gate-1 references are produced by the template pipeline on the simulated TCGA outcome; gates pooling the two external cohorts and a "missing proxy" scenario are added
# R4_SMOKE=1 -> 3 runs per scenario, 10 references
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("R4_SMOKE", "0") == "1"
OUTD <- file.path(B, if (SMOKE) "N11_out/n2_gates_smoke" else "N11_out/n2_gates"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
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


## ---- Published signatures (as in ring2_proxy.R / hiend_m8_formal.R: reported coefficient x TCGA sd, >= 2 genes in Uall) ----
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
L <- L[sig_id %in% I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id & !is.na(coef)]
sdT <- apply(cohorts$TCGA$x[Uall, colnames(cohorts$TCGA$z)], 1, sd)
pub <- lapply(split(L, L$sig_id), function(d) { d <- d[gene_hgnc %in% Uall]; d <- d[!duplicated(gene_hgnc)]; setNames(d$coef * sdT[d$gene_hgnc], d$gene_hgnc) })
pub <- pub[vapply(pub, length, 1L) >= 2]; stopifnot(length(pub) == 148)
if (SMOKE) pub <- pub[1:3]
Yreal <- lapply(setNames(names(cohorts), names(cohorts)), function(k) cohorts[[k]]$y[CC[[k]]])
one <- function(i) { nm <- names(pub)[i]; b <- pub[[nm]]; seed <- 31e6 + 1e3 * i
  r <- tryCatch({ set.seed(seed); refs <- lapply(seq_len(NNULL), function(t) fitsig_cap(sample(U, NPOOL), Zc$TCGA, Yreal$TCGA, length(b)))
    judge2(b, Yreal, refs) }, error = function(e) data.table(err = conditionMessage(e)))
  data.table(src = nm, k_pub = length(b), r) }
cl <- makeCluster(as.integer(Sys.getenv("NCORE", "6")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table); library(glmnet) })))
clusterExport(cl, c("pub", "Yreal", "U", "Zc", "D", "CC", "BASE", "EXT", "NNULL", "NPOOL", "PD", "cidx", "lrtp", "uni_p", "fitsig_cap", "wscore", "pooled_C", "pooled_lrt", "judge2", "one"))
R <- rbindlist(parLapply(cl, seq_along(pub), one), fill = TRUE); stopCluster(cl)
fwrite(R, file.path(O, "published_gates.tsv"), sep = "\t")
nerr <- if ("err" %in% names(R)) sum(!is.na(R$err)) else 0; tm(sprintf("rows %d errors %d", nrow(R), nerr))
d <- R[if ("err" %in% names(R)) is.na(err) else TRUE]; n <- nrow(d)
cnt <- function(lab, x) { ci <- binom.test(sum(x), n)$conf.int; data.table(criterion = lab, n = n, positive = sum(x), rate = mean(x), CP_lo = ci[1], CP_hi = ci[2]) }
S <- rbind(cnt("template_T0", d$T0), cnt("gate1_separate", d$T1), cnt("gate2_separate", d$T2), cnt("both_separate", d$T1 & d$T2),
           cnt("gate1_pooled", d$T1p), cnt("gate2_pooled", d$T2p), cnt("both_pooled", d$T1p & d$T2p))
fwrite(S, file.path(O, "published_gates_summary.tsv"), sep = "\t"); print(S); tm("DONE")
