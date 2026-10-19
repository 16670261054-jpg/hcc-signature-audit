# Exploratory pilot: genome-wide gating (continue/stop rule fixed before running)
# Per-gene multivariable Cox (Wald); gate 1 = TCGA adjusted P genome-wide BH-FDR < 0.05; gate 2 = same direction and P < 0.05 in GSE14520 and LIRI-JP
# Covariates: stage + proliferation (mean z of E2F and G2M) + ESTIMATE Stromal/Immune (TumorPurity only for Affymetrix -> two scores in all cohorts) + age + sex
# Models: main; S1 without immune; S2 with OXPHOS; S3 stage only
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/n16_gate"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %5.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## ---- Age / sex ----
ids <- lapply(cohorts, function(c) colnames(c$z))
for (k in names(cohorts)) stopifnot(length(cohorts[[k]]$time) == ncol(cohorts[[k]]$z), length(cohorts[[k]]$stage) == ncol(cohorts[[k]]$z),
  identical(colnames(cohorts[[k]]$z), colnames(cohorts[[k]]$x)))
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

## ---- ESTIMATE (computed within each cohort) ----
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
cov <- rbindlist(lapply(names(cohorts), function(k) { c <- cohorts[[k]]
  data.table(cohort = k, n = length(ids[[k]]), events = sum(c$event), stage_NA = sum(is.na(c$stage)), age_NA = sum(is.na(c$age)),
             sex_NA = sum(!c$sex %in% c("M", "F")), r_prolif_immune = cor(c$prolif, c$immune), r_prolif_stromal = cor(c$prolif, c$stromal)) }))
fwrite(cov, file.path(O, "covariate_audit.tsv"), sep = "\t"); print(cov); tm("covariates ready")

## ---- Per-gene Cox ----
mods <- list(main = c("stage", "prolif", "stromal", "immune", "age", "sex"), S1_noimmune = c("stage", "prolif", "stromal", "age", "sex"),
             S2_plusoxphos = c("stage", "prolif", "stromal", "immune", "age", "sex", "oxphos"), S3_stageonly = c("stage"))
D <- lapply(cohorts, function(c) { sx <- ifelse(c$sex %in% c("M", "F"), c$sex, NA)
  data.frame(time = c$time, event = c$event, stage = factor(c$stage, levels = c("I-II", "III-IV")), prolif = c$prolif, oxphos = c$oxphos,
             stromal = c$stromal, immune = c$immune, age = c$age, sex = factor(sx, levels = c("F", "M"))) })
Z <- lapply(cohorts, function(c) c$z[U, ])
fit1 <- function(g) { out <- list()
  for (k in names(D)) for (m in names(mods)) { d <- D[[k]]; d$g <- Z[[k]][g, ]
    f <- tryCatch(coxph(as.formula(paste("Surv(time, event) ~ g +", paste(mods[[m]], collapse = " + "))), data = d), error = function(e) NULL)
    cf <- if (is.null(f)) rep(NA, 3) else summary(f)$coefficients["g", c("coef", "se(coef)", "Pr(>|z|)")]
    out[[length(out) + 1]] <- list(gene = g, cohort = k, model = m, logHR = cf[1], se = cf[2], P = cf[3], n = if (is.null(f)) NA else f$n) }
  rbindlist(out) }
cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table) })))
clusterExport(cl, c("D", "Z", "mods", "fit1"))
R <- rbindlist(parLapply(cl, U, fit1)); stopCluster(cl)
fwrite(R, file.path(O, "gene_cox_all.tsv.tmp"), sep = "\t"); file.rename(file.path(O, "gene_cox_all.tsv.tmp"), file.path(O, "gene_cox_all.tsv"))
tm(sprintf("cox done: %d rows, NA %d", nrow(R), sum(is.na(R$P))))

## ---- Gates ----
W <- dcast(R, gene + model ~ cohort, value.var = c("logHR", "se", "P"))
W[, FDR_TCGA := p.adjust(P_TCGA, "BH"), by = model]
W[, gate1 := FDR_TCGA < 0.05]
W[, gate2 := sign(logHR_GSE14520) == sign(logHR_TCGA) & P_GSE14520 < 0.05 & sign(logHR_LIRI) == sign(logHR_TCGA) & P_LIRI < 0.05]
W[, pass := gate1 & gate2 %in% TRUE]
mt <- function(b, s) { w <- 1 / s^2; m <- sum(w * b) / sum(w); q <- sum(w * (b - m)^2); list(meta_logHR = m, meta_se = sqrt(1 / sum(w)), Q = q, I2 = max(0, (q - 2) / q)) }
W[, c("meta_logHR", "meta_se", "Q", "I2") := mt(c(logHR_TCGA, logHR_GSE14520, logHR_LIRI), c(se_TCGA, se_GSE14520, se_LIRI)), by = .(gene, model)]
W[, meta_HR := exp(meta_logHR)]
fwrite(W, file.path(O, "gene_gate_wide.tsv"), sep = "\t")
cnt <- W[, .(n_genes = .N, gate1 = sum(gate1, na.rm = TRUE), gate2_given_gate1 = sum(gate1 & gate2, na.rm = TRUE),
             pass = sum(pass, na.rm = TRUE), pass_risk = sum(pass & logHR_TCGA > 0, na.rm = TRUE), pass_protective = sum(pass & logHR_TCGA < 0, na.rm = TRUE)), by = model]
fwrite(cnt, file.path(O, "gate_counts.tsv"), sep = "\t"); print(cnt)

## ---- Genes passing the gates: PH assumption test (main model) ----
P1 <- W[model == "main" & pass == TRUE][order(P_TCGA)]
if (nrow(P1)) {
  zph <- rbindlist(lapply(P1$gene, function(g) rbindlist(lapply(names(D), function(k) { d <- D[[k]]; d$g <- Z[[k]][g, ]
    f <- coxph(as.formula(paste("Surv(time, event) ~ g +", paste(mods$main, collapse = " + "))), data = d)
    z <- cox.zph(f)$table; data.table(gene = g, cohort = k, zph_P_gene = z["g", "p"], zph_P_global = z["GLOBAL", "p"]) }))))
  fwrite(zph, file.path(O, "pass_main_coxzph.tsv"), sep = "\t")
  P1 <- merge(P1, dcast(zph, gene ~ cohort, value.var = "zph_P_gene"), by = "gene", suffixes = c("", "_zph"))
}
fwrite(P1[order(P_TCGA)], file.path(O, "gate_pass_main.tsv"), sep = "\t")
tm("DONE")
