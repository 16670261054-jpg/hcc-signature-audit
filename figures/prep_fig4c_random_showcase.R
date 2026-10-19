# Fig 4c preparation: standard "template article" display for one size-matched random signature (display only, no calls)
# Signature = the same random signature chosen for Fig 5 (Ak|S004|rep 1: CEP55 / TRIM36 / FAM107A / MAST2, coefficients from signature_genes_all_arms.tsv),
#   chosen by the fixed rule in fig5_singlecell.R; not reselected here.
# Display content and definitions (as in the audited template):
#   KM: TCGA training set, GSE14520, LIRI-JP, high/low groups by the median risk score of each cohort, log-rank;
#   time-dependent ROC (timeROC, Kaplan-Meier weights): TCGA 1/3/5 years;
#   multivariable Cox (TCGA): risk score (standardized) + age + sex + stage (as m1_indep in ring3_lib.R).
suppressPackageStartupMessages({ library(survival); library(data.table); library(timeROC) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
OD <- file.path(B, "figures/data"); O <- file.path(OD, "_load_cohorts"); dir.create(O, showWarnings = FALSE)  # load_cohorts.R writes cohort_meta to O
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
GA <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))
g <- GA[arm == "Ak" & src == "S004" & rep == 1]
b <- setNames(g$coef, g$gene); stopifnot(length(b) == 4, all(names(b) %in% rownames(cohorts$TCGA$z)))

km <- list(); lr <- list()
for (cn in c("TCGA", "GSE14520", "LIRI")) {
  C <- cohorts[[cn]]; s <- score(b, C$z); grp <- factor(ifelse(s > median(s), "High risk", "Low risk"), levels = c("Low risk", "High risk"))
  f <- survfit(C$y ~ grp); sf <- summary(f, censored = TRUE)
  st <- data.table(cohort = cn, strata = sub("grp=", "", as.character(sf$strata)), time = sf$time, surv = sf$surv,
                   n_risk = sf$n.risk, n_event = sf$n.event, n_censor = sf$n.censor)
  st0 <- data.table(cohort = cn, strata = levels(grp), time = 0, surv = 1, n_risk = as.integer(table(grp)), n_event = 0, n_censor = 0)
  km[[cn]] <- rbind(st0, st)
  sd <- survdiff(C$y ~ grp); hr <- summary(coxph(C$y ~ grp))$conf.int
  lr[[cn]] <- data.table(cohort = cn, n = length(s), events = sum(C$event), logrank_P = 1 - pchisq(sd$chisq, 1),
                         HR_high_vs_low = hr[1, 1], HR_lo = hr[1, 3], HR_hi = hr[1, 4], C_index = cidx(s, C$y))
}
fwrite(rbindlist(km), file.path(OD, "fig4c_km_curves.tsv"), sep = "\t")
fwrite(rbindlist(lr), file.path(OD, "fig4c_km_logrank.tsv"), sep = "\t")

## timeROC (TCGA, months)
C <- cohorts$TCGA; s <- score(b, C$z)
tr <- timeROC(T = C$time, delta = C$event, marker = s, cause = 1, weighting = "marginal", times = c(12, 36, 60), iid = FALSE)
roc <- rbindlist(lapply(seq_along(tr$times), function(i) data.table(year = tr$times[i] / 12, FPR = tr$FP[, i], TPR = tr$TP[, i])))
fwrite(roc, file.path(OD, "fig4c_timeroc_curves.tsv"), sep = "\t")
fwrite(data.table(year = tr$times / 12, AUC = tr$AUC, n_at_risk = vapply(tr$times, function(t) sum(C$time >= t), 1L)),
       file.path(OD, "fig4c_timeroc_auc.tsv"), sep = "\t")

## Multivariable Cox (TCGA; same definitions as m1_indep)
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE); ids <- colnames(C$z)
d <- data.frame(Risk = as.numeric(scale(s)), Age = as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids, ph[[1]])]) / 10,
                Sex = factor(ph$gender[match(ids, ph[[1]])]), Stage = factor(C$stage, levels = c("I-II", "III-IV")))
m <- coxph(C$y ~ Risk + Age + Sex + Stage, data = d); sm <- summary(m)
fwrite(data.table(term = rownames(sm$conf.int), HR = sm$conf.int[, 1], lo = sm$conf.int[, 3], hi = sm$conf.int[, 4],
                  P = sm$coefficients[, 5], n = sm$n, events = sm$nevent),
       file.path(OD, "fig4c_multicox.tsv"), sep = "\t")
cat("PREP4C OK\n")
