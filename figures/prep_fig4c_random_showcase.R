# Fig 4c 出图准备：一个尺寸匹配随机签名走一遍"模板论文"的标准展示（只为出图展示，不进任何判读）
# 签名 = Fig 5 选出的同一个随机签名（Ak|S004|rep 1：CEP55 / TRIM36 / FAM107A / MAST2，系数取自 signature_genes_all_arms.tsv），
#   选例规则在 fig5_singlecell.R 写死，此处不再挑选。
# 模板展示内容与口径（照抄被审计模板）：
#   KM：TCGA 训练集、GSE14520、LIRI-JP，按各自队列风险分中位数分高/低组，log-rank；
#   时间依赖 ROC（timeROC，Kaplan–Meier 权重）：TCGA 1/3/5 年；
#   多因素 Cox（TCGA）：风险分（标准化）+ 年龄 + 性别 + 分期（同 ring3_lib.R m1_indep）。
suppressPackageStartupMessages({ library(survival); library(data.table); library(timeROC) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
OD <- file.path(B, "N28_fig/data"); O <- file.path(OD, "_load_cohorts"); dir.create(O, showWarnings = FALSE)  # load_cohorts.R 会往 O 写 cohort_meta
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

## timeROC（TCGA，月）
C <- cohorts$TCGA; s <- score(b, C$z)
tr <- timeROC(T = C$time, delta = C$event, marker = s, cause = 1, weighting = "marginal", times = c(12, 36, 60), iid = FALSE)
roc <- rbindlist(lapply(seq_along(tr$times), function(i) data.table(year = tr$times[i] / 12, FPR = tr$FP[, i], TPR = tr$TP[, i])))
fwrite(roc, file.path(OD, "fig4c_timeroc_curves.tsv"), sep = "\t")
fwrite(data.table(year = tr$times / 12, AUC = tr$AUC, n_at_risk = vapply(tr$times, function(t) sum(C$time >= t), 1L)),
       file.path(OD, "fig4c_timeroc_auc.tsv"), sep = "\t")

## 多因素 Cox（TCGA；与 m1_indep 同口径）
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE); ids <- colnames(C$z)
d <- data.frame(Risk = as.numeric(scale(s)), Age = as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids, ph[[1]])]) / 10,
                Sex = factor(ph$gender[match(ids, ph[[1]])]), Stage = factor(C$stage, levels = c("I-II", "III-IV")))
m <- coxph(C$y ~ Risk + Age + Sex + Stage, data = d); sm <- summary(m)
fwrite(data.table(term = rownames(sm$conf.int), HR = sm$conf.int[, 1], lo = sm$conf.int[, 3], hi = sm$conf.int[, 4],
                  P = sm$coefficients[, 5], n = sm$n, events = sm$nevent),
       file.path(OD, "fig4c_multicox.tsv"), sep = "\t")
cat("PREP4C OK\n")
