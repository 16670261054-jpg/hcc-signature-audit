# 图 5 示例随机签名（Ak|S004|rep1）在全部尺寸匹配随机签名（Ak）中的位置：只做核查，不改选例
# 口径与 prep_fig4c_random_showcase.R 完全相同（中位数切分 KM HR、log-rank、C 指数、TCGA 多因素 Cox）
suppressPackageStartupMessages({ library(survival); library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N28_fig/data/_load_cohorts"); source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
GA <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))[arm == "Ak"]
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE); C0 <- cohorts$TCGA; ids <- colnames(C0$z)
cov <- data.frame(Age = as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids, ph[[1]])]) / 10,
                  Sex = factor(ph$gender[match(ids, ph[[1]])]), Stage = factor(C0$stage, levels = c("I-II", "III-IV")))
one <- function(g) {
  b <- setNames(g$coef, g$gene); out <- list(k = length(b))
  for (cn in c("TCGA", "GSE14520", "LIRI")) {
    C <- cohorts[[cn]]; s <- score(b, C$z); grp <- factor(s > median(s), levels = c(FALSE, TRUE))
    out[[paste0("HR_", cn)]] <- summary(coxph(C$y ~ grp))$conf.int[1, 1]
    out[[paste0("P_", cn)]] <- 1 - pchisq(survdiff(C$y ~ grp)$chisq, 1)
    out[[paste0("C_", cn)]] <- cidx(s, C$y)
  }
  d <- cbind(Risk = as.numeric(scale(score(b, C0$z))), cov); m <- summary(coxph(C0$y ~ Risk + Age + Sex + Stage, data = d))
  out$multiHR <- m$conf.int["Risk", 1]; out$multiP <- m$coefficients["Risk", 5]; as.data.table(out)
}
R <- GA[, one(.SD), by = .(src, rep)]
fwrite(R, file.path(B, "N28_fig/data/fig5_showcase_all_Ak.tsv"), sep = "\t")
ex <- R[src == "S004" & rep == 1]; stopifnot(nrow(ex) == 1)
met <- c("HR_TCGA", "HR_GSE14520", "HR_LIRI", "C_TCGA", "C_GSE14520", "C_LIRI", "multiHR")
pct <- rbindlist(lapply(met, function(v) data.table(metric = v, example = ex[[v]], pct_rank = mean(R[[v]] <= ex[[v]]),
  median_all = median(R[[v]]), q25 = quantile(R[[v]], .25), q75 = quantile(R[[v]], .75), n = nrow(R))))
R[, all3_P05 := P_TCGA < .05 & P_GSE14520 < .05 & P_LIRI < .05]; R[, two_ext_P05 := P_GSE14520 < .05 & P_LIRI < .05]
R[, tcga_gse_P05 := P_TCGA < .05 & P_GSE14520 < .05]
pct <- rbind(pct, data.table(metric = c("frac_TCGA_P05", "frac_TCGA_GSE_P05", "frac_twoext_P05", "frac_all3_P05", "frac_multiP05"),
  example = c(ex$P_TCGA < .05, ex$P_TCGA < .05 & ex$P_GSE14520 < .05, ex$P_GSE14520 < .05 & ex$P_LIRI < .05, ex$P_TCGA < .05 & ex$P_GSE14520 < .05 & ex$P_LIRI < .05, ex$multiP < .05),
  pct_rank = NA, median_all = c(mean(R$P_TCGA < .05), mean(R$tcga_gse_P05), mean(R$two_ext_P05), mean(R$all3_P05), mean(R$multiP < .05)),
  q25 = NA, q75 = NA, n = nrow(R)), fill = TRUE)
fwrite(pct, file.path(B, "N28_fig/data/fig5_showcase_typicality.tsv"), sep = "\t")
cat("TYPICALITY OK n =", nrow(R), "\n")
