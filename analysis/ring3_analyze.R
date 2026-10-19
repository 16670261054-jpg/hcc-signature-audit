# 按 N14 预注册规则汇总第 3 环（只用标准函数：binom.test / fisher.test）
suppressPackageStartupMessages({ library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out/ring3_full")
R <- fread(file.path(O, "results_m06.tsv")); R7 <- fread(file.path(O, "results_m7.tsv"))
cov <- fread(file.path(O, "published_coverage.tsv"))
full_ids <- cov[full == TRUE]$sig_id
mods <- c("M0_pos", "M1_pos", "M2_pos", "M4_pos", "M5a_pos", "M5b_pos", "M6_pos")
long <- rbind(melt(R, id.vars = c("arm", "src", "kind", "rep", "k"), measure.vars = mods, variable.name = "module", value.name = "pos"),
              R7[, .(arm, src, kind, rep, k, module = "M7_pos", pos = M7_pos)])
long[, module := as.character(module)]

summ <- function(d, label) {
  d[, {
    nn <- sum(kind == "null"); xn <- sum(pos[kind == "null"]); nr <- sum(kind == "real"); xr <- sum(pos[kind == "real"])
    bt <- binom.test(xn, nn, p = 0.05); ci <- binom.test(xn, nn)$conf.int
    ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
    .(analysis = label, n_real = nr, pos_real = xr, rate_real = xr / nr,
      n_null = nn, pos_null = xn, FPR_null = xn / nn, FPR_CI_lo = ci[1], FPR_CI_hi = ci[2],
      binom_P_vs_0.05 = bt$p.value,
      RA_typeI_inflated = ci[1] > 0.075,
      fisher_OR = unname(ft$estimate), fisher_P = ft$p.value,
      RB_real_no_advantage = ft$p.value >= 0.05)
  }, by = module]
}
out <- rbind(
  summ(long[arm == "A" & src %in% full_ids], "A_published_fullcoverage(main)"),
  summ(long[arm == "A"], "A_published_all(sensitivity)"),
  summ(long[arm == "B"], "B_mechanism_genesets_templatepipeline"))
fwrite(out, file.path(O, "ring3_decision.tsv"), sep = "\t")
print(out[, .(analysis, module, rate_real = round(rate_real, 3), FPR_null = round(FPR_null, 3),
              CI = sprintf("%.3f-%.3f", FPR_CI_lo, FPR_CI_hi), RA_typeI_inflated, fisher_P = signif(fisher_P, 3), RB_real_no_advantage)])

# 补充：M7 挑最优相对单用 LASSO 的乐观偏差（外部平均 C）
bias <- R7[, .(n = .N, best_mean_c = mean(M7_best_mean_c, na.rm = TRUE), lasso_mean_c = mean(M7_lasso_mean_c, na.rm = TRUE)), by = .(arm, kind)]
bias[, optimism := best_mean_c - lasso_mean_c]; fwrite(bias, file.path(O, "m7_selection_optimism.tsv"), sep = "\t"); print(bias)
# 补充：已发表签名在外部的 C 分布（M0）
pc <- R[arm == "A" & src %in% full_ids, .(n = .N, med_c_GSE = median(M0_c_GSE), med_c_LIRI = median(M0_c_LIRI)), by = kind]
fwrite(pc, file.path(O, "m0_c_distribution.tsv"), sep = "\t"); print(pc)
cat("DONE\n")
