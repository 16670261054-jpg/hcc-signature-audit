# A 组按预注册对照重新判读：真实已发表签名（来自 ring3_full）vs 两种模板流程随机对照（fitsig_k 主分析 / fitsig_pool100 敏感性）
suppressPackageStartupMessages({ library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out")
cov <- fread(file.path(O, "ring3_full/published_coverage.tsv")); full_ids <- cov[full == TRUE]$sig_id
mods <- c("M0_pos", "M1_pos", "M2_pos", "M4_pos", "M5a_pos", "M5b_pos", "M6_pos")
lng <- function(R, R7) { x <- rbind(melt(R, id.vars = c("arm", "src", "kind", "rep", "k"), measure.vars = mods, variable.name = "module", value.name = "pos"),
                                    R7[, .(arm, src, kind, rep, k, module = "M7_pos", pos = M7_pos)]); x[, module := as.character(module)]; x }
real <- lng(fread(file.path(O, "ring3_full/results_m06.tsv")), fread(file.path(O, "ring3_full/results_m7.tsv")))[arm == "A" & kind == "real"]
summ <- function(d, label) d[, {
  nn <- sum(kind == "null"); xn <- sum(pos[kind == "null"]); nr <- sum(kind == "real"); xr <- sum(pos[kind == "real"])
  ci <- binom.test(xn, nn)$conf.int; ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  .(analysis = label, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn, FPR_CI_lo = ci[1], FPR_CI_hi = ci[2],
    binom_P_vs_0.05 = binom.test(xn, nn, p = 0.05)$p.value, RA_typeI_inflated = ci[1] > 0.075,
    fisher_OR = unname(ft$estimate), fisher_P = ft$p.value, RB_real_no_advantage = ft$p.value >= 0.05) }, by = module]
out <- list()
for (nm in c("fitsig_k", "fitsig_pool100")) {
  nl <- lng(fread(file.path(O, paste0("ring3_armA_", nm, "/results_m06.tsv"))), fread(file.path(O, paste0("ring3_armA_", nm, "/results_m7.tsv"))))
  d <- rbind(real, nl)
  out[[length(out) + 1]] <- summ(d[src %in% full_ids], paste0(nm, "|fullcoverage"))
  out[[length(out) + 1]] <- summ(d, paste0(nm, "|all"))
}
out <- rbindlist(out); fwrite(out, file.path(O, "ring3_armA_decision.tsv"), sep = "\t")
cat("DONE\n")
