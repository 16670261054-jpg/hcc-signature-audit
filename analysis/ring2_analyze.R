# 第 2 环汇总（预注册 N17_ring2_preregistration.md D 节；含变更 Ak）
#   R²：真实 vs 随机，每队列 Mann–Whitney U（双侧）；每签名 Phipson–Smyth 经验 P（随机 R² >= 真实）
#   残余：两个外部队列 LRT 都 <0.05 的比例，真实 vs 随机 Fisher（双侧）；随机阳性率 Clopper–Pearson 95% CI；OXPHOS 敏感性同法
#   比较：A 真实 vs A 随机（fitsig_k，主）/ A 真实 vs Ak（尺寸匹配，敏感性）/ B 真实 vs B 随机
#   另报（预注册外、仅描述）：A 组限定三队列全覆盖签名
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out/ring2")
R <- fread(file.path(O, "ring2_rows.tsv"))
full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
real_A <- R[arm == "A" & kind == "real"]; null_A <- R[arm == "A" & kind == "null"]; null_Ak <- R[arm == "Ak"]
real_B <- R[arm == "B" & kind == "real"]; null_B <- R[arm == "B" & kind == "null"]
comps <- list(A_vs_fitsig_k = list(real_A, null_A), A_vs_Ak_sizematched = list(real_A, null_Ak), B_vs_fitsig_k = list(real_B, null_B),
              Afull_vs_fitsig_k = list(real_A[src %in% full_ids], null_A[src %in% full_ids]),
              Afull_vs_Ak_sizematched = list(real_A[src %in% full_ids], null_Ak[src %in% full_ids]))

## ---- R² 与代理相关 ----
r2 <- rbindlist(lapply(names(comps), function(nm) { re <- comps[[nm]][[1]]; nu <- comps[[nm]][[2]]
  rbindlist(lapply(c("TCGA", "GSE14520", "LIRI"), function(k) { a <- re[cohort == k]; b <- nu[cohort == k]
    w <- wilcox.test(a$R2_proxy, b$R2_proxy, exact = FALSE)
    data.table(comparison = nm, cohort = k, n_real = nrow(a), n_null = nrow(b),
               R2_median_real = median(a$R2_proxy), R2_median_null = median(b$R2_proxy), MWU_P = w$p.value,
               k_model_median_real = median(a$k_model), k_model_median_null = median(b$k_model),
               rho_prolif_median_real = median(a$rho_prolif), rho_prolif_median_null = median(b$rho_prolif),
               rho_stromal_median_real = median(a$rho_stromal), rho_immune_median_real = median(a$rho_immune),
               rho_oxphos_median_real = median(a$rho_oxphos),
               absrho_prolif_ge0.5_real = mean(abs(a$rho_prolif) >= 0.5), absrho_prolif_ge0.5_null = mean(abs(b$rho_prolif) >= 0.5)) })) }))
fwrite(r2, file.path(O, "ring2_R2_summary.tsv"), sep = "\t")

## ---- 每签名经验 P（R²）----
emp <- rbindlist(lapply(names(comps)[1:3], function(nm) { re <- comps[[nm]][[1]]; nu <- comps[[nm]][[2]]
  re[, .(comparison = nm, src, cohort, R2_real = R2_proxy,
         emp_P_R2 = mapply(function(s, k, v) { x <- nu[src == s & cohort == k]$R2_proxy; (sum(x >= v) + 1) / (length(x) + 1) }, src, cohort, R2_proxy))] }))
fwrite(emp, file.path(O, "ring2_R2_empirical_per_signature.tsv"), sep = "\t")
emp_sum <- emp[, .(n = .N, frac_empP_lt0.05 = mean(emp_P_R2 < 0.05), median_empP = median(emp_P_R2)), by = .(comparison, cohort)]
fwrite(emp_sum, file.path(O, "ring2_R2_empirical_summary.tsv"), sep = "\t")

## ---- 残余预后信号 ----
both <- function(d, col) { w <- dcast(d[cohort != "TCGA"], arm + src + kind + rep ~ cohort, value.var = col)
  w[, pos := as.integer(GSE14520 < 0.05 & LIRI < 0.05)]; w }
res <- rbindlist(lapply(names(comps), function(nm) rbindlist(lapply(c("lrt_p", "lrt_p_oxphos"), function(col) {
  a <- both(comps[[nm]][[1]], col); b <- both(comps[[nm]][[2]], col)
  xr <- sum(a$pos, na.rm = TRUE); nr <- sum(!is.na(a$pos)); xn <- sum(b$pos, na.rm = TRUE); nn <- sum(!is.na(b$pos))
  ci <- binom.test(xn, nn)$conf.int; ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(comparison = nm, base_model = ifelse(col == "lrt_p", "stage+prolif+stromal+immune+age+sex", "+oxphos"),
             n_real = nr, pos_real = xr, rate_real = xr / nr, n_null = nn, pos_null = xn, rate_null = xn / nn,
             null_CP_lo = ci[1], null_CP_hi = ci[2], fisher_OR = unname(ft$estimate), fisher_P = ft$p.value,
             any_ext_pos_real = mean(a$GSE14520 < 0.05 | a$LIRI < 0.05, na.rm = TRUE),
             any_ext_pos_null = mean(b$GSE14520 < 0.05 | b$LIRI < 0.05, na.rm = TRUE)) }))))
fwrite(res, file.path(O, "ring2_residual_summary.tsv"), sep = "\t")

## ---- 调整前后 HR 衰减（描述）----
hr <- rbindlist(lapply(names(comps)[1:3], function(nm) { re <- comps[[nm]][[1]][cohort != "TCGA"]; nu <- comps[[nm]][[2]][cohort != "TCGA"]
  f <- function(d) d[, .(median_abslogHR_unadj = median(abs(log(hr_unadj))), median_abslogHR_adj = median(abs(log(hr_adj))),
                         median_attenuation = median(1 - abs(log(hr_adj)) / abs(log(hr_unadj)))), by = cohort]
  rbind(cbind(comparison = nm, group = "real", f(re)), cbind(comparison = nm, group = "null", f(nu))) }))
fwrite(hr, file.path(O, "ring2_HR_attenuation.tsv"), sep = "\t")
cat("DONE\n")
