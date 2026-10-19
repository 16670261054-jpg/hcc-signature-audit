# Proxy analysis summary (analysis plan fixed before running; including the Ak amendment)
#   R²: real vs random, Mann-Whitney U per cohort (two-sided); per-signature Phipson-Smyth empirical P (random R² >= real)
#   Residual: proportion with LRT < 0.05 in both external cohorts, real vs random Fisher (two-sided); random positive rate Clopper-Pearson 95% CI; OXPHOS sensitivity likewise
#   Comparisons: A real vs A random (fitsig_k) / A real vs Ak (size-matched) / B real vs B random
#   Also reported (descriptive): arm A restricted to signatures fully covered in all three cohorts
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out/ring2")
R <- fread(file.path(O, "ring2_rows.tsv"))
full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
real_A <- R[arm == "A" & kind == "real"]; null_A <- R[arm == "A" & kind == "null"]; null_Ak <- R[arm == "Ak"]
real_B <- R[arm == "B" & kind == "real"]; null_B <- R[arm == "B" & kind == "null"]
comps <- list(A_vs_fitsig_k = list(real_A, null_A), A_vs_Ak_sizematched = list(real_A, null_Ak), B_vs_fitsig_k = list(real_B, null_B),
              Afull_vs_fitsig_k = list(real_A[src %in% full_ids], null_A[src %in% full_ids]),
              Afull_vs_Ak_sizematched = list(real_A[src %in% full_ids], null_Ak[src %in% full_ids]))

## ---- R² and proxy correlations ----
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

## ---- Per-signature empirical P (R²) ----
emp <- rbindlist(lapply(names(comps)[1:3], function(nm) { re <- comps[[nm]][[1]]; nu <- comps[[nm]][[2]]
  re[, .(comparison = nm, src, cohort, R2_real = R2_proxy,
         emp_P_R2 = mapply(function(s, k, v) { x <- nu[src == s & cohort == k]$R2_proxy; (sum(x >= v) + 1) / (length(x) + 1) }, src, cohort, R2_proxy))] }))
fwrite(emp, file.path(O, "ring2_R2_empirical_per_signature.tsv"), sep = "\t")
emp_sum <- emp[, .(n = .N, frac_empP_lt0.05 = mean(emp_P_R2 < 0.05), median_empP = median(emp_P_R2)), by = .(comparison, cohort)]
fwrite(emp_sum, file.path(O, "ring2_R2_empirical_summary.tsv"), sep = "\t")

## ---- Residual prognostic signal ----
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

## ---- HR attenuation after adjustment (descriptive) ----
hr <- rbindlist(lapply(names(comps)[1:3], function(nm) { re <- comps[[nm]][[1]][cohort != "TCGA"]; nu <- comps[[nm]][[2]][cohort != "TCGA"]
  f <- function(d) d[, .(median_abslogHR_unadj = median(abs(log(hr_unadj))), median_abslogHR_adj = median(abs(log(hr_adj))),
                         median_attenuation = median(1 - abs(log(hr_adj)) / abs(log(hr_unadj)))), by = cohort]
  rbind(cbind(comparison = nm, group = "real", f(re)), cbind(comparison = nm, group = "null", f(nu))) }))
fwrite(hr, file.path(O, "ring2_HR_attenuation.tsv"), sep = "\t")
cat("DONE\n")
