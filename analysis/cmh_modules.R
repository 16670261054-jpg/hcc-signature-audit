# 16 modules, published vs size-matched random: Cochran-Mantel-Haenszel test stratified by matched group (each published signature + its random signatures), Holm correction
# Positive rates and CIs are unchanged; only the between-group test is replaced
suppressPackageStartupMessages(library(data.table))
O <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_out")
fixM0 <- function(d) { if ("M0_c_GSE" %in% names(d)) d[, M0_pos := as.integer(M0_p_GSE < 0.05 & M0_p_LIRI < 0.05 & M0_c_GSE > 0.5 & M0_c_LIRI > 0.5)]; d }
pub06 <- fixM0(fread(file.path(O, "ring3_armA_realrescale/results_m06_real_rescaled.tsv"))[arm == "A" & kind == "real"])
nul06 <- fixM0(rbindlist(lapply(list.files(file.path(O, "ring3_armA_fitsig_pool100_k/chunks"), "^m06_", full.names = TRUE), fread))[arm == "A" & kind == "null"])
pub7 <- fread(file.path(O, "ring3_full/results_m7.tsv"))[arm == "A" & kind == "real"]
nul7 <- fread(file.path(O, "ring3_armA_fitsig_pool100_k/results_m7.tsv"))[arm == "A" & kind == "null"]
hi <- function(f, cols) { d <- fread(file.path(O, f)); list(pub = d[arm == "A" & kind == "real"], nul = d[arm == "Ak"]) }
m3 <- hi("hiend_m3_m12/m3_m12_per_signature.tsv"); m8 <- hi("hiend_m8/m8_per_signature.tsv"); m10 <- hi("hiend_m10/m10_per_signature.tsv"); m13 <- hi("hiend_m13_m15/m13_m15_per_signature.tsv")
spec <- list(M0 = list(pub06, nul06, "M0_pos"), M1 = list(pub06, nul06, "M1_pos"), M2 = list(pub06, nul06, "M2_pos"), M4 = list(pub06, nul06, "M4_pos"),
  M5a = list(pub06, nul06, "M5a_pos"), M5b = list(pub06, nul06, "M5b_pos"), M6 = list(pub06, nul06, "M6_pos"), M7 = list(pub7, nul7, "M7_pos"),
  M3 = list(m3$pub, m3$nul, "M3_pos"), M12 = list(m3$pub, m3$nul, "M12_pos"), M8a = list(m8$pub, m8$nul, "M8a_pos"), M8b = list(m8$pub, m8$nul, "M8b_pos"),
  M10 = list(m10$pub, m10$nul, "M10_pos"), M13 = list(m13$pub, m13$nul, "M13_pos"), M14 = list(m13$pub, m13$nul, "M14_pos"), M15 = list(m13$pub, m13$nul, "DCA_pos"))
res <- rbindlist(lapply(names(spec), function(m) { sp <- spec[[m]]; v <- sp[[3]]
  P <- sp[[1]][!is.na(get(v)), .(src, xr = get(v))]; N <- sp[[2]][!is.na(get(v)), .(nn = .N, xn = sum(get(v))), by = src]
  d <- merge(P, N, by = "src")
  a <- array(0, c(2, 2, nrow(d))); for (i in seq_len(nrow(d))) a[, , i] <- matrix(c(d$xr[i], 1 - d$xr[i], d$xn[i], d$nn[i] - d$xn[i]), 2)
  inf <- apply(a, 3, function(x) all(rowSums(x) > 0) & all(colSums(x) > 0))
  ct <- if (sum(inf) >= 1) tryCatch(mantelhaen.test(a[, , inf, drop = FALSE], exact = FALSE), error = function(e) NULL) else NULL
  data.table(module = m, n_pub = nrow(d), n_rand = sum(d$nn), rate_pub = mean(d$xr), rate_rand = sum(d$xn) / sum(d$nn),
             n_informative_strata = sum(inf), CMH_OR = if (is.null(ct)) NA_real_ else unname(ct$estimate), CMH_P = if (is.null(ct)) 1 else ct$p.value) }))
res[, CMH_P_holm := p.adjust(CMH_P, "holm")]
fwrite(res, file.path(O, "cmh_modules.tsv"), sep = "\t")
F4 <- fread(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N28_fig/source_data/Fig4a_module_rates.tsv"))
chk <- merge(res[, .(module, rate_pub, rate_rand)], F4[, .(module, rate_real, rate_null)], by = "module")
stopifnot(nrow(chk) == 16, max(abs(chk$rate_pub - chk$rate_real)) < 1e-9)
fwrite(chk[, .(module, d_pub = rate_pub - rate_real, d_rand = rate_rand - rate_null)], file.path(O, "cmh_modules_consistency.tsv"), sep = "\t")
print(res[order(CMH_P)]); print(chk[, .(module, d_rand = round(rate_rand - rate_null, 6))])
