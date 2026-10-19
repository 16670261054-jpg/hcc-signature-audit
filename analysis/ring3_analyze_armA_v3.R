# v3: adds the size-matched control fitsig_pool100_k; output written to ring3_armA_decision_v3_<weight>.tsv
# Arm A calls against the predefined controls: real published signatures (from ring3_full) vs template-pipeline random controls
suppressPackageStartupMessages({ library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out")
cov <- fread(file.path(O, "ring3_full/published_coverage.tsv")); full_ids <- cov[full == TRUE]$sig_id
mods <- c("M0_pos", "M1_pos", "M2_pos", "M4_pos", "M5a_pos", "M5b_pos", "M6_pos")
# M0 requires direction: C > 0.5 and two-sided Noether P < 0.05 in both external cohorts (survcomp default is two-sided)
fixM0 <- function(d) { if (all(c("M0_c_GSE", "M0_c_LIRI", "M0_p_GSE", "M0_p_LIRI") %in% names(d)))
  d[, M0_pos := as.integer(M0_p_GSE < 0.05 & M0_p_LIRI < 0.05 & M0_c_GSE > 0.5 & M0_c_LIRI > 0.5)]; d }
lng <- function(R, R7) { R <- fixM0(copy(R)); x <- rbind(melt(R, id.vars = c("arm", "src", "kind", "rep", "k"), measure.vars = mods, variable.name = "module", value.name = "pos"),
                                    R7[, .(arm, src, kind, rep, k, module = "M7_pos", pos = M7_pos)]); x[, module := as.character(module)]; x }
WS <- Sys.getenv("R3_WEIGHT", "rescaled")   # rescaled = coef x sd (main); z = original z-scale weights (sensitivity)
m06 <- if (WS == "rescaled") fread(file.path(O, "ring3_armA_realrescale/results_m06_real_rescaled.tsv")) else fread(file.path(O, "ring3_full/results_m06.tsv"))
real <- lng(m06, fread(file.path(O, "ring3_full/results_m7.tsv")))[arm == "A" & kind == "real"]
summ <- function(d, label) d[, {
  nn <- sum(kind == "null"); xn <- sum(pos[kind == "null"]); nr <- sum(kind == "real"); xr <- sum(pos[kind == "real"])
  ci <- binom.test(xn, nn)$conf.int; ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  .(analysis = label, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn, FPR_CI_lo = ci[1], FPR_CI_hi = ci[2],
    binom_P_vs_0.05 = binom.test(xn, nn, p = 0.05)$p.value, RA_typeI_inflated = ci[1] > 0.075,
    fisher_OR = unname(ft$estimate), fisher_P = ft$p.value, RB_real_no_advantage = ft$p.value >= 0.05) }, by = module]
out <- list()
for (nm in c("fitsig_k", "fitsig_pool100", "fitsig_pool100_k")) {
  nl <- lng(fread(file.path(O, paste0("ring3_armA_", nm, "/results_m06.tsv"))), fread(file.path(O, paste0("ring3_armA_", nm, "/results_m7.tsv"))))
  d <- rbind(real, nl)
  out[[length(out) + 1]] <- summ(d[src %in% full_ids], paste0(nm, "|fullcoverage"))
  out[[length(out) + 1]] <- summ(d, paste0(nm, "|all"))
}
out <- rbindlist(out); fwrite(out, file.path(O, paste0("ring3_armA_decision_v3_", WS, ".tsv")), sep = "\t")
cat("DONE\n")
