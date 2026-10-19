# Table S3 master table of numbers: every number quoted in the text, legends and abstract is taken from here
# One row per citable quantity: id, claim, figure panel, meaning, raw value, displayed value, source file, extraction rule
# All values are read from disk at run time; displayed values are rounded only here
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
rd <- function(f, dir = SD) fread(file.path(dir, f), sep = "\t", quote = "")
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
f3 <- function(x) sprintf("%.2f", x)
pfmt <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), sprintf("%.2g", p))   # two significant digits in the text
N <- list(); add <- function(id, claim, fig, what, value, display, file, rule) {
  stopifnot(length(value) == 1, !is.na(value))
  N[[length(N) + 1]] <<- data.table(id, claim, fig, what, value = signif(value, 8), display, source = file, rule) }

## ---- Corpus and cohorts ----
CM <- rd("cohort_meta.tsv", file.path(OUT, "robust_v2"))
for (cn in CM$cohort) { add(paste0("n_", cn), "Data", "Fig1a", paste(cn, "patients"), CM[cohort == cn]$n, CM[cohort == cn]$n, "N11_out/robust_v2/cohort_meta.tsv", "n")
  add(paste0("ev_", cn), "Data", "Fig1a", paste(cn, "deaths"), CM[cohort == cn]$events, CM[cohort == cn]$events, "N11_out/robust_v2/cohort_meta.tsv", "events") }
add("n_bulk_total", "Data", "Abstract", "Patients in three bulk cohorts", sum(CM$n), sum(CM$n), "N11_out/robust_v2/cohort_meta.tsv", "sum(n)")
LS <- rd("Fig4c_landscape_per_signature_modules.tsv")
add("n_sig_scorable", "Data", "Fig1a", "Published signatures scorable", nrow(LS), nrow(LS), "Fig4c_landscape_per_signature_modules.tsv", "nrow")
MF <- rd("Fig1b_module_frequency.tsv")
add("n_fulltext", "Data", "Fig1b", "Full-text papers audited for downstream analyses", MF$n_fulltext[1], MF$n_fulltext[1], "Fig1b_module_frequency.tsv", "n_fulltext")
for (m in c("timeROC/ROC", "KM", "CIBERSORT/immune", "nomogram", "qPCR/IHC wet", "single-cell"))
  add(paste0("corpus_", gsub("[^A-Za-z]", "", m)), "Data", "Fig1b", paste("Share of papers using", m), MF[module == m]$frac,
      sprintf("%d/%d (%s)", MF[module == m]$n_papers, MF[module == m]$n_fulltext, pct(MF[module == m]$frac, 0)), "Fig1b_module_frequency.tsv", paste0("module == '", m, "'"))
UM <- fread(file.path(DAT, "sc_umap.tsv.gz"), select = c("patient"))
add("sc_cells", "Data", "Fig6a", "scRNA-seq cells after QC", nrow(UM), format(nrow(UM), big.mark = ","), "N28_fig/data/sc_umap.tsv.gz", "nrow")
add("sc_patients", "Data", "Fig6a", "scRNA-seq patients", uniqueN(UM$patient), uniqueN(UM$patient), "N28_fig/data/sc_umap.tsv.gz", "uniqueN(patient)")

## ---- Design parameters (parsed directly from the analysis script source to avoid manual copying) ----
RP <- readLines(file.path(B, "analysis/robust_proxy.R"), encoding = "UTF-8")
gnum <- function(key) { m <- regmatches(RP, regexec(paste0(key, "\", \"([0-9]+)\""), RP)); m <- Filter(length, m); stopifnot(length(m) >= 1); as.numeric(m[[1]][2]) }
add("D_n_splits", "Design", "Fig2a", "TCGA train/test splits", gnum('N11_S'), gnum('N11_S'), "N11_scripts/robust_proxy.R", "N11_S default")
add("D_null_inner", "Design", "Fig2a", "Random gene sets per split (internal null)", gnum('N11_NR'), gnum('N11_NR'), "N11_scripts/robust_proxy.R", "N11_NR default")
add("D_null_ext", "Design", "Fig2b-c", "Random gene sets (external null)", gnum('N11_NE'), gnum('N11_NE'), "N11_scripts/robust_proxy.R", "N11_NE default")
stopifnot(any(grepl("round(nT * 2 / 3)", RP, fixed = TRUE)))
add("D_train_frac", "Design", "Fig2a", "Training fraction per split", 2 / 3, "two-thirds", "N11_scripts/robust_proxy.R", "round(nT * 2 / 3)")
add("D_universe", "Design", "Fig2", "Genes in random-sampling pool", CM$n_genes_universe[1], format(CM$n_genes_universe[1], big.mark = ","), "N11_out/robust_v2/cohort_meta.tsv", "n_genes_universe")
stopifnot(nrow(rd("Fig2a_inner_splits.tsv")) == 8 * gnum('N11_S'))
## ---- C1 ----
IS <- rd("Fig2a_inner_splits_summary.tsv")
add("C1_inner_empP_min", "C1", "Fig2a", "Median internal empirical P, lowest of 8 sets", min(IS$median_emp_p), f3(min(IS$median_emp_p)), "Fig2a_inner_splits_summary.tsv", "min(median_emp_p)")
add("C1_inner_empP_max", "C1", "Fig2a", "Median internal empirical P, highest of 8 sets", max(IS$median_emp_p), f3(max(IS$median_emp_p)), "Fig2a_inner_splits_summary.tsv", "max(median_emp_p)")
add("C1_n_sets", "C1", "Fig2", "Cell-death gene sets tested", nrow(IS), nrow(IS), "Fig2a_inner_splits_summary.tsv", "nrow")
EX <- rd("Fig2b_external.tsv")
add("C1_ext_above_p95", "C1", "Fig2b", "Set x cohort pairs with external C above random 95th percentile", sum(EX$c_ext > EX$null_p95), sprintf("%d/%d", sum(EX$c_ext > EX$null_p95), nrow(EX)), "Fig2b_external.tsv", "sum(c_ext > null_p95)")
PR <- rd("Fig2b_proliferation_reference.tsv")
for (cn in PR$cohort) add(paste0("C1_prolif_C_", cn), "C1", "Fig2b", paste("Proliferation score alone, Harrell's C,", cn), PR[cohort == cn]$c, f3(PR[cohort == cn]$c), "Fig2b_proliferation_reference.tsv", paste0("cohort == '", cn, "'"))
for (cn in PR$cohort) add(paste0("C1_sets_below_prolif_", cn), "C1", "Fig2b", paste("Cell-death sets with external C below proliferation score alone,", cn),
  sum(EX[cohort == cn]$c_ext < PR[cohort == cn]$c), sprintf("%d/%d", sum(EX[cohort == cn]$c_ext < PR[cohort == cn]$c), nrow(EX[cohort == cn])), "Fig2b_external.tsv + Fig2b_proliferation_reference.tsv", "c_ext < prolif c")
AN <- rd("Fig2c_adjusted_null.tsv")
add("C1_adjnull_frac_min", "C1", "Fig2c", "Random sets still P<0.05 after proxy adjustment, lowest cohort x set", min(AN$null_frac_Padj_lt05), pct(min(AN$null_frac_Padj_lt05)), "Fig2c_adjusted_null.tsv", "min(null_frac_Padj_lt05)")
add("C1_adjnull_frac_max", "C1", "Fig2c", "Random sets still P<0.05 after proxy adjustment, highest cohort x set", max(AN$null_frac_Padj_lt05), pct(max(AN$null_frac_Padj_lt05)), "Fig2c_adjusted_null.tsv", "max(null_frac_Padj_lt05)")
add("C1_adj_specific", "C1", "Fig2c", "Set x cohort pairs with residual signal beyond random (emp P<0.05)", sum(AN$emp_p_absz < 0.05), sprintf("%d/%d", sum(AN$emp_p_absz < 0.05), nrow(AN)), "Fig2c_adjusted_null.tsv", "sum(emp_p_absz < 0.05)")

## ---- C2 ----
ST <- rd("Fig3ab_stats.tsv")
for (cn in c("TCGA", "GSE14520", "LIRI")) for (mt in c("R2_proxy", "rho_prolif")) {
  r <- ST[metric == mt & cohort == cn & grp2 == "Size-matched random"]
  add(sprintf("C2_%s_pub_%s", mt, cn), "C2", "Fig3a-b", paste("Published median", mt, cn), r$median_pub, f3(r$median_pub), "Fig3ab_stats.tsv", "median_pub")
  add(sprintf("C2_%s_Ak_%s", mt, cn), "C2", "Fig3a-b", paste("Size-matched random median", mt, cn), r$median_null, f3(r$median_null), "Fig3ab_stats.tsv", "median_null")
  add(sprintf("C2_%s_P_%s", mt, cn), "C2", "Fig3a-b", paste("MWU P published vs size-matched,", mt, cn), r$P, pfmt(r$P), "Fig3ab_stats.tsv", "P") }
for (cn in c("TCGA", "GSE14520", "LIRI")) for (mt in c("R2_proxy", "rho_prolif")) {
  r <- ST[metric == mt & cohort == cn & grp2 == "Input-matched random"]
  add(sprintf("C2_%s_In_%s", mt, cn), "C2", "Fig3a-b", paste("Input-matched random median", mt, cn), r$median_null, f3(r$median_null), "Fig3ab_stats.tsv", "median_null, Input-matched")
  add(sprintf("C2_%s_InP_%s", mt, cn), "C2", "Fig3a-b", paste("MWU P published vs input-matched,", mt, cn), r$P, pfmt(r$P), "Fig3ab_stats.tsv", "P, Input-matched") }
RS <- rd("Fig3c_residual.tsv")
for (g in names(COL_GRP)[1:3]) for (md in unique(RS$model)) { r <- RS[grp == g & model == md]; tag <- fifelse(grepl("OXPHOS", md), "oxphos", "base")
  add(sprintf("C2_resid_%s_%s", gsub("[^A-Za-z]", "", g), tag), "C2", "Fig3c", paste("Residual signal in both external cohorts,", g, md), r$rate,
      sprintf("%s (%d/%d)", pct(r$rate), r$pos, r$n), "Fig3c_residual.tsv", "pos/n")
  if (!is.na(r$fisher_P_vs_published)) add(sprintf("C2_resid_P_%s_%s", gsub("[^A-Za-z]", "", g), tag), "C2", "Fig3c", paste("Fisher P vs published,", g, md),
      r$fisher_P_vs_published, pfmt(r$fisher_P_vs_published), "Fig3c_residual.tsv", "fisher_P_vs_published") }
HM <- rd("Fig3d_HR_medians.tsv")
for (cn in unique(HM$cohort_l)) { r <- HM[cohort_l == cn & grp == "Published"]
  add(paste0("C2_absloghr_unadj_", gsub("-", "", cn)), "C2", "Fig3d", paste("Published median |log HR| unadjusted,", cn), r$u, f3(r$u), "Fig3d_HR_medians.tsv", "u")
  add(paste0("C2_absloghr_adj_", gsub("-", "", cn)), "C2", "Fig3d", paste("Published median |log HR| adjusted,", cn), r$a, f3(r$a), "Fig3d_HR_medians.tsv", "a")
}
## Single definition of attenuation: median of per-signature attenuation (proxy-analysis values on disk), not the ratio of medians (they differ: GSE 20% vs 24%)
HA <- rd("ring2_HR_attenuation.tsv", file.path(OUT, "ring2"))[comparison == "A_vs_Ak_sizematched"]
for (cn in c("GSE14520", "LIRI")) for (g in c("real", "null")) { r <- HA[cohort == cn & group == g]
  add(sprintf("C2_attenuation_%s_%s", g, cn), "C2", "Fig3d", sprintf("Median per-signature attenuation of |log HR|, %s, %s", fifelse(g == "real", "published", "size-matched random"), cn),
      r$median_attenuation, pct(r$median_attenuation, 0), "N11_out/ring2/ring2_HR_attenuation.tsv", "median_attenuation") }

## ---- C3 ----
M4 <- rd("Fig4a_module_rates.tsv")
for (i in seq_len(nrow(M4))) { r <- M4[i]
  add(paste0("C3_", r$module, "_rand"), "C3", "Fig4a", paste("Size-matched random positive rate,", r$label), r$rate_null, sprintf("%s (95%% CI %s-%s)", pct(r$rate_null), pct(r$lo), pct(r$hi)), "Fig4a_module_rates.tsv", "rate_null; lo; hi")
  add(paste0("C3_", r$module, "_pub"), "C3", "Fig4a", paste("Published positive rate,", r$label), r$rate_real, pct(r$rate_real), "Fig4a_module_rates.tsv", "rate_real")
  add(paste0("C3_", r$module, "_P"), "C3", "Fig4a", paste("Fisher P published vs size-matched,", r$label), r$P, pfmt(r$P), "Fig4a_module_rates.tsv", "P") }
add("C3_n_modules_over_bradley", "C3", "Fig4a", "Modules whose random positive-rate CI lower bound exceeds 7.5%", sum(M4$lo > BRADLEY), sprintf("%d/%d", sum(M4$lo > BRADLEY), nrow(M4)), "Fig4a_module_rates.tsv", "sum(lo > 0.075)")
add("C3_n_modules_pub_higher", "C3", "Fig4a", "Modules where published > random (P<0.05)", sum(M4$P < 0.05 & M4$rate_real > M4$rate_null), sum(M4$P < 0.05 & M4$rate_real > M4$rate_null), "Fig4a_module_rates.tsv", "P<0.05 & rate_real>rate_null")
add("C3_n_modules_pub_lower", "C3", "Fig4a", "Modules where published < random (P<0.05)", sum(M4$P < 0.05 & M4$rate_real < M4$rate_null), sum(M4$P < 0.05 & M4$rate_real < M4$rate_null), "Fig4a_module_rates.tsv", "P<0.05 & rate_real<rate_null")
## Positive rates of input-matched random signatures in the 16 modules (non-specificity does not depend on size selection)
r3i <- fread(file.path(OUT, "ring3_armA_decision_v3_rescaled.tsv"))[analysis == "fitsig_k|all", .(module = sub("_pos", "", module), lo = FPR_CI_lo, rate = FPR_null)]
pk <- function(f, mods) fread(file.path(OUT, f))[comparison == "A_vs_fitsig_k" & module %in% mods, .(module = sub("_pos$", "", module), lo = FPR_CP_lo, rate = FPR_null)]
IM <- rbind(r3i, pk("hiend_m3_m12/m3_m12_decision.tsv", c("M3_pos", "M12_pos")), pk("hiend_m8/m8_decision.tsv", c("M8a_pos", "M8b_pos")),
            pk("hiend_m13_m15/m13_m15_decision.tsv", c("M13_pos", "M14_pos", "DCA_pos")),
            fread(file.path(OUT, "hiend_m10/m10_decision.tsv"))[comparison == "A_vs_fitsig_k" & filter == "all", .(module = "M10", lo = FPR_CP_lo, rate = FPR_null)])
stopifnot(nrow(IM) == 16)
add("C3_n_modules_over_bradley_input", "C3", "Fig4a", "Modules whose input-matched random positive-rate CI lower bound exceeds 7.5%", sum(IM$lo > BRADLEY), sprintf("%d/%d", sum(IM$lo > BRADLEY), nrow(IM)), "decision tables, A_vs_fitsig_k", "sum(lo > 0.075)")
add("C3_M10_rand_input", "C3", "Fig4a", "Input-matched random positive rate, Geneformer", IM[module == "M10"]$rate, pct(IM[module == "M10"]$rate), "hiend_m10/m10_decision.tsv", "A_vs_fitsig_k")
write_sd(IM, "Fig4a_input_matched_rates")
## Proxies recomputed after removing signature genes (n4_loo_proxy.R)
LOO <- rd("loo_summary.tsv", file.path(OUT, "n4_loo")); LOC <- rd("loo_consistency_check.tsv", file.path(OUT, "n4_loo"))
stopifnot(LOC$n_matched == LOC$n_rows, LOC$max_abs_diff < 1e-6)
add("N4_robust_n", "Sens", "Supp", "Cohort x metric combinations where order and P<0.05 hold after leave-signature-genes-out", sum(LOO$robust), sprintf("%d/%d", sum(LOO$robust), nrow(LOO)), "N11_out/n4_loo/loo_summary.tsv", "sum(robust)")
for (cn in c("TCGA", "GSE14520", "LIRI")) { r <- LOO[cohort == cn & metric == "R2_proxy_loo"]
  add(paste0("N4_R2_pub_", cn), "Sens", "Supp", paste("Published median proxy R2 after leave-out,", cn), r$pub, f3(r$pub), "N11_out/n4_loo/loo_summary.tsv", "pub")
  add(paste0("N4_R2_size_", cn), "Sens", "Supp", paste("Size-matched median proxy R2 after leave-out,", cn), r$sizematched, f3(r$sizematched), "N11_out/n4_loo/loo_summary.tsv", "sizematched") }
## Selected simulation (ring4v2_selected_sim.R) + published signatures through the gates (n2_published_gates.R)
S1 <- rd("sim_power.tsv", file.path(OUT, "ring4v2")); P2 <- rd("published_gates_summary.tsv", file.path(OUT, "n2_gates"))
sv <- function(sc, bt, cr) S1[scen == sc & abs(beta - bt) < 1e-9 & criterion == cr]
for (cr in c("T0", "FW", "FWp", "T1p", "T2p")) { r <- sv("selected", 0, cr)
  add(paste0("N1_null_", cr), "C6", "Fig7b", paste("Selected-null simulation positive rate,", cr), r$rate, sprintf("%s (95%% CI %s-%s)", pct(r$rate), pct(r$CP_lo), pct(r$CP_hi)), "N11_out/ring4v2/sim_power.tsv", paste("selected, beta 0,", cr)) }
for (cr in c("FWp", "T2p", "FW")) { r <- sv("misspec_oxphos", 0, cr)
  add(paste0("N1_misspec_", cr), "C6", "Fig7b", paste("OXPHOS-misspecified null positive rate,", cr), r$rate, sprintf("%s (95%% CI %s-%s)", pct(r$rate), pct(r$CP_lo), pct(r$CP_hi)), "N11_out/ring4v2/sim_power.tsv", paste("misspec, beta 0,", cr)) }
for (bt in c(0.2, 0.4, 0.6)) for (cr in c("T0", "FW", "FWp")) { r <- sv("selected", bt, cr)
  add(sprintf("N1_power_%s_b%02d", cr, round(10 * bt)), "C6", "Fig7b", sprintf("Selected simulation positive rate, %s, HR/SD %.2f", cr, exp(bt)), r$rate, pct(r$rate), "N11_out/ring4v2/sim_power.tsv", sprintf("selected, beta %.1f, %s", bt, cr)) }
for (cr in P2$criterion) { r <- P2[criterion == cr]
  add(paste0("N2_", cr), "C6", "Fig7c", paste("Published signatures passing,", cr), r$rate, sprintf("%d/%d (%s)", r$positive, r$n, pct(r$rate)), "N11_out/n2_gates/published_gates_summary.tsv", cr) }
fp <- fisher.test(matrix(c(P2[criterion == "both_pooled"]$positive, P2[criterion == "both_pooled"]$n - P2[criterion == "both_pooled"]$positive,
                           sv("selected", 0, "FWp")$positive, sv("selected", 0, "FWp")$n - sv("selected", 0, "FWp")$positive), 2))$p.value
add("N2_vs_null_P", "C6", "Fig7c", "Fisher P, published passing both pooled gates vs selected-null simulation", fp, pfmt(fp), "n2_gates + ring4v2", "fisher")
FN <- rd("Fig4b_funnel_per_signature.tsv")
add("C3_funnel_M1", "C3", "Fig4b", "Published passing independent prognostic factor", sum(FN$M1_pos), sprintf("%d/%d", sum(FN$M1_pos), nrow(FN)), "Fig4b_funnel_per_signature.tsv", "sum(M1_pos)")
add("C3_funnel_M0", "C3", "Fig4b", "Published passing external validation", sum(FN$M0_pos), sprintf("%d/%d", sum(FN$M0_pos), nrow(FN)), "Fig4b_funnel_per_signature.tsv", "sum(M0_pos)")
add("C3_funnel_resid", "C3", "Fig4b", "Published surviving proxy adjustment", sum(FN$resid), sprintf("%d/%d", sum(FN$resid), nrow(FN)), "Fig4b_funnel_per_signature.tsv", "sum(resid)")
add("C3_funnel_all3", "C3", "Fig4b", "Published passing all three", sum(FN$M0_pos & FN$M1_pos & FN$resid), sprintf("%d/%d", sum(FN$M0_pos & FN$M1_pos & FN$resid), nrow(FN)), "Fig4b_funnel_per_signature.tsv", "all three")
add("C3_landscape_median_npos", "C3", "Fig4c", "Median modules positive per published signature (of 16)", median(LS$n_pos), median(LS$n_pos), "Fig4c_landscape_per_signature_modules.tsv", "median(n_pos)")
KM <- rd("Fig5a_km_logrank.tsv")
for (cn in KM$cohort) { r <- KM[cohort == cn]
  add(paste0("C3_showcase_HR_", cn), "C3", "Fig5a", paste("Random signature, high vs low HR,", cn), r$HR_high_vs_low, sprintf("%.2f (%.2f-%.2f)", r$HR_high_vs_low, r$HR_lo, r$HR_hi), "Fig5a_km_logrank.tsv", "HR")
  add(paste0("C3_showcase_P_", cn), "C3", "Fig5a", paste("Random signature, log-rank P,", cn), r$logrank_P, pfmt(r$logrank_P), "Fig5a_km_logrank.tsv", "logrank_P") }
AU <- rd("Fig5b_timeroc_auc.tsv")
for (y in AU$year) add(paste0("C3_showcase_AUC_", y, "y"), "C3", "Fig5b", paste0("Random signature, ", y, "-year AUC (TCGA)"), AU[year == y]$AUC, f3(AU[year == y]$AUC), "Fig5b_timeroc_auc.tsv", "AUC")
MC <- rd("Fig5c_multicox.tsv")[term == "Risk"]
add("C3_showcase_multiHR", "C3", "Fig5c", "Random signature, multivariable HR per SD (TCGA)", MC$HR, sprintf("%.2f (%.2f-%.2f)", MC$HR, MC$lo, MC$hi), "Fig5c_multicox.tsv", "term == 'Risk'")
add("C3_showcase_multiP", "C3", "Fig5c", "Random signature, multivariable P (TCGA)", MC$P, pfmt(MC$P), "Fig5c_multicox.tsv", "term == 'Risk'")
## ---- Fig S5: position of the example random signature among 3,700 size-matched random signatures (check_fig5_showcase_typicality.R) ----
TY <- rd("FigS5_typicality_summary.tsv"); ty <- function(m) TY[metric == m]
stopifnot(abs(ty("HR_TCGA")$example - KM[cohort == "TCGA"]$HR_high_vs_low) < 1e-10,
          abs(ty("multiHR")$example - MC$HR) < 1e-10)   # same signature and definitions as Fig 5
EC <- rd("Fig6c_example_choice.tsv"); stopifnot(EC$example_src == "S004")
add("C3_showcase_n_candidates", "C3", "Fig5", "Candidates meeting the example-selection rule", EC$n_candidates, EC$n_candidates, "Fig6c_example_choice.tsv", "n_candidates")
add("C3_showcase_k_median", "C3", "Fig5", "Published median genes per signature used in example rule", EC$k_median, EC$k_median, "Fig6c_example_choice.tsv", "k_median")
add("C3_showcase_n_random", "C3", "FigS5", "Size-matched random signatures compared with the example", ty("HR_TCGA")$n, format(ty("HR_TCGA")$n, big.mark = ","), "FigS5_typicality_summary.tsv", "n")
for (m in c("HR_TCGA", "HR_GSE14520", "HR_LIRI", "C_TCGA", "C_GSE14520", "C_LIRI", "multiHR"))
  add(paste0("C3_showcase_pctl_", m), "C3", "FigS5", paste("Example random signature, percentile among size-matched random,", m), ty(m)$pct_rank, sprintf("%.0f", 100 * ty(m)$pct_rank), "FigS5_typicality_summary.tsv", "pct_rank")
for (m in c("HR_TCGA", "HR_GSE14520", "HR_LIRI", "multiHR"))
  add(paste0("C3_showcase_median_", m), "C3", "FigS5", paste("Size-matched random signatures, median,", m), ty(m)$median_all, f3(ty(m)$median_all), "FigS5_typicality_summary.tsv", "median_all")
for (m in c("frac_TCGA_P05", "frac_twoext_P05", "frac_all3_P05", "frac_multiP05"))
  add(paste0("C3_showcase_", m), "C3", "FigS5", paste("Size-matched random signatures,", m), ty(m)$median_all, pct(ty(m)$median_all), "FigS5_typicality_summary.tsv", "median_all (= fraction)")

## ---- Single-cell level ----
F8 <- rd("Fig6f_M8_rates.tsv")
for (i in seq_len(nrow(F8))) { r <- F8[i]
  add(sprintf("C4_%s_%s", sub("_pos", "", r$module), gsub("[^A-Za-z]", "", r$grp)), "C4", "Fig6f", paste(r$mod_l, r$grp), r$rate, pct(r$rate), "Fig6f_M8_rates.tsv", "rate") }
DE <- rd("Fig6g_geneformer_deciles.tsv")
add("C4_gf_n_sig", "C4", "Fig6g", "Geneformer genes FDR<0.05", sum(DE$n_sig), sprintf("%d/%s", sum(DE$n_sig), format(sum(DE$n), big.mark = ",")), "Fig6g_geneformer_deciles.tsv", "sum(n_sig)/sum(n)")
add("C4_gf_top_decile_share", "C4", "Fig6g", "Share of significant genes in top abundance decile", DE[NDet_decile == 10]$share, pct(DE[NDet_decile == 10]$share, 0), "Fig6g_geneformer_deciles.tsv", "share, decile 10")
add("C4_gf_zero_deciles", "C4", "Fig6g", "Abundance deciles with zero significant genes", sum(DE$n_sig == 0), sum(DE$n_sig == 0), "Fig6g_geneformer_deciles.tsv", "sum(n_sig == 0)")
C1g <- rd("m10_abund_check1.tsv", file.path(OUT, "hiend_m10"))[measure == "N_Detections"]
add("C4_gf_ndet_sig", "C4", "Fig6g", "Median detections, significant genes", C1g$median_sig, C1g$median_sig, "N11_out/hiend_m10/m10_abund_check1.tsv", "N_Detections median_sig")
add("C4_gf_ndet_nonsig", "C4", "Fig6g", "Median detections, non-significant genes", C1g$median_nonsig, C1g$median_nonsig, "N11_out/hiend_m10/m10_abund_check1.tsv", "N_Detections median_nonsig")
GM <- rd("Fig6g_geneformer_abundance_matched.tsv")
add("C4_gf_pub", "C4", "Fig6g", "Geneformer positive, published", GM[grp == "Published"]$rate, pct(GM[grp == "Published"]$rate), "Fig6g_geneformer_abundance_matched.tsv", "rate")
add("C4_gf_matched", "C4", "Fig6g", "Geneformer positive, abundance-matched random", GM[grp != "Published"]$rate, pct(GM[grp != "Published"]$rate), "Fig6g_geneformer_abundance_matched.tsv", "rate")
add("C4_gf_matched_P", "C4", "Fig6g", "Fisher P published vs abundance-matched", GM[grp != "Published"]$P, pfmt(GM[grp != "Published"]$P), "Fig6g_geneformer_abundance_matched.tsv", "P")
EQ <- rd("Fig6h_knk_official_equivalence.tsv")
add("C4_knk_fc_rho", "C4", "Fig6h", "scTenifoldKnk shared vs official network, FC Spearman (ALAD)", as.numeric(EQ[grepl("Spearman", check)]$value), EQ[grepl("Spearman", check)]$value, "Fig6h_knk_official_equivalence.tsv", "Spearman row")

## ---- C5 ----
DS <- rd("Fig4d_dsigdb_rates.tsv")
for (g in DS$grp) { r <- DS[grp == g]; add(paste0("C5_dsigdb_", gsub("[^A-Za-z]", "", g)), "C5", "Fig4d", paste("DSigDB positive,", g), r$rate, sprintf("%s (%d/%d)", pct(r$rate), r$pos, r$n), "Fig4d_dsigdb_rates.tsv", "pos/n")
  if (!is.na(r$P)) add(paste0("C5_dsigdb_P_", gsub("[^A-Za-z]", "", g)), "C5", "Fig4d", paste("Fisher P published vs", g), r$P, pfmt(r$P), "Fig4d_dsigdb_rates.tsv", "P") }
AL <- rd("m12_annotation_level.tsv", file.path(OUT, "hiend_m12_annot"))
add("C5_ann_pub", "C5", "Fig4e", "Median drugs annotated per gene, published signatures", AL$median_real, sprintf("%.1f", AL$median_real), "N11_out/hiend_m12_annot/m12_annotation_level.tsv", "median_real")
add("C5_ann_pool", "C5", "Fig4e", "Median drugs annotated per gene, sampling pool", AL$median_U_gene, AL$median_U_gene, "N11_out/hiend_m12_annot/m12_annotation_level.tsv", "median_U_gene")

## ---- C6 ----
SV <- rd("Fig7b_sim_verdict.tsv"); SP <- rd("Fig7b_sim_power.tsv")
add("C6_typeI_framework", "C6", "Fig7b", "Framework false-positive rate at beta = 0", SV$typeI_framework_beta0, sprintf("%s (0/200)", pct(SV$typeI_framework_beta0)), "Fig7b_sim_verdict.tsv", "typeI_framework_beta0")
add("C6_typeI_template", "C6", "Fig7b", "Template univariable false-positive rate at beta = 0", SV$typeI_template_beta0, pct(SV$typeI_template_beta0), "Fig7b_sim_verdict.tsv", "typeI_template_beta0")
r <- SP[criterion == "framework_T1andT2" & beta == max(beta)]
add("C6_power_max", "C6", "Fig7b", sprintf("Framework power at HR/SD %.2f", r$HR_per_SD), r$rate, pct(r$rate), "Fig7b_sim_power.tsv", "framework, max beta")
add("C6_power_max_HR", "C6", "Fig7b", "Largest simulated HR per SD", r$HR_per_SD, f3(r$HR_per_SD), "Fig7b_sim_power.tsv", "max beta")
r1 <- SP[criterion == "T1_random_null" & beta == max(beta)]; r2 <- SP[criterion == "T2_proxy_LRT" & beta == max(beta)]
add("C6_power_gate1", "C6", "Fig7b", "Gate 1 power at largest effect", r1$rate, pct(r1$rate), "Fig7b_sim_power.tsv", "T1, max beta")
add("C6_power_gate2", "C6", "Fig7b", "Gate 2 power at largest effect", r2$rate, pct(r2$rate), "Fig7b_sim_power.tsv", "T2, max beta")
add("C6_pub_adjHR_GSE", "C6", "Fig7b", "Published median adjusted HR, GSE14520", SV$published_adjHR_median_GSE, f3(SV$published_adjHR_median_GSE), "Fig7b_sim_verdict.tsv", "published_adjHR_median_GSE")
add("C6_pub_adjHR_LIRI", "C6", "Fig7b", "Published median adjusted HR, LIRI-JP", SV$published_adjHR_median_LIRI, f3(SV$published_adjHR_median_LIRI), "Fig7b_sim_verdict.tsv", "published_adjHR_median_LIRI")
r02 <- SP[criterion == "framework_T1andT2" & abs(HR_per_SD - 1.2214) < 0.01]; t02 <- SP[criterion == "T0_template_uni" & abs(HR_per_SD - 1.2214) < 0.01]
add("C6_power_HR122_framework", "C6", "Fig7b", "Framework positive rate at HR/SD 1.22", r02$rate, pct(r02$rate), "Fig7b_sim_power.tsv", "framework, beta 0.2")
add("C6_power_HR122_template", "C6", "Fig7b", "Template positive rate at HR/SD 1.22", t02$rate, pct(t02$rate), "Fig7b_sim_power.tsv", "T0, beta 0.2")
E5 <- rd("Fig7c_example_5gene.tsv")[set == "main_TAF9"]
add("C6_example_R2", "C6", "Fig7c", "Worked example, proxy R2 (TCGA)", E5$R2_proxy_TCGA, f3(E5$R2_proxy_TCGA), "Fig7c_example_5gene.tsv", "R2_proxy_TCGA")
for (v in c("empP_GSE", "empP_LIRI", "lrt_GSE", "lrt_LIRI")) add(paste0("C6_example_", v), "C6", "Fig7c", paste("Worked example,", v), E5[[v]], pfmt(E5[[v]]), "Fig7c_example_5gene.tsv", v)

## ---- Sensitivity / supplementary figures ----
CC <- rd("FigS1b_148_vs_95_consistency.tsv")
add("S1_n_same_verdict", "Sens", "FigS1b", "Modules with same verdict at n = 148 and n = 95", sum(CC$agree == "Same verdict"), sprintf("%d/%d", sum(CC$agree == "Same verdict"), nrow(CC)), "FigS1b_148_vs_95_consistency.tsv", "agree")
D14 <- rd("FigS2b_m14_prolif_decision.tsv")
add("S2_raw_pub", "Sens", "FigS2b", "Consensus subtypes positive, published, unadjusted", D14$raw_rate_pub, pct(D14$raw_rate_pub), "FigS2b_m14_prolif_decision.tsv", "raw_rate_pub")
add("S2_raw_Ak", "Sens", "FigS2b", "Consensus subtypes positive, size-matched, unadjusted", D14$raw_rate_Ak, pct(D14$raw_rate_Ak), "FigS2b_m14_prolif_decision.tsv", "raw_rate_Ak")
add("S2_adj_pub", "Sens", "FigS2b", "Consensus subtypes positive, published, proliferation-adjusted", D14$adj_rate_pub, pct(D14$adj_rate_pub), "FigS2b_m14_prolif_decision.tsv", "adj_rate_pub")
add("S2_adj_Ak", "Sens", "FigS2b", "Consensus subtypes positive, size-matched, proliferation-adjusted", D14$adj_rate_Ak, pct(D14$adj_rate_Ak), "FigS2b_m14_prolif_decision.tsv", "adj_rate_Ak")
add("S2_adj_P", "Sens", "FigS2b", "Fisher P after proliferation adjustment", D14$adj_fisher_P, pfmt(D14$adj_fisher_P), "FigS2b_m14_prolif_decision.tsv", "adj_fisher_P")
GF4 <- rd("FigS4a_gene_frequency.tsv")[order(-n_sig)][1:25]
add("S4_top25_prolif", "Sens", "FigS4a", "Hallmark E2F/G2M genes among 25 most reused", sum(GF4$prolif), sprintf("%d/25", sum(GF4$prolif)), "FigS4a_gene_frequency.tsv", "top 25 by n_sig")
add("S4_top_gene_n", "Sens", "FigS4a", paste("Signatures containing", GF4$gene[1]), GF4$n_sig[1], sprintf("%s: %d signatures, %d death types", GF4$gene[1], GF4$n_sig[1], GF4$n_mech[1]), "FigS4a_gene_frequency.tsv", "top 1")

## ---- Methods: corpus flow (parsed row by row from the count table in N12_signatures/README.md) and cohorts ----
CORP_EN <- list(corp_dedup = "Unique PMIDs after deduplication", corp_title = "Records after title screening", corp_cand = "Candidate articles",
                corp_rescue_total = "Candidates after manual review", corp_pmcid = "Articles with a PMCID", corp_abs_only = "PMC abstract only",
                corp_fulltext = "Full text retrieved", corp_coef = "Signatures with per-gene coefficients", corp_primary = "Signatures in the primary analysis",
                corp_nocoef = "Articles without per-gene coefficients", corp_rule_excl = "Excluded by prespecified rules")   # English labels for output; README rows are still matched by their original (Chinese) names

RM <- readLines(file.path(B, "N12_signatures/README.md"), encoding = "UTF-8")
rowv <- function(key) { l <- grep(key, RM, value = TRUE, fixed = TRUE)[1]; stopifnot(!is.na(l)); cells <- trimws(strsplit(l, "|", fixed = TRUE)[[1]])
  x <- gsub("[^0-9]", " ", sub("（.*", "", cells[3])); as.numeric(tail(strsplit(trimws(x), " +")[[1]], 1)) }
for (kv in list(c("corp_dedup", "去重 PMID"), c("corp_title", "标题过滤"), c("corp_cand", "→ 候选"), c("corp_rescue_total", "人工复核补回"),
                c("corp_pmcid", "有 PMCID"), c("corp_abs_only", "PMC 只有摘要"), c("corp_fulltext", "读到全文"), c("corp_coef", "抽到逐基因系数"),
                c("corp_primary", "主表可用"), c("corp_nocoef", "没找到逐基因系数"), c("corp_rule_excl", "规则排除"))) {
  v <- rowv(kv[2]); add(kv[1], "Methods", "Methods corpus", CORP_EN[[kv[1]]], v, v, "N12_signatures/README.md", paste("flow-count row:", CORP_EN[[kv[1]]])) }
stopifnot(rowv("读到全文") == rowv("抽到逐基因系数") + rowv("没找到逐基因系数") + rowv("规则排除"))
add("corp_rescued", "Methods", "Methods corpus", "Titles added after manual review", rowv("人工复核补回") - rowv("→ 候选"), rowv("人工复核补回") - rowv("→ 候选"), "N12_signatures/README.md", "294 - 281")
CM3 <- rd("cohort_meta.tsv", file.path(OUT, "ring3_full"))
add("U_score", "Methods", "Methods genes", "Genes measured in all three cohorts (scoring set)", CM3$n_genes_all[1], format(CM3$n_genes_all[1], big.mark = ","), "N11_out/ring3_full/cohort_meta.tsv", "n_genes_all")
add("U_pool", "Methods", "Methods genes", "Expressed genes in sampling pool (rings 2-4)", CM3$n_genes_universe[1], format(CM3$n_genes_universe[1], big.mark = ","), "N11_out/ring3_full/cohort_meta.tsv", "n_genes_universe")
CA <- rd("covariate_audit.tsv", file.path(OUT, "ring2"))
for (cn in CA$cohort) add(paste0("ncomplete_", cn), "Methods", "Methods cohorts", paste(cn, "patients with complete base covariates"), CA[cohort == cn]$complete_base, CA[cohort == cn]$complete_base, "N11_out/ring2/covariate_audit.tsv", "complete_base")
for (cn in CM$cohort) add(paste0("stageNA_", cn), "Methods", "Methods cohorts", paste(cn, "patients with missing stage"), CM[cohort == cn]$stage_NA, CM[cohort == cn]$stage_NA, "N11_out/robust_v2/cohort_meta.tsv", "stage_NA")
## ---- Methods: signature sizes and numbers of random signatures ----
PC <- rd("published_coverage.tsv", file.path(OUT, "ring3_full"))[n_in >= 2]
stopifnot(nrow(PC) == 148)
add("sig_k_median", "Methods", "Methods signatures", "Published signature genes in scoring set, median", median(PC$n_in), median(PC$n_in), "N11_out/ring3_full/published_coverage.tsv", "median(n_in), n_in>=2")
add("sig_k_range", "Methods", "Methods signatures", "Published signature genes in scoring set, range", max(PC$n_in), sprintf("%d-%d", min(PC$n_in), max(PC$n_in)), "N11_out/ring3_full/published_coverage.tsv", "range(n_in)")
add("sig_full", "Methods", "Methods signatures", "Published signatures with all genes in scoring set", sum(PC$full), sum(PC$full), "N11_out/ring3_full/published_coverage.tsv", "sum(full)")
add("rand_per_sig", "Methods", "Methods random", "Random signatures per published signature per family", M4[module == "M0"]$n_null / M4[module == "M0"]$n_real, M4[module == "M0"]$n_null / M4[module == "M0"]$n_real, "Fig4a_module_rates.tsv", "n_null/n_real")
add("rand_per_sig_M7", "Methods", "Methods random", "Random signatures per published signature, multi-algorithm module", M4[module == "M7"]$n_null / M4[module == "M7"]$n_real, M4[module == "M7"]$n_null / M4[module == "M7"]$n_real, "Fig4a_module_rates.tsv", "n_null/n_real")
H <- { s0 <- strsplit(readLines(file.path(B, "N11_geneset/h.all.v2024.1.Hs.symbols.gmt"), warn = FALSE), "	"); setNames(lapply(s0, function(x) x[-(1:2)]), vapply(s0, `[`, "", 1)) }
Upool <- NULL  # the intersection of proxy sets with the sampling pool needs U; read from hiend_m8/U_pool.txt
Upool <- readLines(file.path(OUT, "hiend_m8/U_pool.txt")); stopifnot(length(Upool) == CM3$n_genes_universe[1])
add("proxy_prolif_n", "Methods", "Methods proxies", "Proliferation proxy genes (E2F + G2M, in pool)", length(intersect(union(H$HALLMARK_E2F_TARGETS, H$HALLMARK_G2M_CHECKPOINT), Upool)),
    length(intersect(union(H$HALLMARK_E2F_TARGETS, H$HALLMARK_G2M_CHECKPOINT), Upool)), "h.all.v2024.1 + hiend_m8/U_pool.txt", "intersect")
add("proxy_oxphos_n", "Methods", "Methods proxies", "OXPHOS proxy genes (in pool)", length(intersect(H$HALLMARK_OXIDATIVE_PHOSPHORYLATION, Upool)),
    length(intersect(H$HALLMARK_OXIDATIVE_PHOSPHORYLATION, Upool)), "h.all.v2024.1 + hiend_m8/U_pool.txt", "intersect")
## ---- Sample sizes for legends ----
add("n_null_default", "N", "Fig4a", "Size-matched random signatures per module (default)", as.numeric(names(which.max(table(M4$n_null)))), format(as.numeric(names(which.max(table(M4$n_null)))), big.mark = ","), "Fig4a_module_rates.tsv", "mode(n_null)")
add("n_null_M7", "N", "Fig4a", "Size-matched random signatures, multi-algorithm module", M4[module == "M7"]$n_null, format(M4[module == "M7"]$n_null, big.mark = ","), "Fig4a_module_rates.tsv", "M7 n_null")
add("n_null_M14", "N", "Fig4a", "Size-matched random signatures, consensus clustering (clusterable)", M4[module == "M14"]$n_null, format(M4[module == "M14"]$n_null, big.mark = ","), "Fig4a_module_rates.tsv", "M14 n_null")
P3 <- rd("Fig3ab_per_signature.tsv")[cohort == "TCGA", .N, by = grp]
for (g in P3$grp) add(paste0("n_fig3_", gsub("[^A-Za-z]", "", g)), "N", "Fig3", paste("Signatures per cohort,", g), P3[grp == g]$N, format(P3[grp == g]$N, big.mark = ","), "Fig3ab_per_signature.tsv", "TCGA rows by grp")
E6 <- rd("Fig6e_sc_relative_localisation.tsv")[, .(n = uniqueN(sid)), by = grp]
for (g in E6$grp) add(paste0("n_fig6e_", gsub("[^A-Za-z]", "", g)), "N", "Fig6e", paste("Signatures scored in single cells,", g), E6[grp == g]$n, E6[grp == g]$n, "Fig6e_sc_relative_localisation.tsv", "uniqueN(sid)")
add("n_gf_matched", "N", "Fig6g", "Abundance-matched random signatures", GM[grp != "Published"]$n, format(GM[grp != "Published"]$n, big.mark = ","), "Fig6g_geneformer_abundance_matched.tsv", "n")
add("n_sim_runs", "N", "Fig7b", "Simulation runs per effect size", SP$n[1], SP$n[1], "Fig7b_sim_power.tsv", "n")
T3 <- rbindlist(N)
stopifnot(!anyDuplicated(T3$id))
T3[, source := gsub("N28_fig/data/", "figures/data/", gsub("N11_scripts/", "analysis/", source, fixed = TRUE), fixed = TRUE)]
fwrite(T3, file.path(FDIR, "TableS3_numbers.tsv"), sep = "\t")
back <- fread(file.path(FDIR, "TableS3_numbers.tsv"), sep = "\t")
stopifnot(nrow(back) == nrow(T3), isTRUE(all.equal(back$value, T3$value)))
cat("TABLE S3 OK:", nrow(T3), "numbers\n")
