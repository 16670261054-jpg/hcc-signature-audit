# Fig 4  template downstream modules call most random signatures "positive" + full template display of one random signature + DSigDB counterexample
# Reads: ring3_armA_decision_v3_rescaled.tsv (M0-M7, analysis = fitsig_pool100_k|all), hiend_m3_m12/, hiend_m8/,
#         hiend_m13_m15/, hiend_m10/ (rows A_vs_Ak_sizematched), hiend_m12_annot/, figures/data/fig4c_*.tsv
# Statistics: real vs random Fisher's exact (read from disk); random positive rate Clopper-Pearson 95% CI; KM log-rank (read from disk).
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggh4x); library(ggrepel) })

## ---- a module summary ----
r3 <- fread(file.path(OUT, "ring3_armA_decision_v3_rescaled.tsv"))[analysis == "fitsig_pool100_k|all"]
r3 <- r3[, .(module = sub("_pos", "", module), rate_real, n_real, rate_null = FPR_null, n_null, lo = FPR_CI_lo, hi = FPR_CI_hi, P = fisher_P)]
pick <- function(f, mods, cols = c("FPR_CP_lo", "FPR_CP_hi")) {
  d <- fread(f)[comparison == "A_vs_Ak_sizematched" & module %in% mods]
  d[, .(module = sub("_pos$", "", module), rate_real, n_real, rate_null = FPR_null, n_null, lo = get(cols[1]), hi = get(cols[2]), P = fisher_P)]
}
m3 <- pick(file.path(OUT, "hiend_m3_m12/m3_m12_decision.tsv"), c("M3_pos", "M12_pos"))
m8 <- pick(file.path(OUT, "hiend_m8/m8_decision.tsv"), c("M8a_pos", "M8b_pos"))
m13 <- pick(file.path(OUT, "hiend_m13_m15/m13_m15_decision.tsv"), c("M13_pos", "M14_pos", "DCA_pos"))
m13[module == "DCA", module := "M15"]
m10 <- fread(file.path(OUT, "hiend_m10/m10_decision.tsv"))[comparison == "A_vs_Ak_sizematched" & filter == "all",
        .(module = "M10", rate_real, n_real, rate_null = FPR_null, n_null, lo = FPR_CP_lo, hi = FPR_CP_hi, P = fisher_P)]
M <- rbind(r3, m3, m8, m13, m10)
lab <- fread(text = "module|label|category
M0|External validation (C > 0.5, both cohorts)|Prognostic model
M1|Independent prognostic factor (multivariable Cox)|Prognostic model
M7|Multi-algorithm model selection|Prognostic model
M14|Consensus-clustering subtypes differ in survival|Prognostic model
M6|Nomogram adds to stage (nested LRT)|Clinical utility
M15|Decision curve shows net benefit|Clinical utility
M2|Immune infiltration differs by risk group|Immune & pathway
M13|Immune checkpoint expression differs|Immune & pathway
M4|GSEA Hallmark pathways enriched|Immune & pathway
M5a|TMB differs by risk group|Genomic
M5b|TMB x risk groups separate survival|Genomic
M3|Drug sensitivity differs (oncoPredict)|Drug
M12|Drugs enriched (DSigDB)|Drug
M8a|Signature genes mark a cell type (scRNA)|Single cell
M8b|Signature genes differ tumour vs normal (scRNA)|Single cell
M10|Geneformer in silico deletion|Single cell", sep = "|")
M <- merge(M, lab, by = "module"); stopifnot(nrow(M) == 16)
## Between-group test: CMH stratified by matched group + Holm correction across 16 modules; positive rates and CIs unchanged
CMH <- fread(file.path(OUT, "cmh_modules.tsv"))
stopifnot(all(M$module %in% CMH$module), max(abs(merge(M, CMH, by = "module")[, rate_real - rate_pub])) < 1e-9)
M <- merge(M[, !"P"], CMH[, .(module, P = CMH_P_holm, P_unadj = CMH_P, CMH_OR)], by = "module")
M[, star := stars(P)]
write_sd(M[order(-rate_null)], "Fig4a_module_rates")
COL_CAT <- setNames(SET2[1:6], c("Prognostic model", "Clinical utility", "Immune & pathway", "Genomic", "Drug", "Single cell"))
M[, label := factor(label, levels = M[order(rate_null)]$label)]
M[, category := factor(category, levels = names(COL_CAT))]
pa <- ggplot(M, aes(y = label)) +
  annotate("rect", xmin = 0, xmax = BRADLEY, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_tile(aes(x = -0.07, fill = category), width = 0.035, height = 0.8) +
  geom_segment(aes(x = rate_real, xend = rate_null, yend = label), colour = "grey55", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = COL_GRP[["Size-matched random"]], linewidth = 1.6, alpha = 0.35) +
  geom_point(aes(x = rate_null), shape = 21, size = 2.3, fill = COL_GRP[["Size-matched random"]], stroke = 0.25) +
  geom_point(aes(x = rate_real), shape = 21, size = 1.9, fill = COL_GRP[["Published"]], stroke = 0.25) +
  geom_text(aes(x = 1.1, label = star), size = 2, hjust = 0.5) +
  annotate("text", x = 1.1, y = 16.75, label = "vs random", size = 1.8, vjust = 0) +
  scale_fill_manual("Module type", values = COL_CAT) +
  scale_x_continuous("Signatures called positive (template criterion)", labels = scales::percent,
                     breaks = seq(0, 1, 0.25), limits = c(-0.09, 1.17), expand = c(0, 0)) +
  scale_y_discrete(NULL) + coord_cartesian(clip = "off") +
  labs(title = NULL) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 5.5),
        axis.text.y = element_text(size = 6), legend.key.size = unit(2.4, "mm"),
        plot.margin = margin(12, 4, 2, 2)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

## ---- c template display of one random signature ----
KM <- fread(file.path(DAT, "fig4c_km_curves.tsv")); LR <- fread(file.path(DAT, "fig4c_km_logrank.tsv"))
AUC <- fread(file.path(DAT, "fig4c_timeroc_auc.tsv")); ROC <- fread(file.path(DAT, "fig4c_timeroc_curves.tsv"))
MC <- fread(file.path(DAT, "fig4c_multicox.tsv"))
write_sd(LR, "Fig5a_km_logrank"); write_sd(AUC, "Fig5b_timeroc_auc"); write_sd(MC, "Fig5c_multicox")
km_plot <- function(cn) {
  d <- KM[cohort == cn]; l <- LR[cohort == cn]
  ggplot(d, aes(time / 12, surv, colour = strata)) +
    geom_step(linewidth = 0.45) +
    geom_point(data = d[n_censor > 0], shape = 3, size = 0.5, stroke = 0.25) +
    annotate("text", x = Inf, y = 0.98, hjust = 1, vjust = 1, size = 1.9, family = "Arial",
             label = sprintf("log-rank %s\nHR %.2f", stars(l$logrank_P), l$HR_high_vs_low)) +
    scale_colour_manual(NULL, values = COL_RISK) +
    scale_y_continuous("Overall survival", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_continuous("Years", breaks = seq(0, 10, 2)) +
    labs(title = sprintf("%s (n = %d)", c(TCGA = "TCGA-LIHC", GSE14520 = "GSE14520", LIRI = "LIRI-JP")[[cn]], l$n)) +
    theme(legend.position = "none", plot.title = element_text(size = 6.2))
}
km_strip_y <- function(p) p + theme(axis.title.y = element_blank())
pk <- lapply(c("TCGA", "GSE14520", "LIRI"), km_plot); pk[[2]] <- km_strip_y(pk[[2]]); pk[[3]] <- km_strip_y(pk[[3]])
COL_YR <- COL_YEAR
ROC[, yr := factor(year)]; AUC[, yr := factor(year)]
pr <- ggplot(ROC[order(yr, FPR, TPR)], aes(FPR, TPR, colour = yr)) +
  geom_abline(linetype = REF_LTY, colour = COL_REF, linewidth = 0.3) +
  geom_path(linewidth = 0.45) +
  geom_text(data = AUC, aes(x = 0.98, y = 0.05 + 0.11 * (3 - as.numeric(yr)), label = sprintf("%d-yr AUC %.2f", year, AUC)),
            hjust = 1, size = 1.8, family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = COL_YR, guide = "none") +
  scale_x_continuous("1 - specificity", breaks = c(0, 0.5, 1)) + scale_y_continuous("Sensitivity", breaks = c(0, 0.5, 1)) +
  labs(title = "TCGA-LIHC") + theme(plot.title = element_text(size = 6.2))
MC[, term_l := factor(c(Risk = "Risk score (per SD)", Age = "Age (per 10 y)", SexMALE = "Male sex", "StageIII-IV" = "Stage III-IV")[term],
                     levels = rev(c("Risk score (per SD)", "Age (per 10 y)", "Male sex", "Stage III-IV")))]
pf <- ggplot(MC, aes(HR, term_l)) +
  geom_ref_v(1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, colour = "grey30") +
  geom_point(aes(fill = term == "Risk"), shape = 22, size = 1.9, stroke = 0.25) +
  geom_text(aes(x = 3.9, label = stars(P)), size = 1.9, hjust = 1) +
  scale_fill_manual(values = c(`TRUE` = COL_GRP[["Size-matched random"]], `FALSE` = "white"), guide = "none") +
  scale_x_log10("Hazard ratio (95% CI)", breaks = c(0.5, 1, 2, 4), limits = c(0.5, 4)) + scale_y_discrete(NULL) +
  labs(title = "TCGA-LIHC") + theme(plot.title = element_text(size = 6.2), axis.text.y = element_text(size = 5.5))
g4 <- fread(file.path(OUT, "hiend_m8/signature_genes_all_arms.tsv"))[arm == "Ak" & src == "S004" & rep == 1]$gene
klg <- ggplot(data.table(x = c(1, 2.6), s = names(COL_RISK)), aes(x, 1, colour = s)) +
  geom_segment(aes(xend = x + 0.35, yend = 1), linewidth = 0.6) +
  geom_text(aes(x = x + 0.45, label = s), hjust = 0, size = 1.9, colour = "black", family = "Arial") +
  annotate("text", x = 4.1, y = 1, hjust = 0, size = 1.9, family = "Arial", fontface = "italic",
           label = paste(g4, collapse = " + ")) +
  scale_colour_manual(values = COL_RISK, guide = "none") + xlim(0.9, 7.5) + theme_void()
pc <- wrap_plots(list(klg, pk[[1]], pk[[2]], pk[[3]], pr, pf),
                 design = "AAAAA\nBCDEF", heights = c(0.1, 1), widths = c(1, 1, 1, 0.95, 1.3))
## ---- d DSigDB: the only counterexample ----
AD <- fread(file.path(OUT, "hiend_m12_annot/m12_annot_decision.tsv"))[group == "A_published"]
dd <- data.table(grp = c("Published", "Random", "Annotation-\nmatched"),
                 pos = round(c(AD$rate_real[1] * AD$n_real[1], AD[null == "N_plain"]$rate_null * AD[null == "N_plain"]$n_null,
                               AD[null == "N_annot"]$rate_null * AD[null == "N_annot"]$n_null)),
                 n = c(AD$n_real[1], AD[null == "N_plain"]$n_null, AD[null == "N_annot"]$n_null),
                 P = c(NA, AD[null == "N_plain"]$fisher_P, AD[null == "N_annot"]$fisher_P))
## Between-group test: CMH stratified by matched set (P read from N35_paired_tests); positive rates and CIs unchanged
PT <- fread(file.path(B, "N35_paired_tests/paired_tests_results.tsv"), sep = "	", quote = "")
dd[, key := c(NA, "Fig4d|N_plain", "Fig4d|N_annot")][!is.na(key), `:=`(P = PT$new_P[match(key, PT$key)], test = PT$new_test[match(key, PT$key)])]
stopifnot(!anyNA(dd$P[-1])); dd[, key := NULL]
dd[, rate := pos / n][, c("lo", "hi") := as.list(binom.test(pos, n)$conf.int), by = grp][, star := stars(P)]
dd[, grp := factor(grp, levels = grp)]
write_sd(dd, "Fig4d_dsigdb_rates")
AR <- fread(file.path(OUT, "hiend_m12_annot/m12_annot_rows.tsv"))
AR[, grp := factor(fcase(type == "real", "Published", type == "N_plain", "Random", type == "N_annot", "Annotation-\nmatched"), levels = levels(dd$grp))]
AR <- AR[arm == "A" & !is.na(grp)]
write_sd(AR[, .(n = .N, median_mean_ann = median(mean_ann)), by = grp], "Fig4e_dsigdb_annotation_degree")
COL_D <- c(COL_GRP[["Published"]], COL_GRP[["Input-matched random"]], COL_GRP[["Size-matched random"]]); names(COL_D) <- levels(dd$grp)
pd1 <- ggplot(dd, aes(grp, rate)) +
  geom_col(aes(fill = grp), width = 0.62, colour = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.3) +
  geom_text(aes(y = hi + 0.035, label = star), size = 2) +
  scale_fill_manual(values = COL_D, guide = "none") +
  scale_y_continuous("Signatures with ≥1 enriched drug", labels = scales::percent, limits = c(0, 0.8), expand = c(0, 0)) +
  scale_x_discrete(NULL) + labs(title = "DSigDB enrichment") + theme(axis.text.x = element_text(size = 5.5))
set.seed(SEED)
pd2 <- ggplot(AR, aes(grp, mean_ann)) +
  geom_violin(aes(fill = grp), alpha = 0.75, colour = NA, width = 0.8, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = COL_D, guide = "none") +
  scale_y_log10("Drugs annotated per gene\n(signature mean, log scale)") + scale_x_discrete(NULL) +
  labs(title = NULL) + theme(axis.text.x = element_text(size = 5.5))

## ---- Legend + assembly ----
lg <- ggplot(data.table(x = c(1, 4.6), g = c("Published", "Size-matched random")), aes(x, 1, fill = g)) +
  geom_point(shape = 21, size = 2.2, stroke = 0.25) +
  geom_text(aes(x = x + 0.15, label = c("Published signatures (n = 148)", "Size-matched random (95% CI)")), hjust = 0, size = 2, family = "Arial") +
  annotate("rect", xmin = 8.2, xmax = 8.6, ymin = 0.8, ymax = 1.2, fill = "grey90") +
  annotate("text", x = 8.7, y = 1, label = "Acceptable false-positive range (\u2264 7.5%)", hjust = 0, size = 2, family = "Arial") +
  scale_fill_manual(values = COL_GRP, guide = "none") + xlim(0.8, 12.5) + ylim(0.6, 1.4) + theme_void()
right <- wrap_plots(list(pd1 + labs(tag = "c"), pd2), ncol = 1)
top <- wrap_plots(list(pa + labs(tag = "a"), right), widths = c(1.55, 1))
ttl <- ggplot() + annotate("text", x = 0, y = 0, hjust = 0, size = 2.6, fontface = "bold", family = "Arial",
                           label = "One size-matched random signature, run through the template") + xlim(0, 1) + theme_void()
p <- wrap_plots(list(lg, top, wrap_elements(full = ttl), wrap_elements(full = pc) + labs(tag = "b")), ncol = 1, heights = c(0.05, 1, 0.04, 0.42))
save_panel(pa, "Fig4a_module_dumbbell", 115, 105); for (i in 1:3) save_panel(pk[[i]], paste0("Fig4b", i, "_KM_", c("TCGA", "GSE14520", "LIRI")[i]), 55, 55)
save_panel(pr, "Fig4b4_timeROC", 55, 55); save_panel(pf, "Fig4b5_multicox", 70, 55); save_panel(pd1, "Fig4c1_dsigdb_rates", 65, 62); save_panel(pd2, "Fig4c2_dsigdb_annotation", 65, 62)
saveRDS(list(pa = pa, pk = pk, pr = pr, pf = pf, pd1 = pd1, pd2 = pd2, lg = lg, klg = klg), file.path(DAT, "obj_fig4.rds"))  # assembled in assemble_main.R
cat("FIG4 OK\n")
