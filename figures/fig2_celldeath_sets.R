# Fig 2  cell-death gene sets are not better than random gene sets of the same size
# Reads: N11_out/robust_v2/{inner_splits,external,proxy_external,adjusted_null,cohort_meta}.tsv
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggdist); library(ggh4x) })
R <- file.path(OUT, "robust_v2")   # common gene pool (10,651 genes); the earlier version is in N11_out/robust/
SET_LAB <- c(cuproptosis_Tsvetkov19 = "Cuproptosis", disulfidptosis_10 = "Disulfidptosis",
             WP_FERROPTOSIS = "Ferroptosis (WP)", GOBP_FERROPTOSIS = "Ferroptosis (GO)",
             REACTOME_PYROPTOSIS = "Pyroptosis", GOBP_NECROPTOTIC_PROCESS = "Necroptosis",
             REACTOME_REGULATED_NECROSIS = "Regulated necrosis", GOBP_ANOIKIS = "Anoikis")
lab_set <- function(x) { m <- SET_LAB[x]; m[is.na(m)] <- x[is.na(m)]; unname(m) }

## a 50 internal splits: distribution of empirical P (uniform under the null)
I <- fread(file.path(R, "inner_splits.tsv"))
I[, set_l := lab_set(set)]
stopifnot(I[, .N, by = set]$N == 50)
ord <- I[, .(med = median(emp_p)), by = set_l][order(-med)]$set_l
I[, set_l := factor(set_l, levels = ord)]
sa <- I[, .(n_split = .N, k = k[1], median_emp_p = median(emp_p), frac_emp_p_lt05 = mean(emp_p < 0.05),
            median_c_obs = median(c_obs), median_null_med = median(null_med)), by = set_l]
write_sd(sa, "Fig2a_inner_splits_summary"); write_sd(I[, .(set_l, split, k, n_model, c_obs, null_med, null_p95, emp_p)], "Fig2a_inner_splits")
set.seed(SEED); I[, yj := as.numeric(set_l) - 0.17 + runif(.N, -0.06, 0.06)]
pa <- ggplot(I, aes(emp_p, as.numeric(set_l))) +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#FDDBC7", alpha = 0.6) +
  stat_slab(aes(y = set_l), fill = COL_GRP[["Cell-death gene set"]], alpha = 0.55, scale = 0.62,
            density = "bounded", normalize = "groups", position = position_nudge(y = 0.08), colour = NA) +
  geom_point(aes(y = yj), shape = 21, size = 0.7, stroke = 0.15, fill = COL_GRP[["Cell-death gene set"]],
             colour = "white", alpha = 0.9) +
  stat_pointinterval(aes(y = set_l), .width = c(0.5), point_size = 1.2, interval_size = 0.9,
                     position = position_nudge(y = 0.03), colour = "black") +
  geom_ref_v(0.05) +
  scale_x_continuous("Empirical P vs equal-size random gene sets\n(50 TCGA train/test splits)",
                     limits = c(0, 1), breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), labels = c("0", ".05", ".25", ".5", ".75", "1"),
                     expand = c(0.01, 0)) +
  scale_y_discrete(NULL) + labs(title = "TCGA-LIHC, 50 splits") +
  theme(axis.text.y = element_text(size = 6.5))

## b external cohorts: observed C and random null (median -> p95); proliferation score (untrained) as reference
E <- fread(file.path(R, "external.tsv")); P <- fread(file.path(R, "proxy_external.tsv"))
E[, set_l := factor(lab_set(set), levels = ord)]
E[, cohort_l := factor(COH_LAB[cohort], levels = COH_LAB[c("GSE14520", "LIRI")])]
P <- P[proxy == "prolif"][, cohort_l := factor(COH_LAB[cohort], levels = levels(E$cohort_l))]
write_sd(E[, .(set_l, cohort, k, n_model, c_ext, null_med, null_p95, emp_p)], "Fig2b_external")
write_sd(P, "Fig2b_proliferation_reference")
pb <- ggplot(E, aes(y = set_l)) +
  geom_segment(aes(x = null_med, xend = null_p95, yend = set_l), linewidth = 2.6, colour = COL_NULL, lineend = "round") +
  geom_point(aes(x = null_p95), shape = 124, size = 2.6, colour = "grey35") +
  geom_vline(data = P, aes(xintercept = c), colour = COL_PROXY, linewidth = 0.35, linetype = "22") +
  geom_point(aes(x = c_ext), shape = 21, size = 2, fill = COL_GRP[["Cell-death gene set"]], colour = "black", stroke = 0.25) +
  facet_wrap2(~ cohort_l, nrow = 1, strip = strip_cohort(c("GSE14520", "LIRI"))) +
  scale_x_continuous("Harrell's C in external cohort", breaks = seq(0.5, 0.75, 0.05)) +
  scale_y_discrete(NULL) + labs(title = NULL) +
  theme(axis.text.y = element_text(size = 6.5), panel.spacing.x = unit(3, "mm"))

## c after adjustment for proliferation + OXPHOS + stage: |z| vs p95 of random sets adjusted the same way
A <- fread(file.path(R, "adjusted_null.tsv"))
A[, set_l := factor(lab_set(set), levels = ord)][, cohort_l := factor(COH_LAB[cohort], levels = levels(E$cohort_l))]
write_sd(A[, .(set_l, cohort, k, z_adj, P_adj, null_absz_p95, null_frac_Padj_lt05, emp_p_absz)], "Fig2c_adjusted_null")
pc <- ggplot(A, aes(y = set_l)) +
  geom_segment(aes(x = 0, xend = null_absz_p95, yend = set_l), linewidth = 2.6, colour = COL_NULL, lineend = "round") +
  geom_point(aes(x = null_absz_p95), shape = 124, size = 2.6, colour = "grey35") +
  geom_ref_v(qnorm(0.975)) +
  geom_point(aes(x = abs(z_adj)), shape = 21, size = 2, fill = COL_GRP[["Cell-death gene set"]], colour = "black", stroke = 0.25) +
  facet_wrap2(~ cohort_l, nrow = 1, strip = strip_cohort(c("GSE14520", "LIRI"))) +
  scale_x_continuous("|z| after adjustment for proliferation, OXPHOS and stage", breaks = 0:3) +
  scale_y_discrete(NULL) + labs(title = NULL) +
  theme(axis.text.y = element_text(size = 6.5), panel.spacing.x = unit(3, "mm"))

## Legend (one hand-made legend instead of one per panel)
lg <- data.table(x = c(1, 2.2, 3.6, 5.2), lab = c("Cell-death gene set", "Random null (to 95th percentile)",
                                              "Proliferation score alone", "Nominal P = 0.05 (a, c)"))
pl <- ggplot() +
  geom_point(aes(1, 1), shape = 21, size = 2, fill = COL_GRP[["Cell-death gene set"]], stroke = 0.25) +
  geom_segment(aes(x = 2.75, xend = 3.15, y = 1, yend = 1), linewidth = 2.6, colour = COL_NULL, lineend = "round") +
  geom_segment(aes(x = 5.2, xend = 5.5, y = 1, yend = 1), colour = COL_PROXY, linewidth = 0.35, linetype = "22") +
  geom_segment(aes(x = 7.0, xend = 7.3, y = 1, yend = 1), colour = COL_REF, linewidth = 0.3, linetype = REF_LTY) +
  geom_text(aes(x = c(1.12, 3.25, 5.6, 7.4), y = 1, label = lg$lab), hjust = 0, size = 2.1, family = "Arial") +
  xlim(0.9, 8.6) + theme_void()

pa <- pa + labs(tag = "a"); pb <- pb + labs(tag = "b"); pc <- pc + labs(tag = "c")
p <- wrap_plots(A = pa, B = pb, C = pc, L = pl, design = "AB
AC
LL", heights = c(1, 1, 0.07), widths = c(1, 1.25))
save_panel(pa, "Fig2a_internal", 85, 100); save_panel(pb, "Fig2b_external", 100, 60); save_panel(pc, "Fig2c_adjusted", 100, 60)
save_pub(p, "Fig2_celldeath_sets", width_mm = 180, height_mm = 110)
cat("FIG2 OK\n")
