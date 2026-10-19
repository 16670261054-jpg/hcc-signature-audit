# Fig 3  template selection pulls signatures towards the proliferation axis; residual signal after proxy adjustment is not more frequent than for random signatures
# Reads: N11_out/ring2/{ring2_rows,ring2_R2_summary,ring2_residual_summary}.tsv
# Statistics: continuous outcomes = paired Wilcoxon signed-rank, proportions = CMH stratified by matched set; P read from N35_paired_tests/paired_tests_results.tsv (analysis/paired_tests.R);
#       random positive rate with Clopper-Pearson 95% CI.
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggdist); library(ggh4x) })
R2 <- file.path(OUT, "ring2")
X <- fread(file.path(R2, "ring2_rows.tsv"))
X[, grp := ARM2GRP[paste(arm, kind, sep = "_")]]
X <- X[arm %in% c("A", "Ak")]   # main analysis: 148 published vs the random signatures
G3 <- c("Published", "Size-matched random", "Input-matched random")
X[, grp := factor(grp, levels = G3)][, cohort_l := factor(COH_LAB[cohort], levels = COH_LAB)]
stopifnot(X[, .N, by = .(grp, cohort)][grp == "Published"]$N == 148)

## P values of two-group comparisons: R² read from disk; rho_prolif computed here (same MWU, two-sided)
S <- fread(file.path(R2, "ring2_R2_summary.tsv"))
S <- S[comparison %in% c("A_vs_Ak_sizematched", "A_vs_fitsig_k")]
S[, grp2 := fifelse(comparison == "A_vs_Ak_sizematched", "Size-matched random", "Input-matched random")]
rho_p <- CJ(cohort = unique(X$cohort), g2 = G3[2:3])[, .(P = wilcox.test(X[cohort == .BY$cohort & grp == "Published"]$rho_prolif,
                                                                         X[cohort == .BY$cohort & grp == .BY$g2]$rho_prolif)$p.value), by = .(cohort, g2)]
stat <- rbind(S[, .(metric = "R2_proxy", cohort, grp2, P = MWU_P, median_pub = R2_median_real, median_null = R2_median_null)],
              merge(rho_p[, .(cohort, grp2 = as.character(g2), P)],
                    X[, .(median = median(rho_prolif)), by = .(cohort, grp)][, .(cohort, grp = as.character(grp), median)][
                      , .(median_pub = median[grp == "Published"], grp2 = grp, median_null = median), by = cohort][grp2 != "Published"],
                    by = c("cohort", "grp2"))[, .(metric = "rho_prolif", cohort, grp2, P, median_pub, median_null)])
## Between-group test: paired Wilcoxon signed-rank (each published signature vs the median of its matched random signatures), P read from N35_paired_tests; medians unchanged
PT <- fread(file.path(B, "N35_paired_tests/paired_tests_results.tsv"), sep = "	", quote = "")
stat[, P := PT$new_P[match(paste("Fig3ab", metric, cohort, grp2, sep = "|"), PT$key)]]; stopifnot(!anyNA(stat$P))
stat[, `:=`(test = PT$new_test[match(paste("Fig3ab", metric, cohort, grp2, sep = "|"), PT$key)], star = stars(P))]
write_sd(stat, "Fig3ab_stats")
write_sd(X[, .(grp, src, rep, cohort, k_model, R2_proxy, rho_prolif)], "Fig3ab_per_signature")

rain <- function(var, ylab, ttl, ybr) {
  d <- X[, .(grp, cohort_l, cohort, v = get(var))]
  set.seed(SEED); d[, xj := as.numeric(grp) - 0.2 + runif(.N, -0.07, 0.07)]
  st <- stat[metric == var][, cohort_l := factor(COH_LAB[cohort], levels = COH_LAB)]
  st[, `:=`(x0 = 1, x1 = match(grp2, G3))]
  ymax <- max(d$v, na.rm = TRUE); st[, y := ymax + fifelse(x1 == 2, 0.06, 0.14) * diff(range(d$v, na.rm = TRUE))]
  ggplot(d, aes(grp, v)) +
    geom_point(data = d[grp != "Published"], aes(x = xj, colour = grp), size = 0.12, alpha = 0.18, stroke = 0) +
    geom_point(data = d[grp == "Published"], aes(x = xj), shape = 21, fill = COL_GRP[["Published"]], colour = "white", size = 0.75, stroke = 0.1, alpha = 0.9) +
    stat_slab(aes(fill = grp), side = "right", scale = 0.55, alpha = 0.85, colour = NA, normalize = "groups",
              position = position_nudge(x = 0.02)) +
    geom_boxplot(aes(fill = grp), width = 0.13, outlier.shape = NA, linewidth = 0.28, alpha = 0.95,
                 position = position_nudge(x = -0.05), colour = "black") +
    geom_segment(data = st, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.3, inherit.aes = FALSE) +
    geom_text(data = st, aes(x = (x0 + x1) / 2, y = y, label = star), vjust = -0.15, size = 2.2, inherit.aes = FALSE) +
    facet_wrap2(~ cohort_l, nrow = 1, strip = strip_cohort(names(COL_COHORT))) +
    scale_fill_manual(values = COL_GRP, guide = "none") + scale_colour_manual(values = COL_GRP, guide = "none") +
    scale_x_discrete(NULL, labels = c("Published", "Size-matched\nrandom", "Input-matched\nrandom")) +
    scale_y_continuous(ylab, breaks = ybr, expand = expansion(mult = c(0.02, 0.1))) + labs(title = NULL) +
    theme(axis.text.x = element_text(size = 6), panel.spacing.x = unit(3, "mm"))
}
pa <- rain("R2_proxy", expression("Variance explained by proxies ("*italic(R)^2*")"),
           "Proxy loading of signature scores (proliferation, stroma, immune, OXPHOS)", seq(0, 1, 0.25))
pb <- rain("rho_prolif", expression("Spearman "*rho*" with proliferation"),
           "Template selection pulls random signatures towards the proliferation axis", seq(-0.5, 1, 0.25))

## c proportion with LRT < 0.05 in both external cohorts after proxy adjustment
Rs <- fread(file.path(R2, "ring2_residual_summary.tsv"))
Rs <- Rs[comparison %in% c("A_vs_Ak_sizematched", "A_vs_fitsig_k")]
ML <- c("Base: proliferation, stroma,\nimmune, stage, age, sex", "Base + OXPHOS")
Rs[, model := factor(fifelse(base_model == "+oxphos", ML[2], ML[1]), levels = ML)]
cp <- function(x, n) binom.test(x, n)$conf.int
cd <- rbind(Rs[comparison == "A_vs_Ak_sizematched", .(grp = "Published", model, pos = pos_real, n = n_real, P = NA_real_)],
            Rs[, .(grp = fifelse(comparison == "A_vs_Ak_sizematched", "Size-matched random", "Input-matched random"),
                   model, pos = pos_null, n = n_null, P = fisher_P)])
cd[, rate := pos / n][, c("lo", "hi") := as.list(cp(pos, n)), by = .(grp, model)]
## Between-group test: CMH stratified by matched set (P read from N35_paired_tests); positive rates and CIs unchanged
cd[, key := paste("Fig3c", fifelse(model == ML[2], "lrt_p_oxphos", "lrt_p"), grp, sep = "|")]
cd[grp != "Published", `:=`(P = PT$new_P[match(key, PT$key)], test = PT$new_test[match(key, PT$key)])]; stopifnot(!anyNA(cd[grp != "Published"]$P))
cd[, star := stars(P)][, grp := factor(grp, levels = G3)]
write_sd(cd[, .(grp, model, pos, n, rate, CP_lo = lo, CP_hi = hi, P_vs_published = P, test, star)], "Fig3c_residual")
pc <- ggplot(cd, aes(grp, rate, colour = grp)) +
  geom_ref_h(0.05) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6) +
  geom_point(aes(fill = grp), shape = 21, size = 2.2, colour = "black", stroke = 0.25) +
  geom_text(aes(y = hi, label = star), vjust = -0.6, size = 2.1, colour = "black") +
  facet_wrap(~ model, nrow = 1) +
  scale_fill_manual(values = COL_GRP, guide = "none") + scale_colour_manual(values = COL_GRP, guide = "none") +
  scale_x_discrete(NULL, labels = c("Published", "Size-\nmatched", "Input-\nmatched")) +
  scale_y_continuous("Signatures with residual signal\nin both external cohorts", labels = scales::percent,
                     limits = c(0, 0.11), breaks = seq(0, 0.1, 0.025)) +
  labs(title = NULL) +
  theme(axis.text.x = element_text(size = 6))

## d HR attenuation (descriptive): |log HR| of published signatures before and after adjustment, paired
H <- X[grp == "Published" & cohort != "TCGA", .(src, cohort_l, cohort, u = abs(log(hr_unadj)), a = abs(log(hr_adj)))]
Hn <- X[grp == "Size-matched random" & cohort != "TCGA", .(u = median(abs(log(hr_unadj)), na.rm = TRUE),
                                                          a = median(abs(log(hr_adj)), na.rm = TRUE)), by = cohort_l]
Hm <- H[, .(u = median(u, na.rm = TRUE), a = median(a, na.rm = TRUE)), by = cohort_l]
write_sd(H, "Fig3d_HR_pairs"); write_sd(rbind(Hm[, grp := "Published"], Hn[, grp := "Size-matched random"]), "Fig3d_HR_medians")
hl <- function(d) melt(d, measure.vars = c("u", "a"), variable.name = "state")[, x := fifelse(state == "u", 1, 2)]
Hl <- hl(H); Hml <- hl(Hm); Hnl <- hl(Hn)
pd <- ggplot(Hl, aes(x, value)) +
  geom_line(aes(group = src), colour = "#F4A582", linewidth = 0.18, alpha = 0.55) +
  geom_line(data = Hnl, aes(group = cohort_l), colour = COL_GRP[["Size-matched random"]], linewidth = 0.8, linetype = "22") +
  geom_line(data = Hml, aes(group = cohort_l), colour = COL_GRP[["Published"]], linewidth = 0.9) +
  geom_point(data = Hml, shape = 21, fill = COL_GRP[["Published"]], size = 2, stroke = 0.25) +
  facet_wrap2(~ cohort_l, nrow = 1, strip = strip_cohort(c("GSE14520", "LIRI"))) +
  scale_x_continuous(NULL, breaks = 1:2, labels = c("Unadjusted", "Proxy-\nadjusted"), expand = expansion(add = 0.25)) +
  scale_y_continuous("|log HR| per SD", limits = c(0, NA)) +
  labs(title = NULL) +
  theme(axis.text.x = element_text(size = 6), panel.spacing.x = unit(3, "mm"))

## Group legend
lg <- data.table(grp = factor(G3, levels = G3), x = c(1, 4.2, 7.6), lab = c("Published (n = 148)", "Size-matched random", "Input-matched random"))
pl <- ggplot(lg, aes(x, 1, fill = grp)) + geom_point(shape = 21, size = 2.3, stroke = 0.25) +
  geom_text(aes(x = x + 0.2, label = lab), hjust = 0, size = 2.2, family = "Arial") +
  annotate("segment", x = 11.1, xend = 11.7, y = 1, yend = 1, colour = COL_GRP[["Size-matched random"]], linetype = "22", linewidth = 0.6) +
  annotate("text", x = 11.85, y = 1, label = "Size-matched random, median (d)", hjust = 0, size = 2.2, family = "Arial") +
  scale_fill_manual(values = COL_GRP, guide = "none") + xlim(0.5, 16) + theme_void()

pa <- pa + labs(tag = "a"); pb <- pb + labs(tag = "b"); pc <- pc + labs(tag = "c"); pd <- pd + labs(tag = "d")
p <- wrap_plots(L = pl, A = pa, B = pb, C = pc, D = pd, design = "LL\nAA\nBB\nCD", heights = c(0.06, 1, 1, 0.95), widths = c(1, 1))
save_panel(pa, "Fig3a_R2", 180, 62); save_panel(pb, "Fig3b_rho_prolif", 180, 62); save_panel(pc, "Fig3c_residual", 90, 70); save_panel(pd, "Fig3d_HR_attenuation", 90, 70)
save_pub(p, "Fig3_proxy", width_mm = 180, height_mm = 185)
cat("FIG3 OK\n")
