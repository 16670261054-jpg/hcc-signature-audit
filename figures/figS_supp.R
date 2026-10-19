# 补充图 S1–S4（v3 版式；只读已落盘结果，不新增检验）
#  S1 敏感性：三队列全覆盖的 95 个签名重复模块审计（a 哑铃）；b 148 vs 95 口径下已发表阳性率一致性
#  S2 M14 共识聚类：亚型与增殖的分离度（a 雨云）与调整增殖前后阳性率（b）——探索性，预注册 N25
#  S3 M8 单细胞"验证"：模板判据下阳性率随签名基因数 k 的变化（观测 vs 单基因阳性率推得的期望 1-(1-p)^k）
#  S4 基因复用：前 25 基因棒棒糖 + 前 20 基因 × 死亡类型弦图（弦图为 base 图，按 pub-plotting 规则以位图嵌入）
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggdist) })
cp <- function(x, n) as.numeric(binom.test(x, n)$conf.int)

## ---------------- S1 ----------------
lab <- fread(file.path(SD, "Fig4a_module_rates.tsv"))[, .(module, label)]
r3 <- fread(file.path(OUT, "ring3_armA_decision_v3_rescaled.tsv"))[analysis == "fitsig_pool100_k|fullcoverage"]
r3 <- r3[, .(module = sub("_pos", "", module), rate_real, n_real, rate_null = FPR_null, n_null, P = fisher_P)]
pick <- function(f, mods, cmp = "Afull_vs_Ak_sizematched") fread(f)[comparison == cmp & module %in% mods,
  .(module = sub("_pos$", "", module), rate_real, n_real, rate_null = FPR_null, n_null, P = fisher_P)]
M95 <- rbind(r3, pick(file.path(OUT, "hiend_m3_m12/m3_m12_decision.tsv"), c("M3_pos", "M12_pos")),
             pick(file.path(OUT, "hiend_m8/m8_decision.tsv"), c("M8a_pos", "M8b_pos")),
             pick(file.path(OUT, "hiend_m13_m15/m13_m15_decision.tsv"), c("M13_pos", "M14_pos", "DCA_pos")),
             fread(file.path(OUT, "hiend_m10/m10_decision.tsv"))[comparison == "Afull_vs_Ak_sizematched" & filter == "all",
               .(module = "M10", rate_real, n_real, rate_null = FPR_null, n_null, P = fisher_P)])
M95[module == "DCA", module := "M15"]
M95 <- merge(M95, lab, by = "module"); stopifnot(nrow(M95) == 16, all(M95$n_real == 95))
M95[, c("lo", "hi") := as.list(cp(round(rate_null * n_null), n_null)), by = module][, star := stars(P)]
write_sd(M95, "FigS1a_fullcoverage_module_rates")
M95[, label := factor(label, levels = M95[order(rate_null)]$label)]
s1a <- ggplot(M95, aes(y = label)) +
  annotate("rect", xmin = 0, xmax = BRADLEY, ymin = -Inf, ymax = Inf, fill = COL_BAND) +
  geom_segment(aes(x = rate_real, xend = rate_null, yend = label), colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = COL_GRP[["Size-matched random"]], linewidth = 1.4, alpha = 0.35) +
  geom_point(aes(x = rate_null), shape = 21, size = 2.1, fill = COL_GRP[["Size-matched random"]], stroke = 0.25) +
  geom_point(aes(x = rate_real), shape = 21, size = 1.8, fill = COL_GRP[["Published"]], stroke = 0.25) +
  geom_text(aes(x = 1.1, label = star), size = 2) +
  scale_x_continuous("Signatures called positive", labels = scales::percent, breaks = seq(0, 1, 0.25), limits = c(0, 1.15), expand = c(0, 0)) +
  scale_y_discrete(NULL) + coord_cartesian(clip = "off") +
  labs(title = "Full-coverage subset (n = 95)") + theme(axis.text.y = element_text(size = 5.8))
M148 <- fread(file.path(SD, "Fig4a_module_rates.tsv"))[, .(module, r148 = rate_real, P148 = P)]
CC <- merge(M148, M95[, .(module, r95 = rate_real, P95 = P)], by = "module")
CC[, agree := fifelse((P148 < 0.05) == (P95 < 0.05), "Same verdict", "Verdict differs")]
write_sd(CC, "FigS1b_148_vs_95_consistency")
s1b <- ggplot(CC, aes(r148, r95)) +
  geom_abline(slope = 1, intercept = 0, linetype = REF_LTY, colour = COL_REF, linewidth = 0.3) +
  geom_point(aes(fill = agree), shape = 21, size = 2, stroke = 0.25) +
  ggrepel::geom_text_repel(data = CC[agree != "Same verdict"], aes(label = module), size = 2, family = "Arial", seed = SEED) +
  scale_fill_manual(NULL, values = c("Same verdict" = "grey75", "Verdict differs" = RDBU[["r2"]])) +
  scale_x_continuous("Published positive rate, all 148", labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous("Published positive rate, full-coverage 95", labels = scales::percent, limits = c(0, 1)) +
  coord_equal() + theme(legend.position = c(0.05, 0.95), legend.justification = c(0, 1), legend.key.size = unit(2.6, "mm"))
lgS1 <- ggplot(data.table(x = c(1, 3.4), g = c("Published", "Size-matched random")), aes(x, 1, fill = g)) +
  geom_point(shape = 21, size = 2.1, stroke = 0.25) + geom_text(aes(x = x + 0.15, label = g), hjust = 0, size = 2.1, family = "Arial") +
  scale_fill_manual(values = COL_GRP, guide = "none") + xlim(0.8, 6.5) + theme_void()
pS1 <- wrap_plots(list(lgS1, wrap_plots(list(s1a + labs(tag = "a"), s1b + labs(tag = "b")), widths = c(1.5, 1))), ncol = 1, heights = c(0.06, 1))
save_pub(pS1, "FigS1_fullcoverage_sensitivity", width_mm = 180, height_mm = 95)

## ---------------- S2 ----------------
P14 <- fread(file.path(OUT, "hiend_m14_prolif/m14_prolif_rows.tsv"))
P14 <- P14[arm %in% c("A", "Ak")][, grp := factor(fifelse(arm == "A", "Published", "Size-matched random"), levels = names(COL_GRP)[1:2])]
D14 <- fread(file.path(OUT, "hiend_m14_prolif/m14_prolif_decision.tsv"))
write_sd(P14, "FigS2a_m14_prolif_rows"); write_sd(D14, "FigS2b_m14_prolif_decision")
set.seed(SEED); P14[, xj := as.numeric(grp) - 0.2 + runif(.N, -0.06, 0.06)]
s2a <- ggplot(P14, aes(grp, prolif_sep)) +
  geom_point(aes(x = xj, colour = grp), size = 0.35, alpha = 0.4, stroke = 0) +
  stat_slab(aes(fill = grp), side = "right", scale = 0.5, alpha = 0.85, colour = NA, density = "bounded", position = position_nudge(x = 0.03)) +
  geom_boxplot(aes(fill = grp), width = 0.12, outlier.shape = NA, linewidth = 0.25, position = position_nudge(x = -0.05)) +
  annotate("segment", x = 1, xend = 2, y = 1.04, yend = 1.04, linewidth = 0.3) +
  annotate("text", x = 1.5, y = 1.06, label = stars(D14$MWU_P), size = 2.2, vjust = 0) +
  scale_fill_manual(values = COL_GRP, guide = "none") + scale_colour_manual(values = COL_GRP, guide = "none") +
  scale_x_discrete(NULL, labels = c("Published", "Size-matched\nrandom")) +
  scale_y_continuous("Subtype separation by proliferation\n(|AUC - 0.5| x 2)", limits = c(0, 1.12), breaks = seq(0, 1, 0.25))
b2 <- data.table(state = rep(c("Unadjusted", "Adjusted for\nproliferation"), each = 2), grp = rep(names(COL_GRP)[1:2], 2),
                 rate = c(D14$raw_rate_pub, D14$raw_rate_Ak, D14$adj_rate_pub, D14$adj_rate_Ak),
                 n = rep(c(D14$n_pub, D14$n_Ak), 2))
b2[, c("lo", "hi") := as.list(cp(round(rate * n), n)), by = .(state, grp)]
b2[, state := factor(state, levels = unique(state))][, grp := factor(grp, levels = names(COL_GRP)[1:2])]
st2 <- data.table(state = factor(levels(b2$state), levels = levels(b2$state)), lab = stars(c(D14$raw_fisher_P, D14$adj_fisher_P)))
s2b <- ggplot(b2, aes(state, rate, fill = grp)) +
  geom_col(position = position_dodge(0.7), width = 0.62, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), width = 0.15, linewidth = 0.3) +
  geom_text(data = st2, aes(state, 1.04, label = lab), inherit.aes = FALSE, size = 2.2) +
  scale_fill_manual(NULL, values = COL_GRP[1:2]) +
  scale_y_continuous("Subtypes differ in survival", labels = scales::percent, limits = c(0, 1.08), expand = c(0, 0)) +
  scale_x_discrete(NULL) + theme(legend.position = "top", legend.key.size = unit(2.6, "mm"))
pS2 <- wrap_plots(list(s2a + labs(tag = "a"), s2b + labs(tag = "b")), widths = c(1, 1.1))
save_pub(pS2, "FigS2_consensus_proliferation", width_mm = 140, height_mm = 70)

## ---------------- S3 ----------------
PS <- fread(file.path(OUT, "hiend_m8/m8_per_signature.tsv"))[arm %in% c("Ak", "A") & kind == "null"]
PS[, grp := fifelse(arm == "Ak", "Size-matched random", "Input-matched random")]
OB <- PS[k <= 12, .(n = .N, a = sum(M8a_pos), b = sum(M8b_pos)), by = .(grp, k)][n >= 20]
OB <- melt(OB, id.vars = c("grp", "k", "n"), measure.vars = c("a", "b"), variable.name = "crit", value.name = "pos")
OB[, rate := pos / n][, c("lo", "hi") := as.list(cp(pos, n)), by = .(grp, k, crit)]
GL <- fread(file.path(OUT, "hiend_m8/m8_gene_level_rates.tsv"))
EX <- CJ(k = 1:12, crit = c("a", "b"))[, p := fifelse(crit == "a", GL$p_gene_M8a_in_U, GL$p_gene_M8b_in_U)][, rate := 1 - (1 - p)^k]
CRL <- c(a = "Cell-type marker (M8a)", b = "Tumour vs normal (M8b)")
OB[, crit_l := factor(CRL[as.character(crit)], levels = CRL)]; EX[, crit_l := factor(CRL[crit], levels = CRL)]
write_sd(OB, "FigS3_m8_rate_by_k_observed"); write_sd(EX, "FigS3_m8_rate_by_k_expected")
pS3 <- ggplot(OB, aes(k, rate)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BRADLEY, fill = COL_BAND) +
  geom_line(data = EX, colour = "grey30", linewidth = 0.4, linetype = "22") +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = grp), linewidth = 0.5, position = position_dodge(0.4)) +
  geom_point(aes(fill = grp), shape = 21, size = 1.7, stroke = 0.25, position = position_dodge(0.4)) +
  facet_wrap(~ crit_l, nrow = 1) +
  scale_fill_manual(NULL, values = COL_GRP[2:3]) + scale_colour_manual(values = COL_GRP[2:3], guide = "none") +
  scale_x_continuous("Genes in signature (k)", breaks = seq(1, 12, 1)) +
  scale_y_continuous("Signatures called positive", labels = scales::percent, limits = c(0, 1)) +
  theme(legend.position = "top", legend.key.size = unit(2.6, "mm"))
save_pub(pS3, "FigS3_singlecell_k_dependence", width_mm = 150, height_mm = 70)

## ---------------- S4 ----------------
lol <- readRDS(file.path(DAT, "panel_lollipop.rds"))
img <- png::readPNG(file.path(PANEL, "FigL2_gene_reuse_chord.png"))
chord <- wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
pS4 <- wrap_plots(list(lol + labs(tag = "a"), chord + labs(tag = "b")), widths = c(1, 1.15))
save_pub(pS4, "FigS4_gene_reuse", width_mm = 180, height_mm = 95)
cat("SUPP OK\n")
