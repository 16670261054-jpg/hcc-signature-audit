# 补图 S5：Fig 5 示例随机签名（Ak|S004|rep1）在全部 3,700 个尺寸匹配随机签名中的位置
# 只读盘：data/fig5_showcase_all_Ak.tsv、data/fig5_showcase_typicality.tsv（check_fig5_showcase_typicality.R）
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
R  <- fread(file.path(DAT, "fig5_showcase_all_Ak.tsv")); stopifnot(nrow(R) == 3700)
TY <- fread(file.path(DAT, "fig5_showcase_typicality.tsv"))
COH <- c(TCGA = "TCGA-LIHC", GSE14520 = "GSE14520", LIRI = "LIRI-JP")
L <- rbind(
  melt(R, id.vars = c("src", "rep"), measure.vars = paste0("HR_", names(COH)), variable.name = "v")[, metric := "KM"],
  melt(R, id.vars = c("src", "rep"), measure.vars = paste0("C_", names(COH)), variable.name = "v")[, metric := "C"],
  R[, .(src, rep, v = "multiHR_TCGA", value = multiHR, metric = "MV")])
L[, cohort := factor(COH[sub("^(HR|C|multiHR)_", "", v)], levels = COH)]
EXR <- L[src == "S004" & rep == 1]
EXR[, key := fifelse(metric == "MV", "multiHR", as.character(v))][, pct := TY$pct_rank[match(key, TY$metric)]]
EXR[, lab := sprintf("%.0fth percentile", 100 * pct)][, hj := fifelse(pct > 0.75, 1.05, -0.05)]
stopifnot(nrow(EXR) == 7, !anyNA(EXR$pct))
write_sd(L[, .(src, rep, metric, cohort, value)], "FigS5_random_signature_distribution")
write_sd(EXR[, .(metric, cohort, value, percentile = pct)], "FigS5_example_position")
write_sd(TY, "FigS5_typicality_summary")
COL_EX <- RDBU[["r3"]]
panel <- function(m, xlab, logx = FALSE) {
  d <- L[metric == m]; e <- EXR[metric == m]
  p <- ggplot(d, aes(value)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    geom_histogram(bins = 45, fill = COL_GRP[["Size-matched random"]], colour = "white", linewidth = 0.1, alpha = 0.9) +
    geom_vline(data = e, aes(xintercept = value), colour = COL_EX, linewidth = 0.5) +
    geom_label(data = e, aes(x = value, y = Inf, label = lab, hjust = hj), colour = COL_EX, size = 2, vjust = 1.3,
               fill = "white", label.size = 0, label.padding = unit(0.6, "mm")) +
    ggh4x::facet_wrap2(~ cohort, nrow = 1, scales = "free_x", strip = strip_cohort()) +
    labs(x = xlab, y = "Random signatures")
  if (logx) p <- p + scale_x_log10() + geom_ref_v(1)
  if (m == "C") p <- p + geom_ref_v(0.5)
  p
}
pa <- panel("KM", "Hazard ratio, high vs low score (log scale)", logx = TRUE)
pb <- panel("C", "C-index")
pc <- panel("MV", "Multivariable hazard ratio per s.d. (log scale)", logx = TRUE)
pS5 <- (pa / pb / (pc + plot_spacer() + plot_layout(widths = c(1, 2)))) + plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "a")
save_pub(pS5, "FigS5_showcase_typicality", width_mm = 180, height_mm = 150)
cat("FIGS5 OK\n")
