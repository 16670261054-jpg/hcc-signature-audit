# Panel: Sankey of the 148 published signatures passing successive criteria (descriptive rearrangement, no new tests)
# All three criteria are per-signature results already on disk, for the same 148 signatures:
#   1 external validation (template M0: C significantly > 0.5 in both external cohorts)      <- ring3_armA_realrescale/results_m06_real_rescaled.tsv
#   2 independent prognostic factor (template M1: TCGA multivariable Cox P < 0.05)       <- same file
#   3 still significant after proxy adjustment (proxy analysis: nested LRT < 0.05 in both external cohorts, base model = stage + proliferation + stromal + immune + age + sex) <- ring2/ring2_rows.tsv
# Each stage is annotated with the positive rate of size-matched random signatures at that stage (from the respective decision tables; random batches differ from the published ones, shown for reference only, not connected)
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
r6 <- fread(file.path(OUT, "ring3_armA_realrescale/results_m06_real_rescaled.tsv"))[arm == "A" & kind == "real"]
# M0 requires direction: C > 0.5 and two-sided Noether P < 0.05 in both external cohorts (survcomp default is two-sided)
fixM0 <- function(d) { if (all(c("M0_c_GSE", "M0_c_LIRI", "M0_p_GSE", "M0_p_LIRI") %in% names(d)))
  d[, M0_pos := as.integer(M0_p_GSE < 0.05 & M0_p_LIRI < 0.05 & M0_c_GSE > 0.5 & M0_c_LIRI > 0.5)]; d }
r6 <- fixM0(r6)[, .(src, M0_pos, M1_pos)]
r2 <- fread(file.path(OUT, "ring2/ring2_rows.tsv"))[arm == "A" & kind == "real" & cohort != "TCGA"]
res <- r2[, .(resid = as.integer(all(lrt_p < 0.05) & .N == 2)), by = src]
W <- merge(r6, res, by = "src"); stopifnot(nrow(W) == 148)
## Item-by-item check against the summaries on disk
RS <- fread(file.path(OUT, "ring2/ring2_residual_summary.tsv"))[comparison == "A_vs_Ak_sizematched" & base_model != "+oxphos"]
stopifnot(sum(W$resid) == RS$pos_real)
D3 <- fread(file.path(OUT, "ring3_armA_decision_v3_rescaled.tsv"))[analysis == "fitsig_pool100_k|all"]
stopifnot(abs(mean(W$M0_pos) - D3[module == "M0_pos"]$rate_real) < 1e-9, abs(mean(W$M1_pos) - D3[module == "M1_pos"]$rate_real) < 1e-9)
ST <- c("External\nvalidation", "Independent\nprognostic factor", "Survives proxy\nadjustment")
lab <- function(x) factor(fifelse(x == 1, "Pass", "Fail"), levels = c("Pass", "Fail"))
W[, `:=`(s1 = lab(M0_pos), s2 = lab(M1_pos), s3 = lab(resid))]
fl <- W[, .N, by = .(s1, s2, s3)]
write_sd(W, "Fig4b_funnel_per_signature"); write_sd(fl, "Fig4b_funnel_flows")
rnd <- data.table(stage = factor(ST, levels = ST), rate = c(D3[module == "M0_pos"]$FPR_null, D3[module == "M1_pos"]$FPR_null, RS$rate_null),
                  pub = c(mean(W$M0_pos), mean(W$M1_pos), mean(W$resid)))
write_sd(rnd, "Fig4b_stage_rates")
## Stage order follows the reporting order of template articles: independent prognosis -> external validation -> proxy adjustment
ST <- c("Independent\nprognostic factor", "External\nvalidation", "Survives proxy\nadjustment")
rnd[, stage := factor(c(ST[2], ST[1], ST[3]), levels = ST)]
fl[, path := fcase(s2 == "Pass" & s1 == "Pass" & s3 == "Pass", "All three", s2 == "Pass" & s1 == "Pass", "Prognostic + external", default = "Other")]
COL_PATH <- c("All three" = RDBU[["r3"]], "Prognostic + external" = RDBU[["r1"]], "Other" = "grey80")
COL_ST <- c(Pass = COL_GRP[["Published"]], Fail = "#D1E5F0")
## Nodes (Pass above, Fail below at each stage) and flows (hand-drawn sigmoid, see sankey_polys in 00_config.R)
ax <- list(s2 = 1, s1 = 2, s3 = 3)
nodes <- rbindlist(lapply(names(ax), function(v) fl[, .(n = sum(N)), by = .(st = get(v))][, x := ax[[v]]]))
nodes[, st := factor(st, levels = c("Pass", "Fail"))]; setorder(nodes, x, st)
nodes[, `:=`(ymax = 148 - cumsum(n) + n, ymin = 148 - cumsum(n)), by = x]
lk <- rbind(fl[, .(x0 = 1, x1 = 2, from = s2, to = s1, n = N, fill = path)], fl[, .(x0 = 2, x1 = 3, from = s1, to = s3, n = N, fill = path)])
lk <- lk[, .(n = sum(n)), by = .(x0, x1, from, to, fill)]
## Keep path order stable when splitting by path within the same from -> to
lk[, fill := factor(fill, levels = names(COL_PATH))]; setorder(lk, x0, from, to, fill)
PP <- sankey_polys(lk, c("Pass", "Fail"), w = 0.3)
pA <- ggplot() +
  geom_polygon(data = PP, aes(x, y, group = id, fill = fill), alpha = 0.6, colour = NA) +
  geom_rect(data = nodes, aes(xmin = x - 0.15, xmax = x + 0.15, ymin = ymin, ymax = ymax, fill = st), colour = "white", linewidth = 0.3) +
  geom_text(data = nodes[n >= 10], aes(x, (ymin + ymax) / 2, label = paste0(st, "\n", n)), size = 2, lineheight = 0.85,
            family = "Arial", colour = "white", fontface = "bold") +
  geom_text(data = nodes[n < 10], aes(x + 0.18, (ymin + ymax) / 2, label = paste(st, n)), size = 1.8, family = "Arial", colour = "grey15", hjust = 0) +
  geom_text(data = rnd, aes(x = as.numeric(stage), y = 162, label = sprintf("Published %.0f%%\nRandom %.0f%%", 100 * pub, 100 * rate)),
            size = 1.85, lineheight = 0.9, family = "Arial", colour = "grey20") +
  scale_fill_manual(NULL, values = c(COL_ST, COL_PATH), breaks = names(COL_PATH),
                    labels = c("All three" = "Passes all three", "Prognostic + external" = "Passes first two only", "Other" = "Other paths")) +
  scale_x_continuous(NULL, breaks = 1:3, labels = ST, expand = c(0.12, 0.12)) +
  scale_y_continuous("Published signatures", breaks = c(0, 50, 100, 148), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 172), clip = "off") +
  labs(title = NULL) +
  theme(axis.text.x = element_text(size = 6.2, colour = "black"), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
save_panel(pA, "FigA_funnel_sankey", 95, 78)
saveRDS(pA, file.path(DAT, "panel_funnel.rds"))
cat("FUNNEL OK", sum(W$M0_pos), sum(W$M1_pos), sum(W$resid), "\n")
