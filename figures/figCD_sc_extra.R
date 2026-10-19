# Additional panels C / D (descriptive, no new tests)
#  C single cell: mean score of each signature in 6 cell types, z-scored across cell types within each signature ("relative localisation"), published vs size-matched random (first Ak replicate), split violin
#  D Geneformer: proportion of significant genes by decile of detection count (m10_abund_deciles.tsv)
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(ggdist))

## ---- C ----
U <- fread(file.path(DAT, "sc_umap.tsv.gz"), select = c("cell", "celltype"))
SC <- readRDS(file.path(DAT, "sc_sigscore.rds")); stopifnot(identical(rownames(SC), U$cell))
keep <- grep("^(A\\|.*\\|real\\|0|Ak\\|.*\\|null\\|1)$", colnames(SC), value = TRUE)
ct <- factor(U$celltype, levels = names(COL_CELL))
M <- rowsum(SC[, keep], ct, na.rm = TRUE) / as.vector(table(ct))
Z <- scale(M)                                      # z-score each signature (column) across cell types
D <- as.data.table(as.table(Z)); setnames(D, c("celltype", "sid", "z"))
D <- D[!is.na(z)][, grp := fifelse(grepl("^A\\|", sid), "Published", "Size-matched random")]
D[, celltype := factor(celltype, levels = names(COL_CELL))]
stopifnot(D[, uniqueN(sid), by = grp]$V1 >= 140)
write_sd(D, "Fig6e_sc_relative_localisation")
write_sd(D[, .(n_sig = uniqueN(sid), median_z = median(z)), by = .(grp, celltype)][order(celltype, grp)], "Fig6e_sc_relative_localisation_summary")
D[, xn := as.numeric(celltype)]
set.seed(SEED); D[, xj := xn + fifelse(grp == "Published", -0.14, 0.14) + runif(.N, -0.05, 0.05)]
pC <- ggplot(D, aes(y = z, fill = grp)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  stat_slab(data = D[grp == "Published"], aes(x = xn), side = "left", scale = 0.42, alpha = 0.85, colour = NA,
            normalize = "xy", position = position_nudge(x = -0.03)) +
  stat_slab(data = D[grp != "Published"], aes(x = xn), side = "right", scale = 0.42, alpha = 0.85, colour = NA,
            normalize = "xy", position = position_nudge(x = 0.03)) +
  geom_point(aes(x = xj, colour = grp), size = 0.25, alpha = 0.35, stroke = 0) +
  geom_boxplot(aes(x = xn, group = interaction(celltype, grp)), width = 0.1, outlier.shape = NA, linewidth = 0.25,
               position = position_dodge(0.3), fill = "white") +
  geom_tile(data = unique(D[, .(xn, celltype)]), aes(x = xn, y = -2.35, fill = NULL), fill = COL_CELL[levels(D$celltype)],
            width = 0.7, height = 0.12, inherit.aes = FALSE) +
  scale_fill_manual(NULL, values = COL_GRP[1:2]) + scale_colour_manual(values = COL_GRP[1:2], guide = "none") +
  scale_x_continuous(NULL, breaks = seq_along(levels(D$celltype)), labels = levels(D$celltype), expand = c(0.04, 0.04)) +
  scale_y_continuous("Relative signature score\n(z across cell types)", breaks = -2:2) +
  labs(title = NULL) +
  theme(legend.position = "top", legend.key.size = unit(2.6, "mm"), axis.text.x = element_text(size = 6, colour = "black"))
save_panel(pC, "FigC_sc_split_violin", 110, 70)
saveRDS(pC, file.path(DAT, "panel_sc_violin.rds"))

## ---- D ----
DE <- fread(file.path(OUT, "hiend_m10/m10_abund_deciles.tsv"))
DE[, lab := sprintf("%d-%d", NDet_min, NDet_max)][, lab := factor(lab, levels = lab)]
DE[, share := n_sig / sum(n_sig)]
write_sd(DE, "Fig6g_geneformer_deciles")
pD <- ggplot(DE, aes(lab, frac_sig)) +
  geom_col(aes(fill = frac_sig), width = 0.72, colour = "black", linewidth = 0.2) +
  scale_fill_gradientn(colours = GRAD_HEAT, guide = "none") +
  scale_y_continuous("Genes significant after deletion\n(FDR < 0.05)", labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete("Detections in 1,000 tumour hepatocytes (decile range)") +
  labs(title = NULL,
       subtitle = sprintf("Top decile holds %.0f%% of all significant genes", 100 * DE[NDet_decile == 10]$share)) +
  theme(axis.text.x = element_text(size = 5.5, angle = 35, hjust = 1))
save_panel(pD, "FigD_geneformer_deciles", 85, 65)
saveRDS(pD, file.path(DAT, "panel_gf_deciles.rds"))
cat("CD OK\n")
