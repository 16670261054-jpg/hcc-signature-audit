# Additional panels: descriptive rearrangements of results on disk, no new tests
#  L1 circular audit-landscape heatmap: 148 published signatures x 15 template modules (positive/negative), clustered by positivity profile; outer ring = cell-death type (ggtree fan + geom_fruit)
#  L2 gene-reuse chord diagram: genes most often used by signatures of different cell-death types <-> cell-death type; gene sectors marked by membership of Hallmark E2F/G2M (proliferation)
#  L3 lollipop of gene usage counts (top 25), coloured by proliferation-gene membership
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggtree); library(ggtreeExtra); library(ggnewscale); library(circlize) })
S12 <- file.path(B, "N12_signatures")

## ---- Per-signature module results (148 published) ----
r6 <- fread(file.path(OUT, "ring3_armA_realrescale/results_m06_real_rescaled.tsv"))[arm == "A" & kind == "real"]
# M0 requires direction: C > 0.5 and two-sided Noether P < 0.05 in both external cohorts (survcomp default is two-sided)
fixM0 <- function(d) { if (all(c("M0_c_GSE", "M0_c_LIRI", "M0_p_GSE", "M0_p_LIRI") %in% names(d)))
  d[, M0_pos := as.integer(M0_p_GSE < 0.05 & M0_p_LIRI < 0.05 & M0_c_GSE > 0.5 & M0_c_LIRI > 0.5)]; d }
r6 <- fixM0(r6)
r7 <- fread(file.path(OUT, "ring3_full/results_m7.tsv"))[arm == "A" & kind == "real", .(src, M7_pos)]
m3 <- fread(file.path(OUT, "hiend_m3_m12/m3_m12_per_signature.tsv"))[arm == "A" & kind == "real", .(src, M3_pos, M12_pos)]
m8 <- fread(file.path(OUT, "hiend_m8/m8_per_signature.tsv"))[arm == "A" & kind == "real", .(src, M8a_pos, M8b_pos)]
m10 <- fread(file.path(OUT, "hiend_m10/m10_per_signature.tsv"))[arm == "A" & kind == "real", .(src, M10_pos)]
m13 <- fread(file.path(OUT, "hiend_m13_m15/m13_m15_per_signature.tsv"))[arm == "A" & kind == "real", .(src, M13_pos, M14_pos, M15_pos = DCA_pos)]
W <- Reduce(function(a, b) merge(a, b, by = "src", all = TRUE),
            list(r6[, .(src, k, M0_pos, M1_pos, M2_pos, M4_pos, M5a_pos, M5b_pos, M6_pos)], r7, m3, m8, m10, m13))
stopifnot(nrow(W) == 148)
mods <- setdiff(names(W), c("src", "k"))
## Consistency check: per-signature positive rates must equal the published positive rates of Fig 4a on disk (signatures skipped in M14 kept as NA as in per_signature)
F4 <- fread(file.path(SD, "Fig4a_module_rates.tsv"))
chk <- data.table(module = sub("_pos", "", mods), rate_here = vapply(mods, function(m) mean(W[[m]], na.rm = TRUE), 1))
chk <- merge(chk, F4[, .(module, rate_real)], by = "module")
write_sd(chk, "Fig4c_landscape_consistency_check")
stopifnot(max(abs(chk$rate_here - chk$rate_real)) < 1e-9)

IX <- fread(file.path(S12, "signatures_index.tsv"), encoding = "UTF-8")
mech_map <- function(m) fcase(m == "pyroptosis", "Pyroptosis", m == "cuproptosis", "Cuproptosis", m == "ferroptosis", "Ferroptosis",
                              m == "necroptosis", "Necroptosis", m == "anoikis", "Anoikis", m == "autophagy", "Autophagy",
                              m == "disulfidptosis", "Disulfidptosis", default = "Multiple / other")
COL_MECH <- setNames(SET2, c("Pyroptosis", "Cuproptosis", "Ferroptosis", "Necroptosis", "Anoikis", "Autophagy", "Disulfidptosis", "Multiple / other"))
W <- merge(W, IX[, .(src = sig_id, mech = mech_map(机制), year)], by = "src")
W[, n_pos := rowSums(.SD, na.rm = TRUE), .SDcols = mods]
write_sd(W, "Fig4c_landscape_per_signature_modules")

## ---- L1 circular heatmap ---- (ring numbers and bar-axis ticks drawn vertically to avoid overlap)
X <- as.matrix(W[, ..mods]); X[is.na(X)] <- 0; rownames(X) <- W$src
hc <- hclust(dist(X, method = "binary"), method = "ward.D2")
F4[, mlab := paste0(module, "  ", label)]
mod_ord <- F4[order(-rate_null)]$module            # inner -> outer ring: random positive rate from high to low
long <- melt(W[, c("src", mods), with = FALSE], id.vars = "src", variable.name = "module", value.name = "pos")
SHORT <- c(M0 = "External val.", M1 = "Indep. Cox", M2 = "Immune infil.", M3 = "Drug sens.", M4 = "GSEA", M5a = "TMB diff.", M5b = "TMB x risk",
           M6 = "Nomogram", M7 = "Multi-ML", M8a = "scRNA marker", M8b = "scRNA T vs N", M10 = "Geneformer", M12 = "DSigDB",
           M13 = "Checkpoint", M14 = "Consensus", M15 = "DCA")
RING <- setNames(as.character(seq_along(mod_ord)), mod_ord)
long[, module := factor(RING[sub("_pos", "", module)], levels = RING)]
ring_key <- paste(sprintf("%s %s", RING, SHORT[mod_ord]), collapse = "  |  ")
long[, call := fcase(is.na(pos), "Not run", pos == 1, "Positive", default = "Negative")]
tr <- ggtree(hc, layout = "fan", open.angle = 42, branch.length = "none", linewidth = 0.2)
pL1 <- tr +
  geom_fruit(data = long, geom = geom_tile, mapping = aes(y = src, x = module, fill = call),
             offset = 0.06, pwidth = 1.25, colour = "white", linewidth = 0.12,
             axis.params = list(axis = "x", text.angle = 90, text.size = 1.6, hjust = 1, vjust = 0.5, colour = "grey15")) +
  scale_fill_manual("Template call", values = c(Positive = RDBU[["r2"]], Negative = "#D1E5F0", "Not run" = "grey92")) +
  new_scale_fill() +
  geom_fruit(data = W, geom = geom_tile, mapping = aes(y = src, x = 1, fill = mech), offset = 0.04, pwidth = 0.07, colour = NA) +
  scale_fill_manual("Cell-death type", values = COL_MECH) +
  new_scale_fill() +
  geom_fruit(data = W, geom = geom_col, mapping = aes(y = src, x = n_pos), fill = "grey55", colour = NA,
             offset = 0.04, pwidth = 0.22, orientation = "y",
             axis.params = list(axis = "x", text.size = 1.5, text.angle = 90, hjust = 1, vjust = 0.5, nbreak = 3, colour = "grey30", title = "No. positive", title.size = 1.6, title.height = 0.02),
             grid.params = list(linetype = 3, linewidth = 0.15)) +
  theme(legend.position = "right", legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
        legend.key.size = unit(2.6, "mm"), plot.title = element_text(size = 7.5, face = "bold")) +
  labs(title = "Audit landscape: 148 published signatures x 16 template modules",
       caption = paste0("Rings, inner to outer (random-signature positive rate high to low):
", stringr::str_wrap(ring_key, 150))) +
  theme(plot.caption = element_text(size = 5.2, hjust = 0, colour = "grey25"))
save_panel(pL1, "FigL1_audit_landscape_circular", 150, 135); saveRDS(pL1, file.path(DAT, "panel_landscape.rds"))
write_sd(F4[, .(module, label, rate_null)][order(-rate_null)], "Fig4c_landscape_ring_order")

## ---- L2 / L3 gene reuse ----
L <- fread(file.path(S12, "signatures_long.tsv"), encoding = "UTF-8")
L <- merge(L[sig_id %in% W$src, .(src = sig_id, gene = gene_hgnc)], W[, .(src, mech)], by = "src")
L <- unique(L)
rg <- function(f) { s <- strsplit(readLines(f, warn = FALSE), "\t"); setNames(lapply(s, function(x) x[-(1:2)]), vapply(s, `[`, "", 1)) }
H <- rg(file.path(B, "N11_geneset/h.all.v2024.1.Hs.symbols.gmt"))
prolif <- union(H$HALLMARK_E2F_TARGETS, H$HALLMARK_G2M_CHECKPOINT)
GF <- L[, .(n_sig = uniqueN(src), n_mech = uniqueN(mech)), by = gene][order(-n_sig, gene)]
GF[, prolif := gene %in% prolif]
write_sd(GF, "FigS4a_gene_frequency")
top <- GF[1:25]
top[, gene := factor(gene, levels = rev(gene))]
pL3 <- ggplot(top, aes(n_sig, gene)) +
  geom_segment(aes(x = 0, xend = n_sig, yend = gene, colour = prolif), linewidth = 0.6) +
  geom_point(aes(fill = prolif, size = n_mech), shape = 21, stroke = 0.25) +
  scale_colour_manual(values = c(`TRUE` = RDBU[["r2"]], `FALSE` = "grey70"), guide = "none") +
  scale_fill_manual("Hallmark E2F / G2M", values = c(`TRUE` = RDBU[["r2"]], `FALSE` = "grey85"), labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
  scale_size_continuous("Cell-death types\nusing the gene", range = c(1.2, 3.2), breaks = c(1, 3, 5)) +
  scale_x_continuous("Published signatures containing the gene", expand = expansion(mult = c(0, 0.08))) +
  scale_y_discrete(NULL) + labs(title = NULL) +
  theme(axis.text.y = element_text(face = "italic", size = 6), legend.key.size = unit(2.6, "mm"))
save_panel(pL3, "FigL3_gene_reuse_lollipop", 85, 95); saveRDS(pL3, file.path(DAT, "panel_lollipop.rds"))

## Chord diagram: top 20 genes x cell-death type
tg <- GF[1:20]$gene
lk <- L[gene %in% tg, .N, by = .(mech, gene)]
write_sd(lk, "FigS4b_chord_links")
mech_in <- names(COL_MECH)[names(COL_MECH) %in% lk$mech]
gene_ord <- GF[gene %in% tg][order(-prolif, -n_sig)]$gene
grid_col <- c(COL_MECH[mech_in], setNames(ifelse(gene_ord %in% prolif, RDBU[["r2"]], "grey70"), gene_ord))
draw_chord <- function() {
  circos.clear(); circos.par(gap.after = c(rep(2, length(mech_in) - 1), 10, rep(1.2, length(gene_ord) - 1), 10), start.degree = 90, canvas.xlim = c(-1.25, 1.25), canvas.ylim = c(-1.25, 1.25))
  chordDiagram(as.data.frame(lk[, .(mech, gene, N)]), order = c(mech_in, gene_ord), grid.col = grid_col,
               transparency = 0.35, annotationTrack = "grid", annotationTrackHeight = 0.04,
               preAllocateTracks = list(track.height = 0.2), link.border = NA)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1, s, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = 0.5, font = if (s %in% gene_ord) 3 else 2, family = "Arial") }, bg.border = NA)
  legend("bottomleft", legend = c("Hallmark E2F / G2M gene", "Other gene"), fill = c(RDBU[["r2"]], "grey70"), border = NA,
         bty = "n", cex = 0.5)
  circos.clear()
}
f <- file.path(PANEL, "FigL2_gene_reuse_chord"); w <- 95 / 25.4
svglite::svglite(paste0(f, ".svg"), width = w, height = w); draw_chord(); dev.off()
grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = w, family = "Arial"); draw_chord(); dev.off()
.bitmap_dev(paste0(f, ".tiff"), "tiff", w, w, 600); draw_chord(); dev.off()
.bitmap_dev(paste0(f, ".png"), "png", w, w, 300); draw_chord(); dev.off()
stopifnot(file.size(paste0(f, ".png")) > 1000)
cat("LANDSCAPE OK; proliferation genes among top 25:", sum(top$prolif), "\n")
