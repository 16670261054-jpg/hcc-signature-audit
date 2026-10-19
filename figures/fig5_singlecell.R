# Single-cell figure: the single-cell and perturbation modules do not separate published from random signatures either (M8 single-cell localisation / M10 Geneformer / M9 scTenifoldKnk)
# Reads: figures/data/{sc_umap,sc_sigscore,sc_sigscore_coverage} (prep_sc_embedding.R)
#         N11_out/hiend_m8/{m8_decision,m8_per_signature,m8a_findallmarkers_all,signature_genes_all_arms}.tsv
#         N11_out/hiend_m10/{m10_abund_gene_table,m10_abund_check2}.tsv
#         N11_out/hiend_m9_hvg/ko_speedtest/{ALAD,TINAGL1}.tsv, equivalence_aligned_summary_ALAD.tsv
# Statistics: proportions by Fisher's exact (read from decision tables); random positive rate Clopper-Pearson 95% CI; the rest descriptive.
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(scattermore); library(ggrepel); library(ggh4x) })
M8 <- file.path(OUT, "hiend_m8"); M10 <- file.path(OUT, "hiend_m10"); M9 <- file.path(OUT, "hiend_m9_hvg")

U <- fread(file.path(DAT, "sc_umap.tsv.gz"))
U[, celltype := factor(celltype, levels = names(COL_CELL))]
stopifnot(!anyNA(U$celltype))
umap_axes <- function(p, d = U) {     # standard axes (UMAP 1 / UMAP 2)
  p + labs(x = "UMAP 1", y = "UMAP 2") + coord_cartesian(clip = "off") + theme_pub() +
    theme(legend.key.size = unit(3, "mm"))
}
## Cell-type labels placed at the density mode on a coarse 25-cell grid
lab_pos <- U[, { gx <- cut(UMAP1, 25); gy <- cut(UMAP2, 25); t <- table(gx, gy); w <- which(t == max(t), arr.ind = TRUE)[1, ]
  sel <- as.integer(gx) == w[1] & as.integer(gy) == w[2]; .(x = median(UMAP1[sel]), y = median(UMAP2[sel]), n = .N) }, by = celltype]
write_sd(U[, .N, by = .(celltype, site)][order(celltype, site)], "Fig6a_cells_by_type_site")

## a cell atlas
pa <- ggplot(U, aes(UMAP1, UMAP2)) +
  geom_scattermore(aes(colour = celltype), pointsize = 1.6, pixels = c(1400, 1400), alpha = 0.85) +
  geom_text_repel(data = lab_pos, aes(x, y, label = celltype), size = 2.3, fontface = "bold", family = "Arial",
                  bg.color = "white", bg.r = 0.12, min.segment.length = Inf, seed = SEED) +
  scale_colour_manual(NULL, values = COL_CELL, labels = setNames(sprintf("%s (%s)", lab_pos$celltype, trimws(format(lab_pos$n, big.mark = ","))), as.character(lab_pos$celltype)),
                      guide = guide_legend(override.aes = list(size = 2.2))) +
  labs(title = NULL)
pa <- umap_axes(pa) + theme(legend.position = "right", legend.key.size = unit(3, "mm"))

## b tissue origin
U[, site := factor(site, levels = names(COL_SITE))]
pb <- ggplot(U[sample(.N)], aes(UMAP1, UMAP2)) +
  geom_scattermore(aes(colour = site), pointsize = 1.2, pixels = c(900, 900), alpha = 0.8) +
  scale_colour_manual(NULL, values = COL_SITE, guide = guide_legend(override.aes = list(size = 2.2), ncol = 1)) +
  labs(title = NULL)
pb <- umap_axes(pb) + theme(legend.position = "right", legend.justification = "center",
                             legend.key.spacing.x = unit(0.5, "mm"), legend.text = element_text(size = 5.5))

## c one published signature vs its size-matched random signature: per-cell scores
## Example rule (fixed): published k = median k of published signatures, published signature and its first Ak replicate both M8a-positive, 100% single-cell coverage; lowest src identifier taken
PS <- fread(file.path(M8, "m8_per_signature.tsv")); CV <- fread(file.path(DAT, "sc_sigscore_coverage.tsv"))
kmed <- median(PS[arm == "A" & kind == "real"]$k)
cand <- merge(PS[arm == "A" & kind == "real" & k == kmed & M8a_pos == 1, .(src)],
              PS[arm == "Ak" & rep == 1 & M8a_pos == 1, .(src)], by = "src")
cand <- cand[src %in% CV[arm == "A" & k == k_in_sc]$src & src %in% CV[arm == "Ak" & k == k_in_sc]$src][order(src)]
ex <- cand$src[1]; stopifnot(length(ex) == 1, !is.na(ex))
SC <- readRDS(file.path(DAT, "sc_sigscore.rds")); stopifnot(identical(rownames(SC), U$cell))
GA <- fread(file.path(M8, "signature_genes_all_arms.tsv"))
g_pub <- GA[arm == "A" & kind == "real" & src == ex]$gene; g_rnd <- GA[arm == "Ak" & rep == 1 & src == ex]$gene
U[, s_pub := SC[, sprintf("A|%s|real|0", ex)]][, s_rnd := SC[, sprintf("Ak|%s|null|1", ex)]]
write_sd(data.table(example_src = ex, rule = "k == published median, both M8a-positive, full sc coverage, lowest src id",
                    k_median = kmed, n_candidates = nrow(cand), genes_published = paste(g_pub, collapse = ";"),
                    genes_random = paste(g_rnd, collapse = ";")), "Fig6c_example_choice")
write_sd(U[, .(mean_pub = mean(s_pub), mean_rnd = mean(s_rnd)), by = celltype], "Fig6c_example_mean_by_celltype")
feat <- function(v, ttl, genes) {
  d <- U[order(get(v))]; q <- quantile(d[[v]], 0.99); d[, z := pmin(get(v) / q, 1)]
  p <- ggplot(d, aes(UMAP1, UMAP2)) +
    geom_scattermore(data = d[z == 0], colour = "grey88", pointsize = 1.6, pixels = c(1100, 1100)) +
    geom_scattermore(data = d[z > 0], aes(colour = z), pointsize = 1.6, pixels = c(1100, 1100)) +
    scale_colour_viridis_c(option = "inferno", begin = 0.12, end = 0.98, direction = 1, "Signature\nscore", limits = c(0, 1), breaks = c(0, 1), labels = c("low", "high")) +
    labs(title = ttl, subtitle = stringr::str_wrap(paste(genes, collapse = ", "), 34))
  umap_axes(p, d) + theme(legend.key.width = unit(1.6, "mm"), legend.key.height = unit(3, "mm"),
                          plot.subtitle = element_text(size = 5.2, face = "italic", colour = "grey30"))
}
pc1 <- feat("s_pub", "Published signature", g_pub)
pc2 <- feat("s_rnd", "Size-matched random signature", g_rnd)
pc1 <- pc1 + theme(legend.position = "none")

## d M8a / M8b positive rates (template criteria)
D8 <- fread(file.path(M8, "m8_decision.tsv"))
cp <- function(x, n) as.numeric(binom.test(round(x * n), n)$conf.int)
mk_rates <- function(D, mods, modlab) {
  r <- rbind(D[comparison == "A_vs_Ak_sizematched" & module %in% mods, .(module, grp = "Published", rate = rate_real, n = n_real, P = NA_real_)],
             D[comparison == "A_vs_Ak_sizematched" & module %in% mods, .(module, grp = "Size-matched random", rate = FPR_null, n = n_null, P = fisher_P)],
             D[comparison == "A_vs_fitsig_k" & module %in% mods, .(module, grp = "Input-matched random", rate = FPR_null, n = n_null, P = fisher_P)])
  r[, c("lo", "hi") := as.list(cp(rate, n)), by = .(module, grp)]
  r[, `:=`(star = stars(P), grp = factor(grp, levels = names(COL_GRP)[1:3]), mod_l = factor(modlab[module], levels = modlab))]
  r
}
r8 <- mk_rates(D8, c("M8a_pos", "M8b_pos"), c(M8a_pos = "Cell-type\nmarker", M8b_pos = "Tumour vs\nnormal"))
CMH <- fread(file.path(OUT, "cmh_modules.tsv"))   # size-matched comparison by stratified CMH + Holm (same as Fig 4a)
r8[grp == "Size-matched random", `:=`(P = CMH$CMH_P_holm[match(sub("_pos", "", module), CMH$module)], test = "CMH stratified by matched set, Holm across 16 modules (Fig 4a family)")][, star := stars(P)]
## Input-matched comparison: CMH stratified by matched set (no correction, as in the original design; P read from N35_paired_tests); positive rates and CIs unchanged
PT <- fread(file.path(B, "N35_paired_tests/paired_tests_results.tsv"), sep = "	", quote = "")
r8[grp == "Input-matched random", `:=`(P = PT$new_P[match(paste("Fig6f", module, grp, sep = "|"), PT$key)], test = PT$new_test[match(paste("Fig6f", module, grp, sep = "|"), PT$key)])]
stopifnot(!anyNA(r8[grp != "Published"]$P)); r8[, star := stars(P)]
write_sd(r8, "Fig6f_M8_rates")
rate_plot <- function(r, ttl) ggplot(r, aes(mod_l, rate, colour = grp)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BRADLEY, fill = "grey92") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6, position = position_dodge(0.55)) +
  geom_point(aes(fill = grp), shape = 21, size = 2, colour = "black", stroke = 0.25, position = position_dodge(0.55)) +
  geom_text(aes(y = hi, label = star, group = grp), position = position_dodge(0.55), vjust = -0.5, size = 1.9, colour = "black") +
  scale_fill_manual(NULL, values = COL_GRP) + scale_colour_manual(NULL, values = COL_GRP, guide = "none") +
  scale_y_continuous("Signatures called positive", labels = scales::percent, limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(NULL) + labs(title = NULL) +
  theme(legend.position = "none", axis.text.x = element_text(size = 5.8))
pd <- rate_plot(r8, "scRNA 'validation'")

## e cell type where signatures localise (descriptive): proportion of signatures with >= 1 marker gene of that type (same criterion as M8a)
MK <- fread(file.path(M8, "m8a_findallmarkers_all.tsv"))[abs(avg_log2FC) > 1 & p_val_adj < 0.05]
GS <- GA[(arm == "A" & kind == "real") | arm == "Ak"]
GS[, grp := fifelse(arm == "A", "Published", "Size-matched random")]
hit <- merge(GS[, .(grp, sid = paste(arm, src, rep), gene)], MK[, .(gene, celltype = cluster)], by = "gene", allow.cartesian = TRUE)
nsig <- GS[, .(N = uniqueN(paste(arm, src, rep))), by = grp]
E5 <- hit[, .(n = uniqueN(sid)), by = .(grp, celltype)][nsig, on = "grp"][, frac := n / N]
E5 <- CJ(grp = nsig$grp, celltype = names(COL_CELL))[E5, on = .(grp, celltype), `:=`(n = i.n, frac = i.frac)][is.na(frac), `:=`(n = 0L, frac = 0)]
write_sd(E5, "Fig6d_localisation_by_celltype")
E5[, celltype := factor(celltype, levels = rev(names(COL_CELL)))][, grp := factor(grp, levels = names(COL_GRP)[1:2])]
pe <- ggplot(E5, aes(grp, celltype)) +
  geom_tile(aes(fill = frac), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac), colour = frac > 0.45), size = 2, family = "Arial") +
  geom_point(aes(x = 0.38, colour = NULL), data = unique(E5[, .(celltype)]), shape = 15, size = 2.2,
             colour = COL_CELL[as.character(unique(E5$celltype))]) +
  scale_fill_gradientn("Signatures", colours = GRAD_HEAT, limits = c(0, 1), labels = scales::percent) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(NULL, labels = c("Published", "Random"), expand = expansion(add = c(0.75, 0.05))) +
  scale_y_discrete(NULL) + labs(title = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 6),
        legend.key.width = unit(1.6, "mm"), legend.key.height = unit(3, "mm"))

## f Geneformer in silico deletion: significance is determined by detection count (abundance)
GF <- fread(file.path(M10, "m10_abund_gene_table.tsv"))
GF[, sig := as.logical(sig)][, nl := -log10(pmax(Goal_end_FDR, 1e-300))]
lab_g <- GF[sig == TRUE][order(Goal_end_FDR)][1:6]
write_sd(GF[, .(Gene_name, N_Detections, Goal_end_FDR, sig, mean_log1p_cpm, NDet_decile)], "Fig6g_geneformer_genes")
pf <- ggplot(GF, aes(N_Detections, nl)) +
  geom_point(data = GF[sig == FALSE], colour = "grey65", size = 0.35, alpha = 0.3, stroke = 0, position = position_jitter(width = 0, height = 0.08, seed = SEED)) +
  geom_point(data = GF[sig == TRUE], aes(fill = mean_log1p_cpm), shape = 21, size = 1.3, stroke = 0.15) +
  geom_ref_h(-log10(0.05)) +
  geom_text_repel(data = lab_g, aes(label = Gene_name), size = 1.8, fontface = "italic", family = "Arial",
                  segment.size = 0.2, max.overlaps = 30, seed = SEED, box.padding = 0.15) +
  scale_x_log10("Detections (of 1,000 tumour hepatocytes)") +
  scale_y_continuous(expression(-log[10]*"FDR"), trans = scales::pseudo_log_trans(sigma = 0.5, base = 10), breaks = c(0, 1, 3, 10, 30, 100, 200)) +
  scale_fill_gradientn("Mean\nlog1p CPM", colours = GRAD_HEAT[3:5]) +
  labs(title = NULL) +
  theme(legend.key.width = unit(1.4, "mm"), legend.key.height = unit(2.6, "mm"), legend.position = "inside",
        legend.position.inside = c(0.03, 0.97), legend.justification = c(0, 1), legend.background = element_blank())

## f2 after abundance matching
C2 <- fread(file.path(M10, "m10_abund_check2.tsv"))
D10 <- fread(file.path(M10, "m10_decision.tsv"))
f2 <- data.table(grp = c("Published", "Matched\nrandom"), rate = c(C2$rate_real, C2$rate_null),
                 n = c(C2$n_real, C2$n_null))
f2[, c("lo", "hi") := as.list(cp(rate, n)), by = grp][, grp := factor(grp, levels = grp)]
## Between-group test: CMH stratified by matched set (P read from N35_paired_tests); positive rates and CIs unchanged
PT <- fread(file.path(B, "N35_paired_tests/paired_tests_results.tsv"), sep = "	", quote = "")
f2[, P := c(NA, PT[key == "Fig6g|abund_matched"]$new_P)][, test := c(NA, PT[key == "Fig6g|abund_matched"]$new_test)]; stopifnot(!is.na(f2$P[2]))
write_sd(f2, "Fig6g_geneformer_abundance_matched")
pf2 <- ggplot(f2, aes(grp, rate)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BRADLEY, fill = "grey92") +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = grp), linewidth = 0.6) +
  geom_point(aes(fill = grp), shape = 21, size = 2, stroke = 0.25) +
  annotate("text", x = 1.5, y = 0.19, label = stars(f2$P[2]), size = 2.1) +
  annotate("segment", x = 1, xend = 2, y = 0.175, yend = 0.175, linewidth = 0.3) +
  scale_fill_manual(values = c(COL_GRP[["Published"]], COL_GRP[["Size-matched random"]]), guide = "none") +
  scale_colour_manual(values = c(COL_GRP[["Published"]], COL_GRP[["Size-matched random"]]), guide = "none") +
  scale_y_continuous("Signatures called positive", labels = scales::percent, limits = c(0, 0.2)) +
  scale_x_discrete(NULL) + labs(title = NULL) + theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1))

## g scTenifoldKnk: knocking out ALAD and the random gene TINAGL1 gives the same "significant" genes (shared-network version, descriptive)
KO <- rbindlist(lapply(c("ALAD", "TINAGL1"), function(g) fread(file.path(M9, "ko_speedtest", paste0(g, ".tsv")))[, ko := g]))
KO <- KO[gene != ko][order(ko, -Z)][, rank := seq_len(.N), by = ko]
KO[, sig := p.adj < 0.05]
EQ <- fread(file.path(M9, "equivalence_aligned_summary_ALAD.tsv"))
write_sd(KO[sig == TRUE | rank <= 20, .(ko, gene, Z, FC, p.adj, rank, sig)], "Fig6h_knk_top")
write_sd(EQ, "Fig6h_knk_official_equivalence")
KO[, ko_l := factor(paste0(ko, " KO", c(ALAD = " (sig.)", TINAGL1 = " (rand.)")[ko]),
                    levels = c("ALAD KO (sig.)", "TINAGL1 KO (rand.)"))]
pg <- ggplot(KO[rank <= 300], aes(rank, Z)) +
  geom_line(colour = "grey60", linewidth = 0.3) +
  geom_point(data = KO[sig == TRUE], shape = 21, fill = RDBU[["b2"]], size = 1.8, stroke = 0.25) +
  geom_text_repel(data = KO[sig == TRUE], aes(label = gene), size = 2, fontface = "italic", family = "Arial",
                  nudge_x = 90, direction = "y", segment.size = 0.2, seed = SEED) +
  facet_wrap(~ ko_l, nrow = 1) +
  scale_x_continuous("Gene rank (top 300)", breaks = c(1, 150, 300)) +
  scale_y_continuous("Perturbation Z") + labs(title = NULL) +
  theme(strip.text = element_text(size = 6), panel.spacing.x = unit(2, "mm"))

## Group legend
lg <- data.table(grp = factor(names(COL_GRP)[1:3], levels = names(COL_GRP)), x = c(1, 4, 7.4))
pl <- ggplot(lg, aes(x, 1, fill = grp)) + geom_point(shape = 21, size = 2.2, stroke = 0.25) +
  geom_text(aes(x = x + 0.2, label = c("Published", "Size-matched random", "Input-matched random")), hjust = 0, size = 2.1, family = "Arial") +
  annotate("rect", xmin = 10.6, xmax = 11.2, ymin = 0.85, ymax = 1.15, fill = "grey92") +
  annotate("text", x = 11.35, y = 1, label = "Acceptable false-positive range (\u2264 7.5%, Bradley)", hjust = 0, size = 2.1, family = "Arial") +
  scale_fill_manual(values = COL_GRP, guide = "none") + xlim(0.5, 17.5) + ylim(0.7, 1.3) + theme_void()

tg <- function(p, t) p + labs(tag = t)
## patchwork design needs single-character area names -> letter mapping
pdd <- (tg(pc1, "c") | pc2)
pbot <- (tg(pd, "e") | tg(pf, "f") | pf2 | tg(pg, "g")) + plot_layout(widths = c(1.05, 1.9, 0.75, 1.75))
p <- wrap_plots(list(a = tg(pa, "a"), b = tg(pb, "b"), d = pdd, e = tg(pe, "d"), f = pbot, l = pl),
                design = "aaaaaaabbbbb\naaaaaaabbbbb\nddddddddeeee\nffffffffffff\nllllllllllll",
                heights = c(1.1, 1.1, 1.45, 1.25, 0.1))
save_panel(pa, "Fig5a_umap_celltype", 105, 85); save_panel(pb, "Fig5b_umap_site", 95, 85); save_panel(pc1 + theme(legend.position = "right"), "Fig5c1_umap_published_score", 85, 75)
save_panel(pc2, "Fig5c2_umap_random_score", 85, 75); save_panel(pe, "Fig5d_localisation", 75, 72); save_panel(pd, "Fig5e_sc_rates", 55, 68)
save_panel(pf, "Fig5f1_geneformer", 85, 72); save_panel(pf2, "Fig5f2_geneformer_matched", 42, 68); save_panel(pg, "Fig5g_knk", 85, 68)
saveRDS(list(pa = pa, pb = pb, pc1 = pc1, pc2 = pc2, pe = pe, pd = pd, pf = pf, pf2 = pf2, pg = pg, pl = pl), file.path(DAT, "obj_sc.rds"))  # assembled in assemble_main.R
cat("FIG5 OK  example =", ex, "\n")
