# Fig 1  study design and corpus
# Reads: N12_signatures/{signatures_index,module_frequency_fulltext}.tsv, N11_out/robust/cohort_meta.tsv, figures/data/sc_umap.tsv.gz
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
S12 <- file.path(B, "N12_signatures")

## a design schematic (four stages)
CM <- fread(file.path(OUT, "robust/cohort_meta.tsv"))
UM <- fread(file.path(DAT, "sc_umap.tsv.gz"), select = c("patient"))
IX <- fread(file.path(S12, "signatures_index.tsv"), encoding = "UTF-8")
IX <- IX[group_primary == "是" & 是否有完整系数 == "是"]
inp <- data.table(y = c(3, 2, 1, 0),
  lab = c(sprintf("%d published HCC cell-death\nsignatures (%d scorable)", nrow(IX), 148L),
          "8 cell-death gene sets\n(cuproptosis, ferroptosis, ...)",
          "Random null signatures:\nsize-matched / input-matched",
          sprintf("TCGA-LIHC %d | GSE14520 %d\nLIRI-JP %d | scRNA %s cells",
                  CM[cohort == "TCGA"]$n, CM[cohort == "GSE14520"]$n, CM[cohort == "LIRI"]$n, format(nrow(UM), big.mark = ","))))
COL_RING <- c("#FDDBC7", "#D1E5F0", "#F0F0F0", "#FDDBC7")
ring <- data.table(x = c(2.25, 3.25, 4.25, 5.25), ttl = c("1  Phenomenon", "2  Mechanism", "3  Downstream", "4  Remedy"),
  q = c("Do cell-death sets\nbeat size-matched\nrandom sets?", "How much do\nproliferation,\nstroma, immune,\nOXPHOS explain?",
        "How often do 16\ntemplate modules\ncall random\nsignatures positive?", "Does a two-gate\naudit hold false\npositives at the\nnominal level?"),
  fig = c("Fig. 2", "Fig. 3", "Figs. 4-6", "Fig. 7"), fill = COL_RING)
write_sd(inp, "Fig1a_inputs")
pa <- ggplot() +
  geom_label(data = inp, aes(x = -0.65, y = y * 0.6 + 0.12, label = lab), hjust = 0, size = 1.8, lineheight = 0.9,
             family = "Arial", fill = "grey97", linewidth = 0.2, label.padding = unit(1.1, "mm"), label.r = unit(0.8, "mm")) +
  annotate("segment", x = 1.45, xend = 1.72, y = 1.02, yend = 1.02, linewidth = 0.45, arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_tile(data = ring, aes(x = x, y = 1.02, fill = fill), width = 0.9, height = 1.75, colour = NA) +
  geom_text(data = ring, aes(x = x, y = 1.66, label = ttl), size = 2.1, fontface = "bold", family = "Arial") +
  geom_text(data = ring, aes(x = x, y = 1.05, label = q), size = 1.85, lineheight = 0.9, family = "Arial") +
  geom_text(data = ring, aes(x = x, y = 0.38, label = fig), size = 1.85, fontface = "italic", family = "Arial", colour = "grey30") +
  annotate("segment", x = ring$x[1:3] + 0.455, xend = ring$x[2:4] - 0.455, y = 1.02, yend = 1.02, linewidth = 0.35,
           arrow = arrow(length = unit(0.9, "mm"), type = "closed")) +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.7, 5.75), ylim = c(0.0, 2.05), clip = "off") +
  theme_void_pub() + labs(title = NULL)

## b use of downstream modules in the full texts (modules audited here shown dark)
MF <- fread(file.path(S12, "module_frequency_fulltext.tsv"))
aud <- c("timeROC/ROC", "KM", "CIBERSORT/immune", "GSEA/GSVA", "checkpoint", "nomogram", "consensus clustering/subtype",
         "calibration", "TMB/mutation", "drug sens (oncoPredict/pRRophetic)", "single-cell", "ML multi-algorithm", "DCA")
MF[, audited := module %in% aud]
lab_b <- c("timeROC/ROC" = "Time-dependent ROC", KM = "Kaplan-Meier", "CIBERSORT/immune" = "Immune infiltration",
           "GSEA/GSVA" = "GSEA / GSVA", checkpoint = "Immune checkpoints", nomogram = "Nomogram",
           "qPCR/IHC wet" = "qPCR / IHC (wet lab)", "consensus clustering/subtype" = "Consensus clustering",
           calibration = "Calibration curve", "TMB/mutation" = "TMB / mutation", "drug sens (oncoPredict/pRRophetic)" = "Drug sensitivity",
           "ceRNA/lncRNA" = "ceRNA / lncRNA", methylation = "Methylation", "single-cell" = "Single-cell RNA-seq",
           "TIDE/IPS" = "TIDE / IPS", CNV = "Copy number", HPA = "Human Protein Atlas", "ML multi-algorithm" = "Multi-algorithm ML",
           DCA = "Decision curve", spatial = "Spatial transcriptomics", docking = "Molecular docking", CellChat = "CellChat",
           pseudotime = "Pseudotime", MR = "Mendelian randomisation")
MF[, lab := lab_b[module]]; stopifnot(!anyNA(MF$lab))
MF[, lab := factor(lab, levels = MF[order(frac)]$lab)]
write_sd(MF[, .(module, lab, n_papers, n_fulltext, frac, audited)], "Fig1b_module_frequency")
pb <- ggplot(MF, aes(frac, lab)) +
  geom_segment(aes(x = 0, xend = frac, yend = lab, colour = audited), linewidth = 0.55) +
  geom_point(aes(fill = audited), shape = 21, size = 1.8, stroke = 0.25) +
  geom_text(aes(x = frac + 0.03, label = n_papers), hjust = 0, size = 1.7, family = "Arial", colour = "grey30") +
  scale_colour_manual(values = c(`TRUE` = RDBU[["b3"]], `FALSE` = "grey80"), guide = "none") +
  scale_fill_manual(NULL, values = c(`TRUE` = RDBU[["b3"]], `FALSE` = "grey88"), labels = c(`TRUE` = "Audited here", `FALSE` = "Not audited")) +
  scale_x_continuous(sprintf("Share of %d full-text papers", MF$n_fulltext[1]), labels = scales::percent, limits = c(0, 1.1),
                     breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_discrete(NULL) + labs(title = NULL) +
  theme(axis.text.y = element_text(size = 5.6), legend.position = "inside", legend.position.inside = c(0.98, 0.04),
        legend.justification = c(1, 0), legend.key.size = unit(2.4, "mm"), legend.background = element_blank())

## c signatures by year and mechanism (stacked streams)
mech_map <- function(m) fcase(m == "pyroptosis", "Pyroptosis", m == "cuproptosis", "Cuproptosis", m == "ferroptosis", "Ferroptosis",
                              m == "necroptosis", "Necroptosis", m == "anoikis", "Anoikis", m == "autophagy", "Autophagy",
                              m == "disulfidptosis", "Disulfidptosis", default = "Multiple / other")
IX[, mech := mech_map(机制)]
YM <- IX[, .N, by = .(year, mech)]
YM <- CJ(year = sort(unique(YM$year)), mech = unique(YM$mech))[YM, on = .(year, mech), N := i.N][is.na(N), N := 0L]
write_sd(YM, "Fig1c_signatures_by_year_mechanism")
COL_MECH <- setNames(SET2, c("Pyroptosis", "Cuproptosis", "Ferroptosis", "Necroptosis", "Anoikis", "Autophagy", "Disulfidptosis", "Multiple / other"))
YM[, mech := factor(mech, levels = names(COL_MECH))]
pc <- ggplot(as.data.frame(YM), aes(year, N, alluvium = mech, stratum = mech, fill = mech)) +
  ggalluvial::geom_flow(width = 0.55, alpha = 0.35, knot.pos = 0.3, colour = NA) +
  ggalluvial::geom_stratum(width = 0.55, colour = "white", linewidth = 0.2) +
  scale_fill_manual(NULL, values = COL_MECH) +
  scale_x_continuous(NULL, breaks = sort(unique(YM$year)), labels = function(x) ifelse(x == max(YM$year), paste0(x, "*"), x)) +
  scale_y_continuous("Signatures published", expand = expansion(mult = c(0, 0.05))) +
  labs(title = NULL) +
  theme(legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.5))

p <- wrap_plots(list(pa + labs(tag = "a"), pb + labs(tag = "b"), pc + labs(tag = "c")), design = "AA\nBC", heights = c(0.6, 1), widths = c(1, 1.15))
save_panel(pa, "Fig1a_design", 180, 60); save_panel(pb, "Fig1b_module_frequency", 85, 85); save_panel(pc, "Fig1c_signatures_by_year", 100, 75)
save_pub(p, "Fig1_design_corpus", width_mm = 180, height_mm = 130)
cat("FIG1 OK\n")
