# Earlier version of the two-gate audit figure: design, type I error and power of the simulated positive control, one published signature as example
# Reads: N11_out/ring4/{sim_power,sim_verdict,nault_5gene}.tsv
# Statistics: positive rates with Clopper-Pearson 95% CI (read from disk); criteria as defined in ring4_positive_control.R.
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
R4 <- file.path(OUT, "ring4")
SP <- fread(file.path(R4, "sim_power.tsv")); SV <- fread(file.path(R4, "sim_verdict.tsv")); NA5 <- fread(file.path(R4, "nault_5gene.tsv"))
write_sd(SP, "Fig7b_sim_power"); write_sd(SV, "Fig7b_sim_verdict"); write_sd(NA5, "Fig7c_example_5gene")

COL_CRIT <- c(T0_template_uni = "grey70", T1_random_null = RDBU[["b1"]], T2_proxy_LRT = RDBU[["b2"]], framework_T1andT2 = RDBU[["r2"]])
LAB_CRIT <- c(T0_template_uni = "Template: univariable Cox P < 0.05", T1_random_null = "Gate 1: beats size-matched random",
              T2_proxy_LRT = "Gate 2: survives proxy adjustment", framework_T1andT2 = "Framework (gate 1 and gate 2)")

## a framework schematic
bx <- data.table(x = c(1, 2.35, 2.35, 3.7), y = c(0, 0.62, -0.62, 0),
                 lab = c("Signature score\n(fixed before\nexternal data)",
                         "Gate 1\nEmpirical P vs 100 size-matched\nrandom signatures (Phipson-Smyth)",
                         "Gate 2\nNested Cox LRT over stage, age, sex,\nproliferation, stroma, immune",
                         "Pass both gates\nin both external\ncohorts"),
                 fill = c("grey95", "#D1E5F0", "#D1E5F0", "#FDDBC7"))
ar <- data.table(x = c(1.32, 1.32, 2.92, 2.92), xend = c(1.78, 1.78, 3.38, 3.38), y = c(0.12, -0.12, 0.5, -0.5), yend = c(0.5, -0.5, 0.12, -0.12))
pa <- ggplot() +
  geom_segment(data = ar, aes(x, y, xend = xend, yend = yend), linewidth = 0.4, arrow = arrow(length = unit(1.3, "mm"), type = "closed")) +
  geom_label(data = bx, aes(x, y, label = lab, fill = fill), size = 1.95, family = "Arial", lineheight = 0.92,
             linewidth = 0.25, label.padding = unit(1.4, "mm"), label.r = unit(1, "mm")) +
  scale_fill_identity() + xlim(0.45, 4.25) + ylim(-1, 1) + theme_void_pub() +
  labs(title = NULL)

## b simulation: type I error and power
SP[, crit := factor(criterion, levels = names(COL_CRIT))]
adj <- c(SV$published_adjHR_median_GSE, SV$published_adjHR_median_LIRI)
pb <- ggplot(SP, aes(HR_per_SD, rate, colour = crit, fill = crit)) +
  annotate("rect", xmin = min(adj), xmax = max(adj), ymin = -Inf, ymax = Inf, fill = "#FDDBC7", alpha = 0.8) +
  annotate("text", x = max(adj) + 0.012, y = 0.97, label = "Published signatures,\nmedian adjusted HR", size = 1.75, lineheight = 0.85, hjust = 0, vjust = 1, family = "Arial") +
  geom_ref_h(0.8) + geom_ref_h(BRADLEY) +
  geom_ribbon(aes(ymin = CP_lo, ymax = CP_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(shape = 21, size = 1.5, colour = "white", stroke = 0.3) +
  annotate("text", x = 1.835, y = 0.8, label = "80% power", hjust = 1, vjust = -0.5, size = 1.8, colour = COL_REF) +
  annotate("text", x = 1.835, y = BRADLEY, label = "7.5%", hjust = 1, vjust = -0.5, size = 1.8, colour = COL_REF) +
  scale_colour_manual(NULL, values = COL_CRIT, labels = LAB_CRIT) + scale_fill_manual(NULL, values = COL_CRIT, labels = LAB_CRIT) +
  scale_x_continuous("Simulated hazard ratio per SD of a true 5-gene signal", breaks = round(sort(unique(SP$HR_per_SD)), 2)) +
  scale_y_continuous("Signatures declared prognostic", labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(title = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 5.5),
        legend.margin = margin(0, 0, 0, 0)) + guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

## c example: a published 5-gene HCC score (equal weights; coefficients not published in the article)
ex <- NA5[set == "main_TAF9"]
cc <- data.table(gate = rep(c("Template", "Gate 1", "Gate 2"), each = 2), cohort = rep(c("GSE14520", "LIRI"), 3),
                 P = c(NA, NA, ex$empP_GSE, ex$empP_LIRI, ex$lrt_GSE, ex$lrt_LIRI))
cc[gate == "Template", pass := ex$T0]
cc[gate != "Template", pass := P < 0.05]
cc[, lab := fifelse(gate == "Template", fifelse(pass, "pass", "fail"), paste0(fifelse(pass, "pass", "fail"), "\n", stars(P)))]
cc[, `:=`(gate = factor(gate, levels = c("Template", "Gate 1", "Gate 2")), cohort_l = factor(COH_LAB[cohort], levels = COH_LAB[c("GSE14520", "LIRI")]))]
write_sd(cc, "Fig7c_gate_table")
pc <- ggplot(cc, aes(cohort_l, gate)) +
  geom_tile(aes(fill = pass), colour = "white", linewidth = 0.8, width = 0.92, height = 0.88) +
  geom_text(aes(label = lab), size = 1.9, lineheight = 0.85, family = "Arial", colour = "white", fontface = "bold") +
  scale_fill_manual(values = COL_PASS, guide = "none") +
  scale_y_discrete(NULL, limits = rev) + scale_x_discrete(NULL, position = "top") +
  labs(title = "Published 5-gene score",
       subtitle = sprintf("Proxy R\u00b2 = %.2f (TCGA)\nFramework: %s", ex$R2_proxy_TCGA, fifelse(ex$T1 & ex$T2, "pass", "not distinguishable from random"))) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 6))

p <- wrap_plots(list(pa + labs(tag = "a"), pb + labs(tag = "b"), pc + labs(tag = "c")), design = "AA\nBC", heights = c(0.62, 1), widths = c(1.45, 1))
save_panel(pa, "Fig6a_framework", 180, 55); save_panel(pb, "Fig6b_power", 105, 85); save_panel(pc, "Fig6c_example", 80, 72)
save_pub(p, "Fig7_framework", width_mm = 180, height_mm = 120)
cat("FIG6 OK\n")
