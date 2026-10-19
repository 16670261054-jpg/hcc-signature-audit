# Fig 7: two-gate audit with pooled external cohorts
#  a framework schematic (pooled version)  b selected null simulation: positive rate of each criterion (including the missing-OXPHOS scenario)  c power curves (selected true signal)  d 148 published signatures through the gates
# Reads: N11_out/ring4v2/sim_power.tsv, N11_out/n2_gates/published_gates_summary.tsv; the earlier fig6_framework.R is kept
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
SP <- fread(file.path(OUT, "ring4v2/sim_power.tsv")); PG <- fread(file.path(OUT, "n2_gates/published_gates_summary.tsv"))
SV <- fread(file.path(OUT, "ring4/sim_verdict.tsv"))
write_sd(SP, "Fig7bc_selected_sim_power"); write_sd(PG, "Fig7d_published_gates")

## a schematic
bx <- data.table(x = c(1, 2.35, 2.35, 3.7), y = c(0, 0.62, -0.62, 0),
                 lab = c("Signature score\n(fixed before\nexternal data)",
                         "Gate 1\nHarrell's C, both external cohorts\npooled, vs 100 size-matched random\nsignatures (empirical P)",
                         "Gate 2\nStratified Cox LRT over stage, age, sex,\nproliferation, stroma, immune",
                         "Pass both\ngates"),
                 fill = c("grey95", "#D1E5F0", "#D1E5F0", "#FDDBC7"))
ar <- data.table(x = c(1.32, 1.32, 2.98, 2.98), xend = c(1.72, 1.72, 3.42, 3.42), y = c(0.12, -0.12, 0.5, -0.5), yend = c(0.5, -0.5, 0.12, -0.12))
pa <- ggplot() +
  geom_segment(data = ar, aes(x, y, xend = xend, yend = yend), linewidth = 0.4, arrow = arrow(length = unit(1.3, "mm"), type = "closed")) +
  geom_label(data = bx, aes(x, y, label = lab, fill = fill), size = 1.95, family = "Arial", lineheight = 0.92,
             linewidth = 0.25, label.padding = unit(1.4, "mm"), label.r = unit(1, "mm")) +
  scale_fill_identity() + xlim(0.45, 4.25) + ylim(-1, 1) + theme_void_pub()

## b selected null: positive rates
CRL <- c(T0 = "Template criterion", T2p = "Gate 2 alone", T1p = "Gate 1 alone", FWp = "Both gates", FW = "Both gates,\ncohorts separately")
nb <- SP[beta == 0 & criterion %in% names(CRL)]
nb[, crit := factor(CRL[criterion], levels = rev(CRL))][, scen_l := factor(fifelse(scen == "selected", "Selected null", "Null with omitted\nOXPHOS effect"), levels = c("Selected null", "Null with omitted\nOXPHOS effect"))]
COL_SC <- c("Selected null" = RDBU[["b2"]], "Null with omitted\nOXPHOS effect" = RDBU[["r1"]])
pb <- ggplot(nb, aes(rate, crit, colour = scen_l)) +
  annotate("rect", xmin = 0.025, xmax = BRADLEY, ymin = -Inf, ymax = Inf, fill = COL_BAND) +
  geom_linerange(aes(xmin = CP_lo, xmax = CP_hi), linewidth = 0.7, position = position_dodge(0.55)) +
  geom_point(aes(fill = scen_l), shape = 21, size = 2, colour = "black", stroke = 0.25, position = position_dodge(0.55)) +
  scale_colour_manual(NULL, values = COL_SC) + scale_fill_manual(NULL, values = COL_SC) +
  scale_x_continuous("Null signatures declared prognostic", labels = scales::percent, breaks = c(0, 0.075, 0.2, 0.4)) +
  scale_y_discrete(NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 5.5), axis.text.y = element_text(size = 6))

## c power
COL_CR <- c("Template criterion" = "grey60", "Both gates" = RDBU[["r2"]], "Both gates,\ncohorts separately" = RDBU[["b1"]])
pw <- SP[scen == "selected" & criterion %in% c("T0", "FWp", "FW")][, crit := factor(CRL[criterion], levels = names(COL_CR))]
adj <- c(SV$published_adjHR_median_GSE, SV$published_adjHR_median_LIRI)
pc <- ggplot(pw, aes(HR_per_SD, rate, colour = crit, fill = crit)) +
  annotate("rect", xmin = min(adj), xmax = max(adj), ymin = -Inf, ymax = Inf, fill = "#FDDBC7", alpha = 0.8) +
  geom_ref_h(0.8) +
  geom_ribbon(aes(ymin = CP_lo, ymax = CP_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.55) + geom_point(shape = 21, size = 1.5, colour = "white", stroke = 0.3) +
  annotate("text", x = max(pw$HR_per_SD), y = 0.8, label = "80%", hjust = 1, vjust = -0.5, size = 1.8, colour = COL_REF) +
  scale_colour_manual(NULL, values = COL_CR) + scale_fill_manual(NULL, values = COL_CR) +
  scale_x_continuous("Hazard ratio per SD of a true 5-gene signal", breaks = round(sort(unique(pw$HR_per_SD)), 2)) +
  scale_y_continuous("Signatures declared prognostic", labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 5.5))

## d published signatures
DL <- c(template_T0 = "Template criterion", gate2_pooled = "Gate 2", gate1_pooled = "Gate 1", both_pooled = "Both gates")
pg <- PG[criterion %in% names(DL)][, lab := factor(DL[criterion], levels = rev(DL))]
null_fwp <- SP[scen == "selected" & beta == 0 & criterion == "FWp"]
pd <- ggplot(pg, aes(rate, lab)) +
  annotate("rect", xmin = null_fwp$CP_lo, xmax = null_fwp$CP_hi, ymin = 0.5, ymax = 1.5, fill = "#D1E5F0") +
  geom_col(fill = COL_GRP[["Published"]], width = 0.6, colour = "black", linewidth = 0.2) +
  geom_errorbarh(aes(xmin = CP_lo, xmax = CP_hi), height = 0.18, linewidth = 0.3) +
  geom_text(aes(x = CP_hi + 0.03, label = sprintf("%d/%d", positive, n)), hjust = 0, size = 1.9, family = "Arial") +
  scale_x_continuous("Published signatures passing (n = 148)", labels = scales::percent, limits = c(0, 0.78), expand = c(0, 0)) +
  scale_y_discrete(NULL) + theme(axis.text.y = element_text(size = 6))

p <- wrap_plots(list(A = pa + labs(tag = "a"), B = pb + labs(tag = "b"), C = pc + labs(tag = "c"), D = pd + labs(tag = "d")),
                design = "AAA\nBCD", heights = c(0.55, 1), widths = c(1.05, 1.05, 0.9))
save_panel(pb, "Fig7b_selected_null", 70, 70); save_panel(pc, "Fig7c_power", 75, 70); save_panel(pd, "Fig7d_published_gates", 65, 60)
save_pub(p, "Fig7_framework", width_mm = 180, height_mm = 120)
cat("FIG7 v2 OK\n")
