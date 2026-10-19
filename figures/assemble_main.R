# Assembly of main Figs. 4-6: reads only the panel objects saved by the panel scripts, no data
#   Fig 4  module audit: a 16-module dumbbell | b circular audit landscape | c stage-by-stage Sankey | d,e DSigDB counterexample
#   Fig 5  one size-matched random signature through the full template: KM x3 | time-dependent ROC | multivariable Cox (single-cell projection in Fig 6c)
#   Fig 6  single-cell level: UMAP atlas/tissue origin | published vs random score projection | relative-localisation split violin | localisation heatmap | scRNA positive rate | Geneformer x3 | scTenifoldKnk
# Run order: fig4_modules.R -> fig5_singlecell.R -> figX_landscape.R -> figA_funnel.R -> figCD_sc_extra.R -> this script
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "figures/00_config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggtree); library(ggtreeExtra); library(ggnewscale) })
F4 <- readRDS(file.path(DAT, "obj_fig4.rds")); SCo <- readRDS(file.path(DAT, "obj_sc.rds"))
pL1 <- readRDS(file.path(DAT, "panel_landscape.rds")); pA <- readRDS(file.path(DAT, "panel_funnel.rds"))
pC <- readRDS(file.path(DAT, "panel_sc_violin.rds")); pD <- readRDS(file.path(DAT, "panel_gf_deciles.rds"))
tg <- function(p, t) p + labs(tag = t)

## ---- Fig 4 ----
W <- function(p) wrap_elements(full = p)
pL1s <- pL1 + labs(title = NULL) +
  theme(legend.position = "right", plot.caption = element_text(size = 4.6), legend.key.size = unit(2.4, "mm"),
        legend.text = element_text(size = 5.2), legend.title = element_text(size = 5.6), plot.tag = element_text(size = 9, face = "bold", family = "Arial"))
de <- wrap_plots(list(tg(F4$pd1, "d"), tg(F4$pd2, "e")), ncol = 1)
p4 <- wrap_plots(list(F4$lg,
                      wrap_plots(list(W(tg(F4$pa, "a")), W(tg(pA, "b"))), widths = c(1.25, 1)),
                      wrap_plots(list(W(tg(pL1s, "c")), W(de)), widths = c(1.6, 0.75))),
                 ncol = 1, heights = c(0.035, 1, 1.12))
save_pub(p4, "Fig4_template_modules", width_mm = 180, height_mm = 215)

## ---- Fig 5 ----
pk <- F4$pk
r1 <- wrap_plots(list(tg(pk[[1]], "a"), pk[[2]], pk[[3]]), nrow = 1)
r2 <- wrap_plots(list(tg(F4$pr, "b"), tg(F4$pf, "c")), nrow = 1, widths = c(0.9, 1.25))   # panel d removed (duplicate of the right part of Fig 6c)
p5 <- wrap_plots(list(F4$klg, wrap_elements(full = r1), wrap_elements(full = r2)), ncol = 1, heights = c(0.07, 1, 1.15))
save_pub(p5, "Fig5_random_signature_showcase", width_mm = 180, height_mm = 125)

## ---- Fig 6 ----
pcc <- (tg(SCo$pc1, "c") | SCo$pc2)
pg3 <- wrap_plots(list(tg(SCo$pf, "g"), pD + labs(title = NULL, subtitle = NULL) + scale_x_discrete("Abundance decile", labels = 1:10) +
                    geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")), vjust = -0.4, size = 1.8) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)), SCo$pf2),
                  widths = c(1.25, 1.1, 0.5))
p6 <- wrap_plots(list(a = tg(SCo$pa, "a"), b = tg(SCo$pb, "b"), c = pcc, d = tg(SCo$pe, "d"),
                      v = tg(pC, "e"), r = tg(SCo$pd, "f"), g = pg3, k = tg(SCo$pg, "h"), l = SCo$pl),
                 design = "aaaaaabbbbbb\naaaaaabbbbbb\nccccccccdddd\nvvvvvvvvrrrr\ngggggggkkkkk\nllllllllllll",
                 heights = c(1, 1, 1.35, 0.95, 1.15, 0.09))
save_pub(p6, "Fig6_singlecell_modules", width_mm = 180, height_mm = 255)
cat("ASSEMBLE OK\n")
