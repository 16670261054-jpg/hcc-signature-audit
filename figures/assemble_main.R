# 主图 4–6 拼版（v2，2026-10-18）：只读各脚本保存的面板对象，不碰数据
#   Fig 4  模块审计：a 16 模块哑铃 | b 环形审计全景 | c 逐层过关桑基 | d,e DSigDB 反例
#   Fig 5  一个尺寸匹配随机签名走完整套模板：KM×3 | 时间依赖 ROC | 多因素 Cox（单细胞投影见 Fig 6c）
#   Fig 6  单细胞层：UMAP 图谱/组织来源 | 已发表 vs 随机评分投影 | 相对定位分半小提琴 | 定位热图 | scRNA 阳性率 | Geneformer×3 | scTenifoldKnk
# 运行顺序：fig4_modules.R → fig5_singlecell.R → figX_landscape.R → figA_funnel.R → figCD_sc_extra.R → 本脚本
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
r2 <- wrap_plots(list(tg(F4$pr, "b"), tg(F4$pf, "c")), nrow = 1, widths = c(0.9, 1.25))   # 10-19：删 d（与 Fig 6c 右图重复）
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
