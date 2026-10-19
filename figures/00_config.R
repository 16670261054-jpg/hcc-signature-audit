# ============================================================
# HCC 细胞死亡签名审计 —— 出图全局配置（N28_fig）
# 结构抄自 Paper_AFB1_HCC_nettox/scripts/00_config.R（已验证可跑），只改项目相关部分。
# 口径：pub-plotting SKILL（分析/出图分离、四格式、theme_pub、配色锁死）
#       + 全局规则 10（图上只打星号，精确 P 进 SourceData）。
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

B    <- Sys.getenv("HCC_AUDIT_ROOT", ".")
OUT  <- file.path(B, "N11_out")
FDIR <- file.path(B, "N28_fig")
FIG  <- file.path(FDIR, "figures")
SD   <- file.path(FDIR, "source_data")
DAT  <- file.path(FDIR, "data")
for (d in c(FIG, SD, DAT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
SEED <- 1234

# ---- 配色（全文锁死；2026-10-18 v2：对齐乳酸化项目 Paper_Lactylation_MASH/00_config.R 的 RdBu 体系）----
# 用户反馈 v1 色相太多（红/藏青/琥珀/紫队列条/品红评分/绿闸门…）→ 收成一套：
#   主轴 = RdBu 五档（#2166AC #4393C3 #92C5DE | #F4A582 #D6604D #B2182B）+ 灰；
#   单细胞细胞类型直接沿用乳酸化 COL_CELLTYPE；连续评分一律 inferno；类别 >4 个时用 Set2 柔色。
RDBU <- c(b3 = "#2166AC", b2 = "#4393C3", b1 = "#92C5DE", r1 = "#F4A582", r2 = "#D6604D", r3 = "#B2182B")
COL_GRP <- c("Published"            = RDBU[["r2"]],   # 被检验的主角 = 红（同乳酸化 MASH 红）
             "Size-matched random"  = RDBU[["b2"]],   # 主对照 = 蓝（同乳酸化 Control 蓝）
             "Input-matched random" = RDBU[["b1"]],   # 次对照 = 浅蓝
             "Cell-death gene set"  = RDBU[["r3"]])   # 另一类真实集合 = 深红（只出现在 Fig 2）
ARM2GRP <- c(A_real = "Published", Ak_null = "Size-matched random",
             A_null = "Input-matched random", B_real = "Cell-death gene set",
             B_null = "Input-matched random")
COL_COHORT <- c(TCGA = "grey95", GSE14520 = "grey95", LIRI = "grey95")   # 队列分面条不再着色
COH_LAB <- c(TCGA = "TCGA-LIHC", GSE14520 = "GSE14520", LIRI = "LIRI-JP")
COL_NULL <- "grey82"      # 零分布 / 背景
COL_REF  <- "grey40"; REF_LTY <- 2
COL_PROXY <- RDBU[["b3"]] # 增殖等代理参照线
BRADLEY  <- 0.075         # Bradley 1978 宽松上限
COL_BAND <- "grey92"      # 可接受假阳性带
COL_CELL <- c(Hepatocyte = "#D6604D", Endothelial = "#4393C3", Myeloid = "#F4A582",
              "T/NK" = "#5AAE61", Fibroblast = "#9970AB", B = "#80B1D3")    # = 乳酸化 COL_CELLTYPE 同名同色
COL_SITE <- c(Normal = RDBU[["b1"]], Tumor = RDBU[["r2"]], PVTT = RDBU[["r3"]], Lymph = "#BDBDBD")
COL_RISK <- c("Low risk" = RDBU[["b2"]], "High risk" = RDBU[["r2"]])
COL_YEAR <- c("1" = RDBU[["r1"]], "3" = RDBU[["r2"]], "5" = RDBU[["r3"]])
COL_PASS <- c(`TRUE` = RDBU[["b3"]], `FALSE` = "grey80")
GRAD_HEAT <- c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")   # 比例/AUC 类单向热图（同乳酸化 Fig3c）
SET2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
scale_inferno <- function(...) scale_colour_viridis_c(option = "inferno", ...)

geom_ref_h <- function(y, ...) geom_hline(yintercept = y, linewidth = 0.3, linetype = REF_LTY, colour = COL_REF, ...)
geom_ref_v <- function(x, ...) geom_vline(xintercept = x, linewidth = 0.3, linetype = REF_LTY, colour = COL_REF, ...)

# ---- 星号（全局规则 10；阈值表写进图注）------------------------------------
stars <- function(p) fcase(is.na(p), "", p < 1e-4, "****", p < 1e-3, "***",
                           p < 1e-2, "**", p < 0.05, "*", default = "ns")
STAR_NOTE <- "*P<0.05, **P<0.01, ***P<0.001, ****P<0.0001; ns, not significant."

# ---- 落盘辅助（规则 3）：图上每个数先写 SourceData --------------------------
write_sd <- function(df, name) { p <- file.path(SD, paste0(name, ".tsv"))
  d <- as.data.frame(df)   # 2026-10-18：标签里的换行符在 quote=FALSE 下会把一行拆成两行 → 落盘前换成空格
  for (j in seq_along(d)) if (is.character(d[[j]]) || is.factor(d[[j]])) d[[j]] <- gsub("[\r\n]+", " ", as.character(d[[j]]))
  fwrite(d, p, sep = "\t", quote = FALSE)
  stopifnot(file.exists(p), nrow(fread(p, sep = "\t", quote = "")) == nrow(d)); invisible(p) }

# ---- theme ---------------------------------------------------------------
theme_pub <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(axis.line  = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.title = element_text(size = base_size),
          axis.text  = element_text(size = base_size - 1, colour = "black"),
          plot.title = element_text(size = base_size - 0.5, face = "bold", hjust = 0.5, margin = margin(0, 0, 2, 0)),
          plot.subtitle = element_text(size = base_size - 1.5, colour = "grey35", hjust = 0.5, face = "italic"),
          legend.key.size = unit(3, "mm"),
          legend.text  = element_text(size = base_size - 1),
          legend.title = element_text(size = base_size - 1),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size, face = "bold", margin = margin(1.2, 1.6, 1.2, 1.6)),
          plot.tag = element_text(size = base_size + 1.5, face = "bold"),
          plot.margin = margin(4, 6, 4, 4),
          panel.grid = element_blank())
}
theme_set(theme_pub())
theme_void_pub <- function() theme_void(base_size = 7, base_family = "Arial") +
  theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6), plot.tag = element_text(size = 9, face = "bold"),
        legend.key.size = unit(3, "mm"))

# 分面条着色（ggh4x，代替 gtable 硬编码）
strip_cohort <- function(levels_present) ggh4x::strip_vanilla()   # v2：队列分面条统一浅灰，不着色
## 单面板导出（乳酸化项目做法：每个面板单独成图，拼版另出）
PANEL <- file.path(FIG, "panels"); dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)
save_panel <- function(p, name, w = 90, h = 70) save_pub(p, name, width_mm = w, height_mm = h, dir = PANEL)

# ---- 设备：本机 ragg/Cairo 全挂 → grDevices cairo ---------------------------
.bitmap_dev <- function(path, kind, w, h, dpi) {
  if (kind == "tiff") grDevices::tiff(path, width = w, height = h, units = "in", res = dpi,
                                      type = "cairo", family = "Arial", compression = "lzw")
  else grDevices::png(path, width = w, height = h, units = "in", res = dpi, type = "cairo", family = "Arial")
}
save_pub <- function(plot, filename, width_mm = 180, height_mm = 120, dpi = 600, dir = FIG) {
  f <- file.path(dir, filename); w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(f, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  .bitmap_dev(paste0(f, ".tiff"), "tiff", w, h, dpi); print(plot); dev.off()
  .bitmap_dev(paste0(f, ".png"), "png", w, h, 300); print(plot); dev.off()
  for (e in c(".svg", ".pdf", ".tiff", ".png")) stopifnot(file.exists(paste0(f, e)), file.size(paste0(f, e)) > 1000)
  invisible(f)
}

# ---- 手工 sigmoid 流线（2026-10-18：本机 ggalluvial 0.12.6 + ggplot2 4.0.3 的 geom_flow 画不出曲线/静默丢层，实测）----
# links: data.table(x0, x1, from, to, n, fill)；levels: 每个轴上层级从上到下的顺序（字符向量，所有轴共用）；w = 柱宽
sankey_polys <- function(links, levels, w = 0.3, npt = 40) {
  L <- copy(links); L[, from := factor(from, levels)][, to := factor(to, levels)]
  out <- list(); k <- 0
  for (ax in unique(L$x0)) {
    d <- L[x0 == ax]; tot <- sum(d$n)
    d <- d[order(from, to)]; d[, `:=`(y0_hi = tot - cumsum(n) + n, y0_lo = tot - cumsum(n))]
    d <- d[order(to, from)]; d[, `:=`(y1_hi = tot - cumsum(n) + n, y1_lo = tot - cumsum(n))]
    for (i in seq_len(nrow(d))) {
      t <- seq(0, 1, length.out = npt); s <- (1 - cos(pi * t)) / 2
      xa <- d$x0[i] + w / 2; xb <- d$x1[i] - w / 2; xx <- xa + (xb - xa) * t
      top <- d$y0_hi[i] + (d$y1_hi[i] - d$y0_hi[i]) * s; bot <- d$y0_lo[i] + (d$y1_lo[i] - d$y0_lo[i]) * s
      k <- k + 1
      out[[k]] <- data.table(id = k, x = c(xx, rev(xx)), y = c(top, rev(bot)), fill = d$fill[i])
    }
  }
  rbindlist(out)
}
