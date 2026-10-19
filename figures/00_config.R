# ============================================================
# HCC cell-death signature audit: global figure settings
# Shared settings for all figure scripts.
# Conventions: analysis and plotting separated, multiple export formats, theme_pub, fixed palette
#       significance shown as asterisks only; exact P values go to Source Data.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

B    <- Sys.getenv("HCC_AUDIT_ROOT", ".")
OUT  <- file.path(B, "N11_out")
FDIR <- file.path(B, "N28_fig")
FIG  <- file.path(FDIR, "figures")
SD   <- file.path(FDIR, "source_data")
DAT  <- file.path(B, "figures", "data")
for (d in c(FIG, SD, DAT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
SEED <- 1234

# ---- Palette (fixed throughout) ----
# One palette:
#   main axis = five RdBu levels (#2166AC #4393C3 #92C5DE | #F4A582 #D6604D #B2182B) + grey;
#   cell types use COL_CELLTYPE; continuous scores use inferno; Set2 when there are more than 4 categories.
RDBU <- c(b3 = "#2166AC", b2 = "#4393C3", b1 = "#92C5DE", r1 = "#F4A582", r2 = "#D6604D", r3 = "#B2182B")
COL_GRP <- c("Published"            = RDBU[["r2"]],   # tested signatures = red
             "Size-matched random"  = RDBU[["b2"]],   # main control = blue
             "Input-matched random" = RDBU[["b1"]],   # secondary control = light blue
             "Cell-death gene set"  = RDBU[["r3"]])   # other real sets = dark red (Fig. 2 only)
ARM2GRP <- c(A_real = "Published", Ak_null = "Size-matched random",
             A_null = "Input-matched random", B_real = "Cell-death gene set",
             B_null = "Input-matched random")
COL_COHORT <- c(TCGA = "grey95", GSE14520 = "grey95", LIRI = "grey95")   # cohort facet strips not coloured
COH_LAB <- c(TCGA = "TCGA-LIHC", GSE14520 = "GSE14520", LIRI = "LIRI-JP")
COL_NULL <- "grey82"      # null distribution / background
COL_REF  <- "grey40"; REF_LTY <- 2
COL_PROXY <- RDBU[["b3"]] # reference lines for proliferation and other proxies
BRADLEY  <- 0.075         # Bradley 1978 liberal upper limit
COL_BAND <- "grey92"      # acceptable false-positive band
COL_CELL <- c(Hepatocyte = "#D6604D", Endothelial = "#4393C3", Myeloid = "#F4A582",
              "T/NK" = "#5AAE61", Fibroblast = "#9970AB", B = "#80B1D3")
COL_SITE <- c(Normal = RDBU[["b1"]], Tumor = RDBU[["r2"]], PVTT = RDBU[["r3"]], Lymph = "#BDBDBD")
COL_RISK <- c("Low risk" = RDBU[["b2"]], "High risk" = RDBU[["r2"]])
COL_YEAR <- c("1" = RDBU[["r1"]], "3" = RDBU[["r2"]], "5" = RDBU[["r3"]])
COL_PASS <- c(`TRUE` = RDBU[["b3"]], `FALSE` = "grey80")
GRAD_HEAT <- c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")   # one-directional heatmaps for proportions/AUC
SET2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
scale_inferno <- function(...) scale_colour_viridis_c(option = "inferno", ...)

geom_ref_h <- function(y, ...) geom_hline(yintercept = y, linewidth = 0.3, linetype = REF_LTY, colour = COL_REF, ...)
geom_ref_v <- function(x, ...) geom_vline(xintercept = x, linewidth = 0.3, linetype = REF_LTY, colour = COL_REF, ...)

# ---- Asterisks (threshold table given in the legends) ------------------------------------
stars <- function(p) fcase(is.na(p), "", p < 1e-4, "****", p < 1e-3, "***",
                           p < 1e-2, "**", p < 0.05, "*", default = "ns")
STAR_NOTE <- "*P<0.05, **P<0.01, ***P<0.001, ****P<0.0001; ns, not significant."

# ---- Writing helper: every number shown in a figure is first written to Source Data --------------------------
write_sd <- function(df, name) { p <- file.path(SD, paste0(name, ".tsv"))
  d <- as.data.frame(df)   # newlines in labels would split a row under quote=FALSE -> replaced by spaces before writing
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

# Coloured facet strips (ggh4x, instead of hard-coded gtable)
strip_cohort <- function(levels_present) ggh4x::strip_vanilla()   # cohort facet strips uniformly light grey, not coloured
## Single-panel export (each panel saved separately; assembly done elsewhere)
PANEL <- file.path(FIG, "panels"); dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)
save_panel <- function(p, name, w = 90, h = 70) save_pub(p, name, width_mm = w, height_mm = h, dir = PANEL)

# ---- Device: grDevices cairo (ragg/Cairo unavailable on the analysis machine) ---------------------------
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

# ---- Hand-drawn sigmoid flows (geom_flow in ggalluvial 0.12.6 + ggplot2 4.0.3 failed to draw curves / silently dropped layers) ----
# links: data.table(x0, x1, from, to, n, fill); levels: top-to-bottom order of levels on each axis (character vector shared by all axes); w = bar width
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
