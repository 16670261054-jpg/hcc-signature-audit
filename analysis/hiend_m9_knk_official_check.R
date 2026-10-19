# M9 等价性验证：官方 scTenifoldKnk() 原样调用 1 个基因，与 HVG 等价版（hiend_m9_knk_hvg.R）逐位比对
# 判据（不设近似阈值）：WT 邻接矩阵与 diffRegulation 全部列 identical()；若不 identical，报告最大绝对差并判"不等价"
# 输入与模板一致：as.data.frame(countMat[unique(c(target, hvgs)), ])；靶基因在 HVG 内 → 与等价版输入相同
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(scTenifoldKnk) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out/hiend_m9_hvg")
KNK_CORES <- as.integer(Sys.getenv("KNK_CORES", "4")); KODIR <- file.path(O, Sys.getenv("KNK_KODIR", "ko"))
g <- Sys.getenv("KNK_CHECK_GENE", fread(file.path(KODIR, "_ko_list_used.tsv"))$gene[1])
inp <- readRDS(file.path(O, "input_hvg5000.rds")); hvgs <- readLines(file.path(O, "hvgs.txt")); stopifnot(g %in% hvgs)
data <- as.data.frame(inp[unique(c(g, hvgs)), ]); rm(inp); gc()
t0 <- Sys.time()
set.seed(123)
res <- scTenifoldKnk(countMatrix = data, gKO = g, qc_mtThreshold = 0.1, qc_minLSize = 500, nc_nNet = 10, nc_nCells = 500, nCores = KNK_CORES)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(res$diffRegulation, file.path(O, sprintf("official_%s_diffRegulation.rds", g)))
WTs <- readRDS(file.path(O, "WT_network.rds")); WTo <- as.matrix(res$tensorNetworks$WT)
DRs <- as.data.frame(readRDS(file.path(KODIR, paste0(g, ".rds")))); DRo <- as.data.frame(res$diffRegulation)
num <- intersect(names(DRs)[sapply(DRs, is.numeric)], names(DRo))
out <- data.table(gene = g, official_minutes = mins,
  WT_dim_equal = identical(dim(WTs), dim(WTo)), WT_names_equal = identical(dimnames(WTs), dimnames(WTo)),
  WT_identical = isTRUE(all(WTs == WTo)), WT_max_abs_diff = if (identical(dim(WTs), dim(WTo))) max(abs(WTs - WTo)) else NA_real_,
  DR_gene_order_equal = identical(DRs$gene, DRo$gene),
  DR_identical = identical(DRs[num], DRo[num]),
  DR_max_abs_diff = if (identical(DRs$gene, DRo$gene)) max(sapply(num, function(k) max(abs(DRs[[k]] - DRo[[k]]), na.rm = TRUE))) else NA_real_,
  n_sig_shared = sum(DRs$p.adj[DRs$gene != g] < 0.05), n_sig_official = sum(DRo$p.adj[DRo$gene != g] < 0.05))
fwrite(out, file.path(O, sprintf("equivalence_check_%s.tsv", g)), sep = "\t"); print(out); cat("DONE\n")
