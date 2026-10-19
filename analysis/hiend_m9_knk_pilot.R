# 高端模块 M9 虚拟敲除（scTenifoldKnk）—— 试跑（2026-10-15）
# 被审计对象 = 课程 265 `37.scTenifoldKnk` 写法（判据照抄被审计模板）：
#   只留疾病组（这里 = Tumor）→ vst 5000 HVG + 靶基因 → scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500)
#   → diffRegulation 去掉靶基因后 p.adj<0.05 的基因数 ≥1 即"敲除有下游效应"（阳性）
# KO 基因选取规则（模板是作者手挑，这里定成可复现规则）：签名内 |coef × TCGA sd| 最大、且在肿瘤细胞里检出的基因
# 试跑：1 个已发表签名（S001）的 KO 基因 + 2 个随机基因（U 池，固定种子，要求肿瘤细胞检出率 ≥1%），记录耗时
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(scTenifoldKnk) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m9"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE)
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds"))
stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, M$Cell[M$site == "Tumor"]]; rownames(X) <- map$hgnc[!map$dropped_dup]
so <- CreateSeuratObject(X, min.cells = 5, min.features = 100); rm(X); gc()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000); hvgs <- VariableFeatures(so)
cm <- GetAssayData(so, layer = "counts"); det <- Matrix::rowMeans(cm > 0)
tm(sprintf("tumor cells %d, genes %d", ncol(cm), nrow(cm)))

L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8")[sig_id == "S001"]
L <- unique(L[, .(gene_hgnc, coef = as.numeric(coef))])
# 试跑只用 |coef| 排序（S001 各基因 sd 量级相近）；正式阶段改 |coef × TCGA sd|
L <- L[gene_hgnc %in% rownames(cm)][order(-abs(coef))]
U <- readLines(file.path(B, "N11_out/hiend_m8/U_pool.txt")); U <- U[U %in% rownames(cm) & det[match(U, rownames(cm))] >= 0.01]
set.seed(20261015); ko <- data.table(kind = c("real_S001", "null_1", "null_2"), gene = c(L$gene_hgnc[1], sample(U, 2)))
ko[, detect_rate := det[gene]]; ko[, in_hvg := gene %in% hvgs]
fwrite(ko, file.path(O, "pilot_ko_genes.tsv"), sep = "\t"); print(ko)

res <- list()
for (i in seq_len(nrow(ko))) {
  g <- ko$gene[i]; f <- file.path(O, sprintf("pilot_diffreg_%s.tsv", g))
  if (file.exists(f)) next
  t1 <- Sys.time(); set.seed(12345)
  dat <- as.matrix(cm[unique(c(g, hvgs)), ])
  r <- scTenifoldKnk(countMatrix = dat, gKO = g, qc_mtThreshold = 0.1, qc_minLSize = 500, nc_nNet = 10, nc_nCells = 500)
  df <- as.data.table(r$diffRegulation); fwrite(df, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  df <- df[gene != g]
  res[[i]] <- data.table(ko[i], n_tested = nrow(df), n_sig = sum(df$p.adj < 0.05), pos = as.integer(sum(df$p.adj < 0.05) >= 1),
                         minutes = round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1))
  fwrite(rbindlist(res), file.path(O, "pilot_summary.tsv"), sep = "\t"); tm(sprintf("KO %s done", g)); print(res[[i]])
}
tm("DONE")
