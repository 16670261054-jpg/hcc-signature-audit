# 高端模块 M9 虚拟敲除（scTenifoldKnk 1.0.3）—— HVG 等价版（2026-10-16，替代 hiend_m9_knk_shared.R）
# 被审计对象 = 课程 265 `37.scTenifoldKnk`：只留疾病组（= Tumor）→ countMat[unique(c(target, hvgs)), ]（vst 5000 HVG）
#   → scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500) → 去靶基因后 p.adj<0.05 的基因数 ≥1 = 阳性
# 等价性论证：只敲 HVG 内的基因 → 模板对每个靶基因的输入都是同一个 5000 基因矩阵 → 同种子下 makeNetworks / tensorDecomposition /
#   strictDirection 产出的 WT 与靶基因无关、逐位相同 → WT 建一次即与官方逐个调用等价；KO 段逐行照抄 scTenifoldKnk() 函数体。
#   验证：hiend_m9_knk_official_check.R 用官方函数原样跑 1 个基因，diffRegulation 必须逐位一致。
# 官方函数体里 `countMatrix[rowSums(countMatrix != 0) >= qc_minCells, ]` 未赋值（不生效），此处同样不做该过滤以保持一致。
# 种子：set.seed(123)（课程 246 26.scTenifoldKnk 写法），放在 scQC 之前，与"set.seed 后调用官方函数"的顺序一致。
# 环境变量：KNK_CORES（默认 4；24 核会爆内存）；KNK_KOFILE（KO 基因清单 tsv，需含 gene 列；缺省=测速用 2 个）；KNK_KODIR
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table)
  library(scTenifoldKnk); library(scTenifoldNet) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out", Sys.getenv("KNK_OUT", "hiend_m9_hvg")); dir.create(O, showWarnings = FALSE, recursive = TRUE)
CT <- Sys.getenv("KNK_CELLTYPE", "")   # 10-17 方案 A：KNK_CELLTYPE=Hepatocyte → 只用肿瘤肝细胞建网（默认空 = 肿瘤样本全部细胞，模板原样）
BUILD_ONLY <- Sys.getenv("KNK_BUILD_ONLY", "0") == "1"
KNK_CORES <- as.integer(Sys.getenv("KNK_CORES", "4"))
KODIR <- file.path(O, Sys.getenv("KNK_KODIR", "ko")); dir.create(KODIR, showWarnings = FALSE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## ---- 1. 肿瘤细胞计数矩阵 + HVG（照抄模板） ----
fi <- file.path(O, "input_hvg5000.rds")
if (!file.exists(fi)) {
  map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
  M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE)
  X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); stopifnot(identical(rownames(X), map$orig))
  cells <- M$Cell[M$site == "Tumor" & (CT == "" | M$celltype == CT)]
  X <- X[!map$dropped_dup, cells]; rownames(X) <- map$hgnc[!map$dropped_dup]
  so <- CreateSeuratObject(X, min.cells = 5, min.features = 100); rm(X); gc()
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000); hvgs <- VariableFeatures(so)
  cm <- GetAssayData(so, layer = "counts"); rm(so); gc()
  inp <- as.matrix(cm[hvgs, ]); rm(cm); gc()
  saveRDS(inp, paste0(fi, ".tmp")); file.rename(paste0(fi, ".tmp"), fi)
  writeLines(hvgs, file.path(O, "hvgs.txt"))
}
inp <- readRDS(fi); hvgs <- readLines(file.path(O, "hvgs.txt"))
tm(sprintf("input %d genes x %d tumor cells", nrow(inp), ncol(inp)))

## ---- 2. 野生型网络（建一次；逐行照抄官方函数体，默认参数照抄） ----
fw <- file.path(O, "WT_network.rds")
if (!file.exists(fw)) {
  set.seed(123)
  countMatrix <- scTenifoldKnk:::scQC(inp, mtThreshold = 0.1, minLSize = 500)
  writeLines(rownames(countMatrix), file.path(O, "genes_after_qc.txt"))
  tm(sprintf("QC: %d genes x %d cells", nrow(countMatrix), ncol(countMatrix)))
  WT <- makeNetworks(X = countMatrix, q = 0.9, nNet = 10, nCells = 500, scaleScores = TRUE, symmetric = FALSE, nComp = 3, nCores = KNK_CORES)
  tm("makeNetworks done")
  WT <- tensorDecomposition(xList = WT, K = 3, maxError = 1e-05, maxIter = 1000, nDecimal = 3)
  WT <- WT$X; WT <- scTenifoldKnk:::strictDirection(WT, lambda = 0); WT <- as.matrix(WT); diag(WT) <- 0; WT <- t(WT)
  saveRDS(WT, paste0(fw, ".tmp")); file.rename(paste0(fw, ".tmp"), fw)
  tm("tensorDecomposition done, WT saved")
}
WT <- readRDS(fw)
if (BUILD_ONLY) { tm("BUILD_ONLY: network ready, stop"); quit(save = "no", status = 0) }

## ---- 3. KO 基因 ----
kf <- Sys.getenv("KNK_KOFILE", "")
if (nzchar(kf)) { ko <- fread(kf) } else {
  old <- fread(file.path(B, "N11_out/hiend_m9/ko_genes.tsv"))   # 仅测速：取旧清单中第一个 HVG 内的已发表基因与随机基因
  ko <- rbind(old[kind == "published" & gene %in% hvgs][1], old[kind == "random" & gene %in% hvgs][1])
}
miss <- setdiff(ko$gene, rownames(WT)); if (length(miss)) tm(paste("not in WT (dropped):", paste(miss, collapse = ",")))
ko <- ko[gene %in% rownames(WT)]
fwrite(ko, file.path(KODIR, "_ko_list_used.tsv"), sep = "\t")

## ---- 4. 逐基因敲除（可续跑） ----
for (i in seq_len(nrow(ko))) {
  g <- ko$gene[i]; f <- file.path(KODIR, paste0(g, ".tsv")); if (file.exists(f)) next
  t1 <- Sys.time()
  KO <- WT; KO[g, ] <- 0
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2, nCores = KNK_CORES)   # 官方 scTenifoldKnk() 内此函数解析到 scTenifoldNet（10-17 用 environment 核实）；dRegulation 解析到 scTenifoldKnk 自己的版本
  DR <- as.data.table(scTenifoldKnk:::dRegulation(MA, g)); rm(KO, MA)
  saveRDS(DR, sub("\\.tsv$", ".rds", f))   # 全精度留档，供与官方逐位比对（fwrite 只写 15 位）
  fwrite(DR, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  tm(sprintf("KO %d/%d %s: %.1f min, n_sig=%d", i, nrow(ko), g, as.numeric(difftime(Sys.time(), t1, units = "mins")), sum(DR[gene != g]$p.adj < 0.05)))
}
tm("DONE")
