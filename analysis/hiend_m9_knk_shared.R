# 高端模块 M9 虚拟敲除（scTenifoldKnk）—— 共享野生型网络版（2026-10-15）
# 被审计对象 = 课程 265 `37.scTenifoldKnk` 写法：只留疾病组（= Tumor）→ vst 5000 HVG + 靶基因
#   → scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500) → 去靶基因后 p.adj<0.05 的基因数 ≥1 = 阳性
# 为什么能共享：scTenifoldKnk() 函数体里 makeNetworks → tensorDecomposition → strictDirection 与 gKO 无关，
#   gKO 只出现在 `KO[gKO, ] <- 0` 及之后的 manifoldAlignment / dRegulation。下面逐行照抄函数体（scTenifoldKnk 1.0.3），
#   野生型网络建一次落盘，每个 KO 基因只做后半段。
# 与模板唯一的差别：模板是 HVG 5000 + 单个靶基因；这里是 HVG 5000 + 全部 KO 基因（已发表签名 KO 基因 + 随机基因）一起建网。
# KO 基因规则（模板是作者手挑，这里定成可复现规则）：签名内 |coef × TCGA sd| 最大、且过 QC 的基因
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(survival)
  library(scTenifoldKnk); library(scTenifoldNet) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m9"); dir.create(file.path(O, "ko"), showWarnings = FALSE, recursive = TRUE)
NRAND <- as.integer(Sys.getenv("KNK_NRAND", "200"))
KNK_CORES <- as.integer(Sys.getenv("KNK_CORES", "4"))   # 10-16：默认 detectCores()=24 时每个 worker 复制约 1 GiB 全局变量 → 32 GB 内存爆掉；核数不影响抽样（抽样在主进程串行）
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## ---- 1. KO 基因清单 ----
fk <- file.path(O, "ko_genes.tsv")
if (!file.exists(fk)) {
  TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
  source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
  sdv <- apply(cohorts$TCGA$x[Uall, ], 1, sd)
  L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
  I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
  L <- L[sig_id %in% I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id & !is.na(coef) & gene_hgnc %in% Uall]
  L <- unique(L, by = c("sig_id", "gene_hgnc")); L[, w := abs(coef * sdv[gene_hgnc])]
  pubko <- L[order(sig_id, -w)][, .(gene_rank = seq_len(.N), gene = gene_hgnc, w), by = sig_id]
  fwrite(pubko, file.path(O, "published_gene_weights.tsv"), sep = "\t")
  writeLines(U, file.path(O, "U_pool.txt"))
  rm(cohorts); gc()
}

## ---- 2. 肿瘤细胞计数矩阵 + HVG（照抄模板） ----
map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE)
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, M$Cell[M$site == "Tumor"]]; rownames(X) <- map$hgnc[!map$dropped_dup]
so <- CreateSeuratObject(X, min.cells = 5, min.features = 100); rm(X); gc()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000); hvgs <- VariableFeatures(so)
cm <- GetAssayData(so, layer = "counts"); rm(so); gc()
tm(sprintf("tumor cells %d, genes %d", ncol(cm), nrow(cm)))

## ---- 3. 野生型网络（建一次） ----
fw <- file.path(O, "WT_network.rds")
if (!file.exists(fk) || !file.exists(fw)) {
  pubko <- fread(file.path(O, "published_gene_weights.tsv")); U <- readLines(file.path(O, "U_pool.txt"))
  # 候选 KO 基因必须过 scQC 后仍在矩阵里：先按模板 QC 过滤细胞，再要求 ≥25 个细胞检出（qc_minCells 默认值）
  cq <- scTenifoldKnk:::scQC(cm[unique(c(hvgs, intersect(c(pubko$gene, U), rownames(cm)))), ], mtThreshold = 0.1, minLSize = 500)
  ok <- rownames(cq)[Matrix::rowSums(cq != 0) >= 25]
  pk <- pubko[gene %in% ok][, .SD[1], by = sig_id]
  set.seed(20261015); rk <- sample(setdiff(intersect(U, ok), pk$gene), NRAND)
  ko <- rbind(data.table(kind = "published", sig_id = pk$sig_id, gene = pk$gene),
              data.table(kind = "random", sig_id = sprintf("R%03d", seq_along(rk)), gene = rk))
  ko[, detect_rate := Matrix::rowMeans(cq[gene, ] > 0)]; ko[, in_hvg := gene %in% hvgs]
  fwrite(ko, fk, sep = "\t")
  netg <- unique(c(hvgs, ko$gene))
  tm(sprintf("network genes %d (HVG %d + KO extra %d); KO genes: published %d, random %d",
             length(netg), length(hvgs), length(netg) - length(hvgs), sum(ko$kind == "published"), sum(ko$kind == "random")))
  # 以下逐行照抄 scTenifoldKnk() 函数体（默认参数照抄）
  countMatrix <- as.matrix(cm[netg, ]); rm(cq); gc()
  countMatrix <- scTenifoldKnk:::scQC(countMatrix, mtThreshold = 0.1, minLSize = 500)
  tm(sprintf("QC: %d genes x %d cells", nrow(countMatrix), ncol(countMatrix)))
  set.seed(12345)
  WT <- makeNetworks(X = countMatrix, q = 0.9, nNet = 10, nCells = 500, scaleScores = TRUE, symmetric = FALSE, nComp = 3,
                     nCores = KNK_CORES)
  tm("makeNetworks done")
  WT <- tensorDecomposition(xList = WT, K = 3, maxError = 1e-05, maxIter = 1000, nDecimal = 3)
  WT <- WT$X; WT <- scTenifoldKnk:::strictDirection(WT, lambda = 0); WT <- as.matrix(WT); diag(WT) <- 0; WT <- t(WT)
  saveRDS(WT, paste0(fw, ".tmp")); file.rename(paste0(fw, ".tmp"), fw)
  tm("tensorDecomposition done, WT saved")
}
WT <- readRDS(fw); ko <- fread(fk)
ko <- ko[gene %in% rownames(WT)]

## ---- 4. 逐基因敲除（可续跑） ----
# 测速（2026-10-16）：KNK_MAXKO=N 只敲前 N 个；KNK_KODIR 指定输出子目录（测速用 ko_speedtest，与正式结果分开）
MAXKO <- as.integer(Sys.getenv("KNK_MAXKO", "0")); KODIR <- file.path(O, Sys.getenv("KNK_KODIR", "ko")); dir.create(KODIR, showWarnings = FALSE)
if (MAXKO > 0) ko <- rbind(ko[kind == "published"][1:ceiling(MAXKO / 2)], ko[kind == "random"][1:floor(MAXKO / 2)])
for (i in seq_len(nrow(ko))) {
  g <- ko$gene[i]; f <- file.path(KODIR, paste0(g, ".tsv")); if (file.exists(f)) next
  t1 <- Sys.time()
  KO <- WT; KO[g, ] <- 0
  MA <- scTenifoldKnk:::manifoldAlignment(WT, KO, d = 2, nCores = KNK_CORES)
  DR <- as.data.table(scTenifoldKnk:::dRegulation(MA, g)); rm(KO, MA)
  fwrite(DR, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  tm(sprintf("KO %d/%d %s (%s): %.1f min, n_sig=%d", i, nrow(ko), g, ko$kind[i],
             as.numeric(difftime(Sys.time(), t1, units = "mins")), sum(DR[gene != g]$p.adj < 0.05)))
}

## ---- 5. 汇总 ----
S <- rbindlist(lapply(seq_len(nrow(ko)), function(i) { g <- ko$gene[i]; f <- file.path(KODIR, paste0(g, ".tsv"))
  if (!file.exists(f)) return(NULL); d <- fread(f)[gene != g]
  data.table(ko[i], n_tested = nrow(d), n_sig = sum(d$p.adj < 0.05), pos = as.integer(sum(d$p.adj < 0.05) >= 1)) }))
fwrite(S, file.path(O, "m9_summary.tsv"), sep = "\t")
tm("DONE")
