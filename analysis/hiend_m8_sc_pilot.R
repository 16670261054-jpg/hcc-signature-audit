# 高端模块 M8 单细胞定位 —— 试跑（2026-10-15）
# 被审计对象 = 课程 265 `12.Seurat` 的写法（判据照抄被审计模板，不是我们的判据）：
#   M8a "细胞类型定位"：FindAllMarkers(only.pos=FALSE, min.pct=0.25, logfc.threshold=1)，
#        |avg_log2FC|>1 & p_val_adj<0.05 → 签名中 ≥1 个基因进 sig.cellMarkers 即阳性
#   M8b "单细胞验证差异"：每个细胞类型 Tumor vs Normal FindMarkers(group.by='group', logfc.threshold=0.1)，
#        两组各 >5 个细胞才比；|avg_log2FC|>1 & p_val_adj<0.05 → ≥1 个签名基因、≥1 个细胞类型即阳性
# 模板的检验对每个基因与所属基因集无关 → 全基因算一次落盘，随机集/真实集查表即可（与逐集重跑等价）
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(survival) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m8"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## 基因符号更新：与 load_cohorts.R 同一张 HGNC 表、同一规则
hg <- fread(file.path(B, "N12_signatures/ref/hgnc.txt"), select = c("symbol", "entrez_id", "status", "prev_symbol"), data.table = FALSE, quote = "")
hg <- hg[hg$status == "Approved", ]
sp <- strsplit(hg$prev_symbol, "|", fixed = TRUE); pv <- data.frame(old = trimws(unlist(sp)), new = rep(hg$symbol, lengths(sp))); pv <- pv[nzchar(pv$old), ]
pv <- pv[!duplicated(pv$old) & !(pv$old %in% pv$old[duplicated(pv$old)]) & !(pv$old %in% hg$symbol), ]
upd <- function(x) { m <- match(x, pv$old); ifelse(is.na(m), x, pv$new[m]) }

## 抽样池 U：与第 3 环同一定义（load_cohorts.R），首次运行时落盘
if (!file.exists(file.path(O, "U_pool.txt"))) {
  TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA"); source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
  writeLines(U, file.path(O, "U_pool.txt")); writeLines(Uall, file.path(O, "Uall_score.txt"))
  rm(list = setdiff(ls(), c("B", "D", "O", "t0", "tm", "upd", "pv", "hg"))); gc()
}

X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); tm(sprintf("counts %d x %d", nrow(X), ncol(X)))
new <- upd(rownames(X)); dup <- duplicated(new)
map <- data.table(orig = rownames(X), hgnc = new, changed = rownames(X) != new, dropped_dup = dup)
fwrite(map, file.path(O, "gene_symbol_map.tsv"), sep = "\t")
X <- X[!dup, ]; rownames(X) <- new[!dup]; tm(sprintf("symbols updated: %d changed, %d duplicates dropped", sum(map$changed), sum(dup)))

M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE); rownames(M) <- M$Cell
stopifnot(all(colnames(X) %in% rownames(M)))
so <- CreateSeuratObject(X, meta.data = M[colnames(X), ], min.cells = 5, min.features = 100); rm(X); gc()
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 100 & percent.mt < 10)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
tm(sprintf("seurat object: %d genes x %d cells", nrow(so), ncol(so)))
fwrite(as.data.table(table(site = so$site, celltype = so$celltype)), file.path(O, "cells_site_celltype.tsv"), sep = "\t")

## M8a
Idents(so) <- "celltype"
mk <- FindAllMarkers(so, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 1)
mk <- as.data.table(mk); fwrite(mk, file.path(O, "m8a_findallmarkers_all.tsv"), sep = "\t")
tm(sprintf("M8a FindAllMarkers rows %d; sig %d", nrow(mk), mk[abs(avg_log2FC) > 1 & p_val_adj < 0.05, .N]))

## M8b
so$group <- paste0(so$site, "_", so$celltype)
res <- list()
for (ct in sort(unique(so$celltype))) {
  a <- paste0("Tumor_", ct); b <- paste0("Normal_", ct)
  if (sum(so$group == a) > 5 && sum(so$group == b) > 5) {
    r <- FindMarkers(so, ident.1 = a, ident.2 = b, group.by = "group", logfc.threshold = 0.1)
    res[[ct]] <- data.table(celltype = ct, gene = rownames(r), r)
    tm(sprintf("M8b %s: tested %d, sig %d", ct, nrow(r), sum(abs(r$avg_log2FC) > 1 & r$p_val_adj < 0.05)))
  }
}
gd <- rbindlist(res); fwrite(gd, file.path(O, "m8b_tumor_vs_normal_all.tsv"), sep = "\t")

## 试跑判读：每个基因是否"阳性"；按尺寸 k 的解析期望 1-(1-p)^k 只作试跑参考，正式阶段用与第 3 环同种子的真实随机集查表
Ugenes <- readLines(file.path(O, "U_pool.txt"))
pos_a <- unique(mk[abs(avg_log2FC) > 1 & p_val_adj < 0.05]$gene); pos_b <- unique(gd[abs(avg_log2FC) > 1 & p_val_adj < 0.05]$gene)
pa <- mean(Ugenes %in% pos_a); pb <- mean(Ugenes %in% pos_b)
out <- data.table(k = c(3, 5, 8, 10, 15, 20), M8a_expected = 1 - (1 - pa)^c(3, 5, 8, 10, 15, 20), M8b_expected = 1 - (1 - pb)^c(3, 5, 8, 10, 15, 20))
out[, `:=`(p_gene_M8a = pa, p_gene_M8b = pb, U_n = length(Ugenes), U_in_sc = sum(Ugenes %in% rownames(so)))]
fwrite(out, file.path(O, "pilot_expected_fpr.tsv"), sep = "\t")
tm("DONE")
