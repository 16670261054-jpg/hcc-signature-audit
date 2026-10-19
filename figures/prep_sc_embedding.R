# Figure preparation: GSE149614 single-cell UMAP embedding + per-cell signature scores (display only, no new criteria or tests)
# QC identical to hiend_m8_sc_pilot.R (same cells as used for the M8 calls); cell types from the original authors' celltype annotation.
# Embedding: LogNormalize -> 2000 HVG -> PCA 30 -> harmony (by sample) -> UMAP.
# Scores: unweighted mean of signature genes in the data layer (colMeans computed directly for the large v5 object, not AddModuleScore).
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(harmony) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
M8 <- file.path(B, "N11_out/hiend_m8")
O <- file.path(B, "figures/data"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
set.seed(1234)

## Same symbol-update table as the per-gene M8 script (reads its saved map, so results are bit-identical)
map <- fread(file.path(M8, "gene_symbol_map.tsv"))
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds"))
stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, ]; rownames(X) <- map$hgnc[!map$dropped_dup]
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE); rownames(M) <- M$Cell
so <- CreateSeuratObject(X, meta.data = M[colnames(X), ], min.cells = 5, min.features = 100); rm(X); gc()
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 100 & percent.mt < 10)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
## Cell-count check: must match cells_site_celltype.tsv written by hiend_m8_sc_pilot.R
chk <- fread(file.path(M8, "cells_site_celltype.tsv"))
now <- as.data.table(table(site = so$site, celltype = so$celltype))
mm <- merge(chk, now, by = c("site", "celltype"))
stopifnot(nrow(mm) == nrow(chk), all(mm$N.x == mm$N.y))
tm(sprintf("cells match pilot: %d", ncol(so)))

so <- FindVariableFeatures(so, nfeatures = 2000)
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunHarmony(so, group.by.vars = "sample", verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:30, seed.use = 1234, verbose = FALSE)
tm("umap done")
E <- as.data.table(Embeddings(so, "umap"), keep.rownames = "cell"); setnames(E, c("cell", "UMAP1", "UMAP2"))
md <- so@meta.data
E[, `:=`(celltype = md$celltype, site = md$site, patient = md$patient, sample = md$sample, nFeature = md$nFeature_RNA)]
fwrite(E, file.path(O, "sc_umap.tsv.gz"), sep = "\t")

## Per-cell signature scores: 148 published (A real) + first replicate of size-matched random Ak per source signature + 8 real B sets
G <- fread(file.path(M8, "signature_genes_all_arms.tsv"))
G <- G[(arm == "A" & kind == "real") | (arm == "Ak" & rep == 1) | (arm == "B" & kind == "real")]
G[, sid := paste(arm, src, kind, rep, sep = "|")]
dat <- GetAssayData(so, assay = "RNA", layer = "data")
sids <- unique(G$sid)
S <- vapply(sids, function(s) { g <- intersect(G[sid == s]$gene, rownames(dat))
  if (!length(g)) return(rep(NA_real_, ncol(dat))); Matrix::colMeans(dat[g, , drop = FALSE]) }, numeric(ncol(dat)))
rownames(S) <- colnames(dat)
saveRDS(S, file.path(O, "sc_sigscore.rds"))
cov <- G[, .(k = .N, k_in_sc = sum(gene %in% rownames(dat))), by = .(sid, arm, src, kind, rep)]
fwrite(cov, file.path(O, "sc_sigscore_coverage.tsv"), sep = "\t")

## Mean expression and expression fraction of signature genes by cell type (for the dot plot)
ug <- intersect(unique(G$gene), rownames(dat))
ct <- so$celltype
AV <- rbindlist(lapply(sort(unique(ct)), function(k) { m <- dat[ug, ct == k, drop = FALSE]
  data.table(celltype = k, gene = ug, avg = Matrix::rowMeans(expm1(m)), pct = Matrix::rowMeans(m > 0)) }))
fwrite(AV, file.path(O, "sc_gene_celltype_avg.tsv"), sep = "\t")
tm("all done")
