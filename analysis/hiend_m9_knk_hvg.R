# Module M9 virtual knockout (scTenifoldKnk 1.0.3): HVG-equivalent implementation (replaces hiend_m9_knk_shared.R)
# Audited = template scTenifoldKnk workflow: disease group only (= Tumor) -> countMat[unique(c(target, hvgs)), ] (vst 5000 HVG)
#   -> scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500) -> positive if >= 1 gene other than the target has p.adj < 0.05
# Equivalence: only genes within the HVGs are knocked out -> the template input is the same 5000-gene matrix for every target -> with the same seed, makeNetworks / tensorDecomposition /
#   strictDirection give a WT network independent of the target and bit-identical -> building WT once is equivalent to calling the reference function per gene; the KO part follows scTenifoldKnk() line by line.
#   Verification: hiend_m9_knk_official_check.R runs the reference function unchanged for one gene; diffRegulation must be bit-identical.
# In the reference function body `countMatrix[rowSums(countMatrix != 0) >= qc_minCells, ]` is not assigned (has no effect); this filter is likewise not applied here.
# Seed: set.seed(123) (template), placed before scQC, matching the order "set.seed, then call the reference function".
# Environment: KNK_CORES (default 4; 24 cores exhaust memory); KNK_KOFILE (KO gene list tsv with a gene column; default = 2 timing genes); KNK_KODIR
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table)
  library(scTenifoldKnk); library(scTenifoldNet) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out", Sys.getenv("KNK_OUT", "hiend_m9_hvg")); dir.create(O, showWarnings = FALSE, recursive = TRUE)
CT <- Sys.getenv("KNK_CELLTYPE", "")   # KNK_CELLTYPE=Hepatocyte -> build the network from tumour hepatocytes only (default empty = all cells of tumour samples, as in the template)
BUILD_ONLY <- Sys.getenv("KNK_BUILD_ONLY", "0") == "1"
KNK_CORES <- as.integer(Sys.getenv("KNK_CORES", "4"))
KODIR <- file.path(O, Sys.getenv("KNK_KODIR", "ko")); dir.create(KODIR, showWarnings = FALSE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## ---- 1. Tumour-cell count matrix + HVG (as in the template) ----
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

## ---- 2. Wild-type network (built once; follows the reference function body line by line, default parameters) ----
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

## ---- 3. KO genes ----
kf <- Sys.getenv("KNK_KOFILE", "")
if (nzchar(kf)) { ko <- fread(kf) } else {
  old <- fread(file.path(B, "N11_out/hiend_m9/ko_genes.tsv"))   # timing only: first published gene and first random gene within the HVGs from the earlier list
  ko <- rbind(old[kind == "published" & gene %in% hvgs][1], old[kind == "random" & gene %in% hvgs][1])
}
miss <- setdiff(ko$gene, rownames(WT)); if (length(miss)) tm(paste("not in WT (dropped):", paste(miss, collapse = ",")))
ko <- ko[gene %in% rownames(WT)]
fwrite(ko, file.path(KODIR, "_ko_list_used.tsv"), sep = "\t")

## ---- 4. Per-gene knockout (resumable) ----
for (i in seq_len(nrow(ko))) {
  g <- ko$gene[i]; f <- file.path(KODIR, paste0(g, ".tsv")); if (file.exists(f)) next
  t1 <- Sys.time()
  KO <- WT; KO[g, ] <- 0
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2, nCores = KNK_CORES)   # inside the reference scTenifoldKnk() this function resolves to scTenifoldNet (checked via environment); dRegulation resolves to scTenifoldKnk's own version
  DR <- as.data.table(scTenifoldKnk:::dRegulation(MA, g)); rm(KO, MA)
  saveRDS(DR, sub("\\.tsv$", ".rds", f))   # full precision kept for bit-wise comparison with the reference (fwrite writes 15 digits)
  fwrite(DR, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  tm(sprintf("KO %d/%d %s: %.1f min, n_sig=%d", i, nrow(ko), g, as.numeric(difftime(Sys.time(), t1, units = "mins")), sum(DR[gene != g]$p.adj < 0.05)))
}
tm("DONE")
