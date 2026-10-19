# Module M9 virtual knockout (scTenifoldKnk): shared wild-type network version
# Audited = template scTenifoldKnk workflow: disease group only (= Tumor) -> vst 5000 HVG + target
#   -> scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500) -> positive if >= 1 gene other than the target has p.adj < 0.05
# Why sharing works: in scTenifoldKnk(), makeNetworks -> tensorDecomposition -> strictDirection do not depend on gKO;
#   gKO appears only in `KO[gKO, ] <- 0` and the subsequent manifoldAlignment / dRegulation. The function body (scTenifoldKnk 1.0.3) is reproduced line by line below;
#   the wild-type network is built once and saved, and each KO gene runs only the second half.
# Only difference from the template: template = 5000 HVG + a single target; here = 5000 HVG + all KO genes (published-signature KO genes + random genes) in one network.
# KO gene rule (chosen by hand in the template; made reproducible here): the gene with the largest |coef x TCGA sd| in the signature that passes QC
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(survival)
  library(scTenifoldKnk); library(scTenifoldNet) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m9"); dir.create(file.path(O, "ko"), showWarnings = FALSE, recursive = TRUE)
NRAND <- as.integer(Sys.getenv("KNK_NRAND", "200"))
KNK_CORES <- as.integer(Sys.getenv("KNK_CORES", "4"))   # with the default detectCores()=24 each worker copies ~1 GiB of globals -> 32 GB memory exhausted; core number does not affect sampling (sampling is serial in the main process)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

## ---- 1. KO gene list ----
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

## ---- 2. Tumour-cell count matrix + HVG (as in the template) ----
map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE)
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, M$Cell[M$site == "Tumor"]]; rownames(X) <- map$hgnc[!map$dropped_dup]
so <- CreateSeuratObject(X, min.cells = 5, min.features = 100); rm(X); gc()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000); hvgs <- VariableFeatures(so)
cm <- GetAssayData(so, layer = "counts"); rm(so); gc()
tm(sprintf("tumor cells %d, genes %d", ncol(cm), nrow(cm)))

## ---- 3. Wild-type network (built once) ----
fw <- file.path(O, "WT_network.rds")
if (!file.exists(fk) || !file.exists(fw)) {
  pubko <- fread(file.path(O, "published_gene_weights.tsv")); U <- readLines(file.path(O, "U_pool.txt"))
  # Candidate KO genes must remain after scQC: filter cells with the template QC first, then require detection in >= 25 cells (qc_minCells default)
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
  # The following reproduces the scTenifoldKnk() function body line by line (default parameters)
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

## ---- 4. Per-gene knockout (resumable) ----
# Timing: KNK_MAXKO=N knocks out only the first N genes; KNK_KODIR sets the output subdirectory (ko_speedtest for timing, separate from formal results)
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

## ---- 5. Summary ----
S <- rbindlist(lapply(seq_len(nrow(ko)), function(i) { g <- ko$gene[i]; f <- file.path(KODIR, paste0(g, ".tsv"))
  if (!file.exists(f)) return(NULL); d <- fread(f)[gene != g]
  data.table(ko[i], n_tested = nrow(d), n_sig = sum(d$p.adj < 0.05), pos = as.integer(sum(d$p.adj < 0.05) >= 1)) }))
fwrite(S, file.path(O, "m9_summary.tsv"), sep = "\t")
tm("DONE")
