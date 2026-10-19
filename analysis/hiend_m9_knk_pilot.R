# Module M9 virtual knockout (scTenifoldKnk): pilot
# Audited = template scTenifoldKnk workflow (criteria reproduced from the audited template):
#   disease group only (here = Tumor) -> vst 5000 HVG + target -> scTenifoldKnk(qc_mtThreshold=0.1, qc_minLSize=500, nc_nNet=10, nc_nCells=500)
#   -> positive ("knockout has downstream effects") if >= 1 gene other than the target has p.adj < 0.05 in diffRegulation
# KO gene rule (chosen by hand in the template; made reproducible here): the gene with the largest |coef x TCGA sd| in the signature that is detected in tumour cells
# Pilot: KO gene of one published signature (S001) + 2 random genes (pool U, fixed seed, detected in >= 1% of tumour cells); timing recorded
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
# Pilot ranks by |coef| only (genes of S001 have similar sd); the formal stage uses |coef x TCGA sd|
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
