# Module M10 Geneformer: input export
# Audited configuration = template Geneformer perturbation pipeline: Pretrained V1-10M, delete all genes, start=Disease -> goal=Control,
#   max_ncells=1000, emb_layer=-1, InSilicoPerturberStats(mode="goal_state_shift") -> Goal_end_FDR
# Mapping: Disease = tumour hepatocytes (Tumor_Hepatocyte), Control = adjacent normal hepatocytes (Normal_Hepatocyte)
# Whole-transcriptome counts (no HVG truncation); donor-balanced sampling; n_counts = whole-transcriptome library size
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m10/gf"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
CAP <- 150L   # max cells per group per donor

map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"))
M <- M[celltype == "Hepatocyte" & site %in% c("Tumor", "Normal")]
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, M$Cell]; rownames(X) <- map$hgnc[!map$dropped_dup]
M[, nCount := Matrix::colSums(X)]; M[, nFeature := Matrix::colSums(X > 0)]
M[, mt := Matrix::colSums(X[grepl("^MT-", rownames(X)), ]) / nCount]
M <- M[nFeature > 100 & mt < 0.1]          # same QC as M8

set.seed(20261015)
S <- M[, .SD[sample(.N, min(.N, CAP))], by = .(site, patient)]
S[, Type := ifelse(site == "Tumor", "Disease", "Control")]
aud <- S[, .N, by = .(Type, patient)][order(Type, -N)]
aud[, share := N / sum(N), by = Type]
fwrite(aud, file.path(O, "export_donor_audit.tsv"), sep = "\t")
X <- X[, S$Cell]; X <- X[Matrix::rowSums(X) > 0, ]
meta <- data.frame(row.names = S$Cell, group = S$Type, Type = S$Type, donor = S$patient, nCount = S$nCount)
writeMM(X, file.path(O, "counts.mtx")); R.utils_ok <- requireNamespace("R.utils", quietly = TRUE)
if (R.utils_ok) R.utils::gzip(file.path(O, "counts.mtx"), overwrite = TRUE) else system2("gzip", c("-f", shQuote(file.path(O, "counts.mtx"))))
fwrite(data.table(gene_symbol = rownames(X)), file.path(O, "features.tsv"), sep = "\t")
write.table(meta, file.path(O, "cell_meta.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cat("cells:", ncol(X), " genes:", nrow(X), "\n"); print(aud[, .(n = sum(N), donors = .N, top_share = max(share)), by = Type])
cat("DONE\n")
