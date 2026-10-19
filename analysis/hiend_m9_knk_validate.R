# Small-scale validation of the M9 shared-network approach on a reduced gene set, answering three questions
#  V1 code path: same input and seed, are reference scTenifoldKnk() and "shared WT + second half" bit-identical
#  V2 actual template practice: set.seed(12345) at the start of each run, input = HVG + that target -> how similar are WT networks across targets
#  V3 is the result "the same whatever is knocked out": pairwise Jaccard of significant-gene sets and FC Spearman for 6 KO genes; reference per-gene and shared versions
#  V4 seed noise: same gene, reference function, different seeds -> variability of the reference itself
# Reduction: network genes = top NG HVGs + 6 KO genes; other parameters as in the template (nNet=10, nCells=500, qc_minLSize=500)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(scTenifoldKnk); library(scTenifoldNet) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
O <- file.path(B, "N11_out/hiend_m9/validate"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
NG <- as.integer(Sys.getenv("KNK_NG", "1000"))
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))

map <- fread(file.path(B, "N11_out/hiend_m8/gene_symbol_map.tsv"))
M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"), data.table = FALSE)
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); stopifnot(identical(rownames(X), map$orig))
X <- X[!map$dropped_dup, M$Cell[M$site == "Tumor"]]; rownames(X) <- map$hgnc[!map$dropped_dup]
so <- CreateSeuratObject(X, min.cells = 5, min.features = 100); rm(X); gc()
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 5000); hv <- VariableFeatures(so)[seq_len(NG)]
cm <- GetAssayData(so, layer = "counts"); rm(so); gc()

ko <- fread(file.path(B, "N11_out/hiend_m9/ko_genes.tsv"))
set.seed(1); kg <- c(ko[kind == "published"][sample(.N, 3)]$gene, ko[kind == "random"][sample(.N, 3)]$gene)
kg <- setdiff(kg, hv)[1:6]; kg <- kg[!is.na(kg)]
fwrite(ko[gene %in% kg], file.path(O, "validate_ko_genes.tsv"), sep = "\t"); tm(paste("KO genes:", paste(kg, collapse = ",")))

official <- function(genes, g, seed) { set.seed(seed)
  r <- scTenifoldKnk(countMatrix = as.matrix(cm[genes, ]), gKO = g, qc_mtThreshold = 0.1, qc_minLSize = 500, nc_nNet = 10, nc_nCells = 500)
  list(WT = as.matrix(r$tensorNetworks$WT), DR = as.data.table(r$diffRegulation)) }
buildWT <- function(genes, seed) { cmx <- scTenifoldKnk:::scQC(as.matrix(cm[genes, ]), mtThreshold = 0.1, minLSize = 500); set.seed(seed)
  W <- makeNetworks(X = cmx, q = 0.9, nNet = 10, nCells = 500, scaleScores = TRUE, symmetric = FALSE, nComp = 3, nCores = parallel::detectCores())
  W <- tensorDecomposition(xList = W, K = 3, maxError = 1e-05, maxIter = 1000, nDecimal = 3)$X
  W <- as.matrix(scTenifoldKnk:::strictDirection(W, lambda = 0)); diag(W) <- 0; t(W) }
koDR <- function(WT, g) { KO <- WT; KO[g, ] <- 0
  as.data.table(scTenifoldKnk:::dRegulation(scTenifoldKnk:::manifoldAlignment(WT, KO, d = 2, nCores = parallel::detectCores()), g)) }
sig <- function(d, g) d[gene != g & p.adj < 0.05]$gene
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
fcr <- function(d1, d2, g1, g2) { m <- merge(d1[!gene %in% c(g1, g2), .(gene, a = FC)], d2[!gene %in% c(g1, g2), .(gene, b = FC)], by = "gene"); cor(m$a, m$b, method = "spearman") }

## V1: bit-wise identity (reference vs shared path, same input HVG + all 6 KO genes, same seed)
gall <- unique(c(hv, kg))
o1 <- official(gall, kg[1], 12345); tm("V1 official done")
W1 <- buildWT(gall, 12345); d1 <- koDR(W1, kg[1]); tm("V1 shared done")
v1 <- data.table(check = c("WT_max_abs_diff", "DR_FC_max_abs_diff", "DR_padj_max_abs_diff", "sig_identical"),
  value = c(max(abs(o1$WT[rownames(W1), colnames(W1)] - W1)),
            max(abs(merge(o1$DR, d1, by = "gene")[, FC.x - FC.y])), max(abs(merge(o1$DR, d1, by = "gene")[, p.adj.x - p.adj.y])),
            as.numeric(setequal(sig(o1$DR, kg[1]), sig(d1, kg[1])))))
fwrite(v1, file.path(O, "V1_code_path_identity.tsv"), sep = "\t"); print(v1)

## V3-shared: knock out 6 genes on the shared WT
DS <- setNames(lapply(kg, function(g) koDR(W1, g)), kg); tm("V3 shared KOs done")
## V2 + V3-official: actual template practice (each gene separately: HVG + that gene, set.seed(12345))
DO <- list(); WO <- list()
for (g in kg) { r <- official(unique(c(hv, g)), g, 12345); DO[[g]] <- r$DR; WO[[g]] <- r$WT; tm(paste("official", g)) }
## V4: same gene, reference function, different seeds
r4 <- official(unique(c(hv, kg[1])), kg[1], 2026); tm("V4 done")

pairs <- CJ(a = kg, b = kg)[a < b]
v3 <- rbind(
  pairs[, .(arm = "shared_WT", a, b, jaccard = mapply(function(x, y) jac(sig(DS[[x]], x), sig(DS[[y]], y)), a, b),
            fc_spearman = mapply(function(x, y) fcr(DS[[x]], DS[[y]], x, y), a, b))],
  pairs[, .(arm = "official_per_gene", a, b, jaccard = mapply(function(x, y) jac(sig(DO[[x]], x), sig(DO[[y]], y)), a, b),
            fc_spearman = mapply(function(x, y) fcr(DO[[x]], DO[[y]], x, y), a, b))])
fwrite(v3, file.path(O, "V3_between_KO_similarity.tsv"), sep = "\t")
# Same gene: shared vs reference per-gene; reference seed 12345 vs seed 2026
v34 <- rbind(
  data.table(compare = "shared_vs_official_same_gene", gene = kg,
             jaccard = sapply(kg, function(g) jac(sig(DS[[g]], g), sig(DO[[g]], g))),
             fc_spearman = sapply(kg, function(g) fcr(DS[[g]], DO[[g]], g, g)),
             pos_shared = sapply(kg, function(g) length(sig(DS[[g]], g)) > 0), pos_official = sapply(kg, function(g) length(sig(DO[[g]], g)) > 0),
             nsig_shared = sapply(kg, function(g) length(sig(DS[[g]], g))), nsig_official = sapply(kg, function(g) length(sig(DO[[g]], g)))),
  data.table(compare = "official_seed12345_vs_seed2026", gene = kg[1], jaccard = jac(sig(DO[[kg[1]]], kg[1]), sig(r4$DR, kg[1])),
             fc_spearman = fcr(DO[[kg[1]]], r4$DR, kg[1], kg[1]), pos_shared = NA, pos_official = length(sig(r4$DR, kg[1])) > 0,
             nsig_shared = NA, nsig_official = length(sig(r4$DR, kg[1]))))
fwrite(v34, file.path(O, "V34_same_gene_agreement.tsv"), sep = "\t")
# V2: correlation of WT networks between targets under the reference per-gene practice (on common genes)
cg <- Reduce(intersect, lapply(WO, rownames))
v2 <- pairs[, .(a, b, wt_cor = mapply(function(x, y) cor(as.vector(WO[[x]][cg, cg]), as.vector(WO[[y]][cg, cg])), a, b))]
fwrite(v2, file.path(O, "V2_official_WT_between_targets.tsv"), sep = "\t")
# Significant genes of each KO (is it always the same set)
fwrite(rbind(rbindlist(lapply(kg, function(g) data.table(arm = "shared_WT", ko = g, gene = sig(DS[[g]], g)))),
             rbindlist(lapply(kg, function(g) data.table(arm = "official_per_gene", ko = g, gene = sig(DO[[g]], g))))),
       file.path(O, "sig_genes_long.tsv"), sep = "\t")
print(v34); print(v3[, .(jac_mean = mean(jaccard, na.rm = TRUE), fc_mean = mean(fc_spearman)), by = arm]); print(summary(v2$wt_cor))
tm("DONE")
