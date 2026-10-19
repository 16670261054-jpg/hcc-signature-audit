# Modules M3 drug sensitivity + M12 DSigDB drug enrichment
# Audited = template workflow; criteria reproduced from the template (not our criteria):
#   M3 (template oncoPredict + boxplot steps):
#       GDSC2_Res=exp(GDSC2_Res); test matrix rowMeans>0.5; calcPhenotype(batchCorrect='eb', powerTransformPhenotype=TRUE,
#       removeLowVaryingGenes=0.2, minNumSamples=10, removeLowVaringGenesFrom='rawData'), set.seed(12345); per-drug Wilcoxon between high/low risk groups.
#       Positive: >= 1 drug with BH-FDR < 0.05 (the template's uncorrected p < 0.001 is reported separately).
#       Fixes: log2(x+1) if the test matrix is not on log scale; keep DRUG_ID column names without stripping suffixes; values are predicted IC50 (higher = more resistant).
#   M12 (template drug-enrichment step): enricher(genes, TERM2GENE=DSigDB[,1:2], minGSSize=10, maxGSSize=500,
#       pvalueCutoff=1, qvalueCutoff=1), positive = >= 1 drug with pvalue < 0.05 & p.adjust < 0.05 (template criterion).
# Signatures = signature_genes_all_arms.tsv written by hiend_m8_formal.R (same set as proxy analysis)
# Calls: random positive rate Clopper-Pearson 95% CI lower limit > 0.075 -> excess type I error; real vs random two-sided Fisher
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/hiend_m3_m12"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
DL <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "external")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))
G[, sid := paste(arm, src, kind, rep, sep = "|")]
sigs <- split(G, G$sid)
meta <- unique(G[, .(sid, arm, src, kind, rep)])

## ================= M3 =================
pred_f <- file.path(O, "DrugPredictions_TCGA.tsv")
if (!file.exists(pred_f)) {
  ex_dir <- file.path(DL, "oncoPredict_osf/extracted")
  if (!file.exists(file.path(ex_dir, "GDSC2_Res.rds"))) {
    dir.create(ex_dir, showWarnings = FALSE)
    for (n in c("GDSC2_Expr (RMA Normalized and Log Transformed).rds", "GDSC2_Res.rds")) {
      unzip(file.path(DL, "oncoPredict_osf/DataFiles.zip"), files = paste0("DataFiles/DataFiles/Training Data/", n), exdir = ex_dir, junkpaths = TRUE)
      stopifnot(file.exists(file.path(ex_dir, n))) }
  }
  trainExpr <- readRDS(file.path(ex_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
  trainRes <- readRDS(file.path(ex_dir, "GDSC2_Res.rds")); trainRes <- exp(trainRes)   # template: GDSC2_Res is ln(IC50); exp() before powerTransform
  X <- cohorts$TCGA$x[, colnames(cohorts$TCGA$z)]
  scale_note <- sprintf("TCGA max=%.2f", max(X, na.rm = TRUE))
  if (max(X, na.rm = TRUE) > 50) X <- log2(X + 1)                                       # ensure log scale
  X <- X[rowMeans(X) > 0.5, ]                                                            # template: data[rowMeans(data)>0.5,]
  fwrite(data.table(note = scale_note, logged = max(cohorts$TCGA$x, na.rm = TRUE) > 50, train_max = max(trainExpr)), file.path(O, "m3_scale_check.tsv"), sep = "\t")
  # oncoPredict 1.3.1: default folder=FALSE does not write DrugPredictions.csv (older versions did); results are only in the return value
  set.seed(12345)
  pm <- oncoPredict::calcPhenotype(trainingExprData = trainExpr, trainingPtype = trainRes, testExprData = X,
                             batchCorrect = "eb", powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2,
                             minNumSamples = 10, printOutput = FALSE, removeLowVaringGenesFrom = "rawData")
  if (is.list(pm)) pm <- pm$DrugPredictions
  stopifnot(is.matrix(pm), nrow(pm) == ncol(X), !anyDuplicated(colnames(pm)))           # keep DRUG_ID names
  pr <- data.table(sample = rownames(pm), pm)
  fwrite(pr, pred_f, sep = "\t"); tm(sprintf("calcPhenotype done: %d samples x %d drugs", nrow(pr), ncol(pr) - 1))
}
pr <- fread(pred_f); smp <- pr[[1]]; DM <- as.matrix(pr[, -1]); rownames(DM) <- smp
Z <- cohorts$TCGA$z; stopifnot(all(colnames(Z) %in% smp)); DM <- DM[colnames(Z), ]
m3_one <- function(d) { b <- setNames(d$coef, d$gene); s <- drop(t(Z[names(b), , drop = FALSE]) %*% b)
  hi <- s > median(s); p <- apply(DM, 2, function(v) wilcox.test(v[hi], v[!hi], exact = FALSE)$p.value)
  q <- p.adjust(p, "BH")
  data.table(sid = d$sid[1], M3_n_fdr05 = sum(q < 0.05), M3_pos = as.integer(any(q < 0.05)), M3_n_p001_template = sum(p < 0.001)) }
cl <- makeCluster(min(20, detectCores() - 2)); clusterExport(cl, c("Z", "DM", "m3_one")); invisible(clusterEvalQ(cl, library(data.table)))
R3 <- rbindlist(parLapply(cl, sigs, m3_one)); stopCluster(cl)
tm("M3 per-signature done")

## ================= M12 DSigDB =================
suppressPackageStartupMessages(library(clusterProfiler))
dr <- fread(file.path(DL, "DSigDB/DSigDB_All_detailed.txt"), sep = "\t", quote = "", select = 1:2, col.names = c("Drug", "Gene"))
dr <- unique(dr[!is.na(Gene) & Gene != ""])
m12_one <- function(d) { kk <- suppressMessages(enricher(unique(d$gene), pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, TERM2GENE = dr))
  r <- if (is.null(kk)) data.frame() else as.data.frame(kk)
  n <- if (nrow(r)) sum(r$pvalue < 0.05 & r$p.adjust < 0.05) else 0L
  data.table(sid = d$sid[1], M12_n_sig = n, M12_pos = as.integer(n > 0)) }
cl <- makeCluster(min(16, detectCores() - 2)); clusterExport(cl, c("dr", "m12_one"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(clusterProfiler) })))
R12 <- rbindlist(parLapply(cl, sigs, m12_one)); stopCluster(cl)
tm("M12 per-signature done")

## ================= Calls =================
S <- Reduce(function(a, b) merge(a, b, by = "sid"), list(meta, R3, R12))
fwrite(S, file.path(O, "m3_m12_per_signature.tsv"), sep = "\t")
full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
comps <- list(A_vs_fitsig_k = list(S[arm == "A" & kind == "real"], S[arm == "A" & kind == "null"]),
              A_vs_Ak_sizematched = list(S[arm == "A" & kind == "real"], S[arm == "Ak"]),
              B_vs_fitsig_k = list(S[arm == "B" & kind == "real"], S[arm == "B" & kind == "null"]),
              Afull_vs_Ak_sizematched = list(S[arm == "A" & kind == "real" & src %in% full_ids], S[arm == "Ak" & src %in% full_ids]))
D <- rbindlist(lapply(names(comps), function(nm) rbindlist(lapply(c("M3_pos", "M12_pos"), function(m) {
  re <- comps[[nm]][[1]][[m]]; nu <- comps[[nm]][[2]][[m]]
  xr <- sum(re); nr <- length(re); xn <- sum(nu); nn <- length(nu); ci <- binom.test(xn, nn)$conf.int
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(comparison = nm, module = m, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn,
             FPR_CP_lo = ci[1], FPR_CP_hi = ci[2], typeI_inflated = ci[1] > 0.075, fisher_OR = unname(ft$estimate), fisher_P = ft$p.value) }))))
fwrite(D, file.path(O, "m3_m12_decision.tsv"), sep = "\t")
print(D); tm("DONE")
