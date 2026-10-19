# Arm A real-signature weights: reported coefficients act on log expression -> equivalent z-scale weight = coef x sd_gene (TCGA log expression)
# Main analysis uses these weights; the original z-scale weights (ring3_full results) are kept as a sensitivity analysis. M7 refits each time and is not rerun.
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/ring3_armA_realrescale"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); source(file.path(B, "analysis/ring3_lib.R"), encoding = "UTF-8")
P <- prep_ring3(B, TG)
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
L <- L[sig_id %in% I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id & !is.na(coef)]
sdT <- apply(cohorts$TCGA$x[Uall, colnames(cohorts$TCGA$z)], 1, sd)
pub <- lapply(split(L, L$sig_id), function(d) { d <- d[gene_hgnc %in% Uall]; d <- d[!duplicated(gene_hgnc)]; setNames(d$coef * sdT[d$gene_hgnc], d$gene_hgnc) })
pub <- pub[vapply(pub, length, 1L) >= 2]
cat("TCGA expression range (log scale check):", round(range(cohorts$TCGA$x[Uall, ]), 2), "\n")
rows <- lapply(names(pub), function(sid) { b <- pub[[sid]]
  r <- tryCatch(c(m0_external(b), m1_indep(b, P), m2_immune(b, P), m4_gsea(b, P, seed = 1), m5_tmb(b, P), m6_nomo(b)), error = function(e) c(err = 1))
  data.table(arm = "A", src = sid, kind = "real", rep = 0L, k = length(b), t(r)) })
R <- rbindlist(rows, fill = TRUE); fwrite(R, file.path(O, "results_m06_real_rescaled.tsv"), sep = "\t")
cat("rows", nrow(R), " errors", if ("err" %in% names(R)) sum(!is.na(R$err)) else 0, "\nDONE\n")
