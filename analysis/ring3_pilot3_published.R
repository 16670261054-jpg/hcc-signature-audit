# Exploratory run of the published-signature path: (1) gene coverage (common universe of the three cohorts) (2) M0-M6 on 5 signatures with reported coefficients, NA check
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/ring3_pilot")
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); source(file.path(B, "analysis/ring3_lib.R"), encoding = "UTF-8")
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8")
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
L[, coef := suppressWarnings(as.numeric(coef))]
cov <- L[, .(n_genes = .N, n_coef_NA = sum(is.na(coef)), n_in_U = sum(gene_hgnc %in% U),
             missing = paste(gene_hgnc[!gene_hgnc %in% U], collapse = ";")), by = .(sig_id, PMID, mechanism)]
cov[, frac := n_in_U / n_genes]
fwrite(cov, file.path(O, "published_coverage.tsv"), sep = "\t")
cat("signatures", nrow(cov), " full_coverage", sum(cov$frac == 1 & cov$n_coef_NA == 0), " frac>=0.8", sum(cov$frac >= 0.8), "\n")
print(table(cut(cov$frac, c(-0.01, 0.5, 0.8, 0.99, 1))))
P <- prep_ring3(B, TG)
pick <- head(cov[frac == 1 & n_coef_NA == 0]$sig_id, 5); rows <- list()
for (sid in pick) {
  d <- L[sig_id == sid & !is.na(coef)]; b <- setNames(d$coef, d$gene_hgnc); b <- b[names(b) %in% U]
  r <- c(m0_external(b), m1_indep(b, P), m2_immune(b, P), m4_gsea(b, P), m5_tmb(b, P), m6_nomo(b))
  rows[[sid]] <- data.table(sig_id = sid, n = length(b), t(r))
}
res <- rbindlist(rows, fill = TRUE); fwrite(res, file.path(O, "pilot3_published_modules.tsv"), sep = "\t"); print(res)
na_cols <- names(res)[vapply(res, function(x) any(is.na(x)), logical(1))]
cat("NA_COLUMNS:", if (length(na_cols)) paste(na_cols, collapse = ",") else "none", "\nDONE\n")
