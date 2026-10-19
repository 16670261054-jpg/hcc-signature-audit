suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
set.seed(1)
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA"); O <- file.path(B, "N11_out/ring3_pilot")
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); source(file.path(B, "analysis/ring3_lib.R"), encoding = "UTF-8")
rg <- function(f) { s <- strsplit(readLines(f, warn = FALSE), "\t"); setNames(lapply(s, function(x) x[-(1:2)]), vapply(s, `[`, "", 1)) }
H <- lapply(rg(file.path(B, "N11_geneset/h.all.v2024.1.Hs.symbols.gmt")), function(g) intersect(g, U))
b <- fitsig(sets$cuproptosis_Tsvetkov19, cohorts$TCGA$z, cohorts$TCGA$y); g <- grp(b)
t0 <- Sys.time(); fit <- eBayes(lmFit(cohorts$TCGA$z, model.matrix(~ g))); tt <- fit$t[, 2]; cat("limma", difftime(Sys.time(), t0, units = "secs"), "\n")
for (i in 1:3) { t0 <- Sys.time(); r <- suppressWarnings(fgsea(H, tt, minSize = 15, maxSize = 500)); cat("fgsea default", i, difftime(Sys.time(), t0, units = "secs"), sum(r$padj < 0.05), "\n") }
for (i in 1:2) { t0 <- Sys.time(); r2 <- suppressWarnings(fgsea(H, tt, minSize = 15, maxSize = 500, eps = 1e-10, nproc = 1)); cat("fgsea eps1e-10 nproc1", i, difftime(Sys.time(), t0, units = "secs"), sum(r2$padj < 0.05), "\n") }
cat("same_sig_calls", identical(r$padj < 0.05, r2$padj < 0.05), "\nDONE\n")
