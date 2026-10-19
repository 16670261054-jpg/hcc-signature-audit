# Exploratory: quick description of ZNF180, the only gene passing the genome-wide gates: correlation with proxies, single-cell detection (annotation only, not a gate)
suppressPackageStartupMessages({ library(survival); library(data.table); library(Matrix) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/n16_gate")
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
g <- "ZNF180"
r <- rbindlist(lapply(names(cohorts), function(k) { c <- cohorts[[k]]; v <- c$z[g, ]
  data.table(cohort = k, mean_expr = mean(c$x[g, ]), r_prolif = cor(v, c$prolif, method = "spearman"), r_oxphos = cor(v, c$oxphos, method = "spearman"),
             p_prolif = cor.test(v, c$prolif, method = "spearman", exact = FALSE)$p.value) }))
fwrite(r, file.path(O, "znf180_proxy_cor.tsv"), sep = "\t"); print(r)
D <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GSE149614")
X <- readRDS(file.path(D, "GSE149614_counts_sparse.rds")); M <- fread(file.path(D, "GSE149614_HCC.metadata.updated.txt.gz"))
stopifnot(g %in% rownames(X)); v <- X[g, M$Cell]
sc <- M[, .(n_cells = .N, detect_rate = mean(v[.I] > 0), mean_umi = mean(v[.I])), by = .(site, celltype)][order(site, celltype)]
fwrite(sc, file.path(O, "znf180_sc_detect.tsv"), sep = "\t"); print(sc[site %in% c("Tumor", "Normal")])
