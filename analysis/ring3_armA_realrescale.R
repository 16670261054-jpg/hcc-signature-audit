# A 组真实签名权重校正：原文系数作用于 log 表达尺度 -> 在 z 尺度上等价权重 = coef × sd_gene（TCGA log 表达）
# 主分析改用此权重；原 z 尺度权重（ring3_full 中的结果）作为敏感性分析保留。M7 每次重新拟合，不受影响，不重跑。
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
