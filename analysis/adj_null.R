# N11 第 1 步补：调整代理后的随机零分布（外部队列）
# 预设判读规则（开跑前写定）：
#   R3 某集合在某外部队列中，签名调整后 Wald |z|（Cox: sig + prolif + oxphos + stage）
#      < 同尺寸随机集同样调整后 |z| 的第 95 百分位 -> 该队列"残余信号不特异"
#   仅当签名在 >=1 个外部队列 R3 不成立时，才把该集合列为"可能有特异残余信号"的候选，进下游机制分析
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(parallel) })
set.seed(20261016)
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
O <- file.path(B, Sys.getenv("N11_OUT", "N11_out/robust")); dir.create(O, showWarnings = FALSE, recursive = TRUE)
TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")

Zt <- cohorts$TCGA$z; yt <- cohorts$TCGA$y
NE <- as.integer(Sys.getenv("N11_NE", "200"))
adjz <- function(s, C) {
  d <- data.frame(sig = as.numeric(scale(s)), prolif = as.numeric(scale(C$prolif)),
                  oxphos = as.numeric(scale(C$oxphos)), stage = factor(C$stage))
  m <- tryCatch(summary(coxph(C$y ~ sig + prolif + oxphos + stage, data = d))$coefficients, error = function(e) NULL)
  if (is.null(m)) return(c(z = NA, p = NA)) else c(z = unname(m["sig", 4]), p = unname(m["sig", 5]))
}
cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(survival); library(glmnet) }) }))
clusterExport(cl, c("fitsig", "score", "adjz", "Zt", "yt", "sets", "U", "cohorts", "NE"))
res <- parLapply(cl, names(sets), function(nm) {
  set.seed(77 + match(nm, names(sets)))
  g <- sets[[nm]]; k <- length(g); b <- fitsig(g, Zt, yt)
  nullb <- lapply(seq_len(NE), function(r) fitsig(sample(U, k), Zt, yt))
  out <- list()
  for (cn in c("GSE14520", "LIRI")) {
    C <- cohorts[[cn]]; o <- adjz(score(b, C$z), C)
    nz <- vapply(nullb, function(bb) abs(adjz(score(bb, C$z), C)["z"]), numeric(1))
    np <- vapply(nullb, function(bb) adjz(score(bb, C$z), C)["p"], numeric(1))
    out[[cn]] <- data.frame(set = nm, cohort = cn, k = k, z_adj = o["z"], P_adj = o["p"],
      null_absz_p95 = unname(quantile(nz, 0.95, na.rm = TRUE)), null_frac_Padj_lt05 = mean(np < 0.05, na.rm = TRUE),
      emp_p_absz = (sum(nz >= abs(o["z"]), na.rm = TRUE) + 1) / (sum(!is.na(nz)) + 1),
      R3_residual_not_specific = abs(o["z"]) < quantile(nz, 0.95, na.rm = TRUE))
  }
  do.call(rbind, out)
})
stopCluster(cl)
res <- rbindlist(res); fwrite(res, file.path(O, "adjusted_null.tsv"), sep = "\t")
print(res[, .(set, cohort, P_adj = signif(P_adj, 3), null_frac_Padj_lt05 = round(null_frac_Padj_lt05, 3), emp_p_absz = round(emp_p_absz, 3), R3_residual_not_specific)])
cat("DONE\n")
