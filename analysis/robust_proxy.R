# N11 稳健性检验：8 个细胞死亡基因集 × 50 次 TCGA 拆分 × 随机零分布 + GSE14520 / LIRI-JP 外部验证 + 增殖/OXPHOS 代理调整
# 只做可行性判断，不进稿。
# ---- 预设判读规则（开跑前写定，结果出来后不改）----
#   R1  内部：50 次拆分中，签名测试集 C-index 相对同尺寸随机零分布的经验 P，取中位数；中位 P >= 0.05 -> "内部不比随机强"
#   R1e 外部：全 TCGA 训练的签名在外部队列的 C-index < 外部随机零分布第 95 百分位 -> 该队列"不比随机强"
#   R2  外部：Cox(签名 + 增殖 + OXPHOS + 分期) 中签名 Wald P >= 0.05 -> 该队列"被代理解释"
#   去留：8 个集合中 >= 5 个满足 R1，且其中多数在两个外部队列里也满足 R1e 或 R2 -> "代理"故事成立，值得做；
#         否则（混合/摇摆）-> 收手
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(parallel) })
set.seed(20261015)
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
O <- file.path(B, Sys.getenv("N11_OUT", "N11_out/robust")); dir.create(O, showWarnings = FALSE, recursive = TRUE)
TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
source(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "analysis/load_cohorts.R"), encoding = "UTF-8")

## ---- 1) 内部：50 次拆分 × (签名 + 100 随机) ----
S <- as.integer(Sys.getenv("N11_S", "50")); NR <- as.integer(Sys.getenv("N11_NR", "100"))
Zt <- cohorts$TCGA$z; yt <- cohorts$TCGA$y; nT <- ncol(Zt)
splits <- lapply(seq_len(S), function(i) sample(nT, round(nT * 2 / 3)))
cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(survival); library(glmnet) }) }))
jobs <- expand.grid(set = names(sets), split = seq_len(S), stringsAsFactors = FALSE)
clusterExport(cl, c("jobs", "fitsig", "score", "cidx", "Zt", "yt", "sets", "splits", "U", "NR", "nT"))
inner <- parLapply(cl, seq_len(nrow(jobs)), function(j) {
  nm <- jobs$set[j]; i <- jobs$split[j]; set.seed(1000 * i + match(nm, names(sets)))
  tr <- splits[[i]]; te <- setdiff(seq_len(nT), tr); g <- sets[[nm]]; k <- length(g)
  b <- fitsig(g, Zt[, tr], yt[tr]); co <- cidx(score(b, Zt[, te]), yt[te])
  nl <- vapply(seq_len(NR), function(r) { bb <- fitsig(sample(U, k), Zt[, tr], yt[tr]); cidx(score(bb, Zt[, te]), yt[te]) }, numeric(1))
  data.frame(set = nm, split = i, k = k, n_model = length(b), c_obs = co, null_med = median(nl),
             null_p95 = unname(quantile(nl, 0.95)), emp_p = (sum(nl >= co) + 1) / (NR + 1),
             genes = paste(names(b), collapse = ";"))
})
inner <- rbindlist(inner); fwrite(inner, file.path(O, "inner_splits.tsv"), sep = "\t")

## ---- 2) 外部：全 TCGA 训练 -> GSE14520 / LIRI；外部零分布 200 随机 ----
NE <- as.integer(Sys.getenv("N11_NE", "200"))
clusterExport(cl, c("cohorts", "NE"))
ext <- parLapply(cl, names(sets), function(nm) {
  set.seed(7 + match(nm, names(sets)))
  g <- sets[[nm]]; k <- length(g); b <- fitsig(g, Zt, yt)
  out <- list()
  nullb <- lapply(seq_len(NE), function(r) fitsig(sample(U, k), Zt, yt))
  for (cn in c("GSE14520", "LIRI")) {
    C <- cohorts[[cn]]; s <- score(b, C$z); co <- cidx(s, C$y)
    nl <- vapply(nullb, function(bb) cidx(score(bb, C$z), C$y), numeric(1))
    d <- data.frame(sig = as.numeric(scale(s)), prolif = as.numeric(scale(C$prolif)), oxphos = as.numeric(scale(C$oxphos)), stage = factor(C$stage))
    m1 <- summary(coxph(C$y ~ sig, data = d))$coefficients
    m2 <- summary(coxph(C$y ~ sig + prolif + oxphos + stage, data = d))$coefficients
    out[[cn]] <- data.frame(set = nm, cohort = cn, k = k, n_model = length(b), genes = paste(names(b), collapse = ";"),
      c_ext = co, null_med = median(nl), null_p95 = unname(quantile(nl, 0.95)), emp_p = (sum(nl >= co) + 1) / (NE + 1),
      HR_uni = m1["sig", 2], P_uni = m1["sig", 5], HR_adj = m2["sig", 2], P_adj = m2["sig", 5],
      n_adj = sum(complete.cases(d)), r_prolif = cor(d$sig, d$prolif), r_oxphos = cor(d$sig, d$oxphos))
  }
  do.call(rbind, out)
})
stopCluster(cl)
ext <- rbindlist(ext); fwrite(ext, file.path(O, "external.tsv"), sep = "\t")

## ---- 3) 代理本身的外部 C-index（方向按全 TCGA 定）----
px <- rbindlist(lapply(c("prolif", "oxphos"), function(p) {
  sg <- sign(coef(coxph(yt ~ cohorts$TCGA[[p]])))
  rbindlist(lapply(c("GSE14520", "LIRI"), function(cn) data.table(proxy = p, cohort = cn, c = cidx(sg * cohorts[[cn]][[p]], cohorts[[cn]]$y))))
}))
fwrite(px, file.path(O, "proxy_external.tsv"), sep = "\t")

## ---- 4) 按预设规则汇总 ----
sm <- inner[, .(k = k[1], c_obs_median = median(c_obs), null_med_median = median(null_med),
                frac_beat_p95 = mean(c_obs > null_p95), emp_p_median = median(emp_p)), by = set]
sm[, R1_inner_not_better := emp_p_median >= 0.05]
e2 <- dcast(ext, set ~ cohort, value.var = c("c_ext", "null_p95", "P_adj"))
e2[, R1e_GSE := c_ext_GSE14520 < null_p95_GSE14520][, R1e_LIRI := c_ext_LIRI < null_p95_LIRI]
e2[, R2_GSE := P_adj_GSE14520 >= 0.05][, R2_LIRI := P_adj_LIRI >= 0.05]
e2[, ext_proxy_both := (R1e_GSE | R2_GSE) & (R1e_LIRI | R2_LIRI)]
sm <- merge(sm, e2, by = "set")
fwrite(sm, file.path(O, "decision_summary.tsv"), sep = "\t")
dec <- data.table(n_sets = nrow(sm), n_R1 = sum(sm$R1_inner_not_better),
                  n_R1_and_ext = sum(sm$R1_inner_not_better & sm$ext_proxy_both))
dec[, go := n_R1 >= 5 & n_R1_and_ext > n_R1 / 2]
fwrite(dec, file.path(O, "decision.tsv"), sep = "\t"); print(dec)
cat("WARNINGS_N=", length(warnings()), "\nDONE\n")
