# N11 试跑：HCC 细胞死亡签名 vs 随机基因集 vs 增殖/OXPHOS 代理（TCGA-LIHC，单队列，只给用户看，不进稿）
# 预设判读规则（跑之前定）：
#   R1 签名测试集 C-index 低于随机同尺寸基因集零分布的第 95 百分位 -> "不比随机强"
#   R2 全队列多因素 Cox（签名 + 增殖 + OXPHOS + 分期）中签名 Wald P >= 0.05 -> "被代理解释"
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
set.seed(20261015)
T <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_out")
O <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_out")
O <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_out")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
cat("OUT=", O, "\n")

ex <- fread(file.path(T, "TCGA-LIHC_gene_expression.tsv.gz"), data.table = FALSE)
rownames(ex) <- ex[[1]]; ex <- as.matrix(ex[, -1])
ex <- ex[, substr(colnames(ex), 14, 15) == "01"]
sv <- fread(file.path(T, "TCGA-LIHC_survival.tsv"), data.table = FALSE)
ph <- fread(file.path(T, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
sv <- sv[sv$sample %in% colnames(ex) & !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS), ]
stcol <- grep("pathologic_stage", colnames(ph), value = TRUE)[1]
st <- ph[[stcol]][match(sv$sample, ph[[1]])]
sv$stage <- ifelse(grepl("IV|III", st), "III-IV", ifelse(grepl("II|I", st), "I-II", NA))
ex <- ex[, sv$sample]
ex <- ex[rowMeans(ex > 1) > 0.5 & apply(ex, 1, sd) > 0, ]
z <- t(scale(t(ex)))
cat("n_samples=", ncol(ex), " n_genes=", nrow(ex), " events=", sum(sv$OS), "\n")

sets <- list(
  cuproptosis = c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1","GLS","CDKN2A",
                  "NFE2L2","NLRP3","ATP7B","ATP7A","SLC31A1","LIPT2","DBT","GCSH","DLST"),
  disulfidptosis = c("GYS1","NDUFS1","OXSM","LRPPRC","NDUFA11","NUBPL","NCKAP1","RPN1","SLC3A2","SLC7A11"))
gmt <- strsplit(readLines(file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_geneset/h.all.v2024.1.Hs.symbols.gmt")), "	")
H <- setNames(lapply(gmt, function(x) x[-(1:2)]), vapply(gmt, `[`, "", 1))
hs <- function(n) intersect(H[[n]], rownames(z))
prolif <- colMeans(z[union(hs("HALLMARK_E2F_TARGETS"), hs("HALLMARK_G2M_CHECKPOINT")), ])
oxphos <- colMeans(z[hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION"), ])

tr <- sample(seq_len(ncol(z)), round(ncol(z) * 2 / 3)); te <- setdiff(seq_len(ncol(z)), tr)
y <- Surv(sv$OS.time, sv$OS)
cidx <- function(score, idx) unname(concordance(y[idx] ~ score, reverse = TRUE)$concordance)

# 模板流程：训练集单因素 Cox P<0.05 预筛 -> LASSO-Cox(cv) -> 风险分；若无基因入选则用全集均值 z
pipe <- function(g) {
  g <- intersect(g, rownames(z)); X <- t(z[g, , drop = FALSE])
  p <- apply(X[tr, , drop = FALSE], 2, function(v) summary(coxph(y[tr] ~ v))$coefficients[5])
  keep <- names(p)[p < 0.05]
  if (length(keep) >= 2) {
    fit <- tryCatch(cv.glmnet(X[tr, keep], y[tr], family = "cox", nfolds = 5), error = function(e) NULL)
    if (!is.null(fit)) {
      b <- as.matrix(coef(fit, s = "lambda.min"))[, 1]; b <- b[b != 0]
      if (length(b) > 0) return(list(score = drop(X[, names(b), drop = FALSE] %*% b), n = length(b), mode = "lasso"))
    }
  }
  if (length(keep) == 1) return(list(score = X[, keep] * sign(coef(coxph(y[tr] ~ X[tr, keep]))), n = 1, mode = "single"))
  s <- rowMeans(X); sg <- sign(coef(coxph(y[tr] ~ s[tr])))
  list(score = s * sg, n = length(g), mode = "mean")
}

res <- list(); nullres <- list()
NR <- 500
for (nm in names(sets)) {
  g <- intersect(sets[[nm]], rownames(z)); k <- length(g)
  obs <- pipe(g); c_obs <- cidx(obs$score[te], te)
  nulls <- vapply(seq_len(NR), function(i) { r <- pipe(sample(rownames(z), k)); cidx(r$score[te], te) }, numeric(1))
  nullres[[nm]] <- data.table(set = nm, iter = seq_len(NR), c_test = nulls)
  d <- data.frame(y = y, sig = as.numeric(scale(obs$score)), prolif = as.numeric(scale(prolif)),
                  oxphos = as.numeric(scale(oxphos)), stage = factor(sv$stage))
  m1 <- coxph(y ~ sig, data = d)
  m2 <- coxph(y ~ sig + prolif + oxphos + stage, data = d)
  s1 <- summary(m1)$coefficients; s2 <- summary(m2)$coefficients
  res[[nm]] <- data.table(set = nm, n_genes_in_set = k, n_genes_in_model = obs$n, mode = obs$mode,
    c_test_signature = c_obs, null_median = median(nulls), null_p95 = quantile(nulls, 0.95),
    empirical_p = (sum(nulls >= c_obs) + 1) / (NR + 1),
    R1_not_better_than_random = c_obs < quantile(nulls, 0.95),
    HR_uni = s1["sig", 2], P_uni = s1["sig", 5],
    HR_adj = s2["sig", 2], P_adj = s2["sig", 5],
    R2_explained_by_proxy = s2["sig", 5] >= 0.05,
    cor_sig_prolif = cor(d$sig, d$prolif), cor_sig_oxphos = cor(d$sig, d$oxphos))
}
# 代理本身在测试集的 C-index（不经任何训练，方向按训练集定）
px <- data.table(proxy = c("proliferation", "oxphos"),
  c_test = c(cidx(prolif[te] * sign(coef(coxph(y[tr] ~ prolif[tr]))), te),
             cidx(oxphos[te] * sign(coef(coxph(y[tr] ~ oxphos[tr]))), te)))
fwrite(rbindlist(res), file.path(O, "pilot_summary.tsv"), sep = "\t")
fwrite(rbindlist(nullres), file.path(O, "pilot_null_cindex.tsv"), sep = "\t")
fwrite(px, file.path(O, "pilot_proxy_cindex.tsv"), sep = "\t")
fwrite(data.table(n = ncol(z), events = sum(sv$OS), n_train = length(tr), n_test = length(te),
  stage_NA = sum(is.na(sv$stage)), NR = NR), file.path(O, "pilot_meta.tsv"), sep = "\t")
cat("DONE\n")
