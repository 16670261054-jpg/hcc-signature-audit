# 第 3 环 A 组随机对照补跑：按预注册改为模板流程（fitsig）。原 ring3_full.R 中 A 组 null 用 coxfit（不筛选）= 偏离预注册
# R3_NULLMODE=fitsig_k（预注册主分析：k 个随机基因走 fitsig） | fitsig_pool100（敏感性：100 个随机基因走 fitsig）
# 判读规则见 N14_ring3_preregistration.md（已定稿，结果出来后不改）
# 设计：
#   A 已发表签名（149 个，原文系数）  vs  每个签名 25 个同尺寸随机集（随机基因 + TCGA 多因素 Cox 系数），M7 每个 10 个
#   B 8 个机制基因集（模板流程 fitsig） vs  每个 200 个同尺寸随机集（同样走 fitsig），M7 每个 25 个
# 分块落盘、可续跑：已存在的块文件直接跳过
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, Sys.getenv("R3_OUT", paste0("N11_out/ring3_armA_", Sys.getenv("R3_NULLMODE", "fitsig_k")))); CH <- file.path(O, "chunks"); dir.create(CH, showWarnings = FALSE, recursive = TRUE)
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); source(file.path(B, "analysis/ring3_lib.R"), encoding = "UTF-8")
RA <- as.integer(Sys.getenv("R3_RA", "25")); RA7 <- as.integer(Sys.getenv("R3_RA7", "10"))
RB <- as.integer(Sys.getenv("R3_RB", "200")); RB7 <- as.integer(Sys.getenv("R3_RB7", "25"))
NM <- Sys.getenv("R3_NULLMODE", "fitsig_k"); stopifnot(NM %in% c("fitsig_k", "fitsig_pool100", "fitsig_pool100_k"))
POOL <- NM %in% c("fitsig_pool100", "fitsig_pool100_k")
# fitsig_pool100_k（2026-10-15 尺寸匹配）：同 pool100 的 100 基因池与种子；lambda.min 模型超过 k 个基因时，
#   沿同一 LASSO 路径取非零系数数 <= k 的最小 lambda（glmnet 路径按自由度截断）；<= k 时与模板完全相同
fitsig_cap <- function(g, Z, y, k) {
  X <- t(Z[g, , drop = FALSE])
  p <- apply(X, 2, function(v) { s <- tryCatch(summary(coxph(y ~ v))$coefficients[5], error = function(e) 1); s })
  kp <- names(p)[p < 0.05]
  if (length(kp) >= 2) {
    f <- tryCatch(suppressWarnings(cv.glmnet(X[, kp], y, family = "cox", nfolds = 5)), error = function(e) NULL)
    if (!is.null(f)) { b <- as.matrix(coef(f, s = "lambda.min"))[, 1]; b <- b[b != 0]
      if (length(b) > k) { gf <- f$glmnet.fit; lam <- gf$lambda[max(which(gf$df <= k))]; b <- as.matrix(coef(gf, s = lam))[, 1]; b <- b[b != 0] }
      if (length(b)) return(b) }
  }
  if (length(kp) == 1) return(setNames(sign(coef(coxph(y ~ X[, kp]))), kp))
  sg <- sign(coef(coxph(y ~ rowMeans(X)))); setNames(rep(sg / length(g), length(g)), g)
}
P <- prep_ring3(B, TG)

## ---- 签名清单 ----
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
keep_ids <- I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id
L <- L[sig_id %in% keep_ids & !is.na(coef)]
cov <- L[, .(n_genes = .N, n_in = sum(gene_hgnc %in% Uall)), by = .(sig_id, mechanism)]
cov[, full := n_in == n_genes]
pub <- lapply(split(L, L$sig_id), function(d) { d <- d[gene_hgnc %in% Uall]; d <- d[!duplicated(gene_hgnc)]; setNames(d$coef, d$gene_hgnc) })
pub <- pub[vapply(pub, length, 1L) >= 2]
fwrite(cov, file.path(O, "published_coverage.tsv"), sep = "\t")

## ---- 任务表 ----
coxfit <- function(g) { X <- t(cohorts$TCGA$z[g, , drop = FALSE]); f <- tryCatch(coxph(cohorts$TCGA$y ~ X), error = function(e) NULL)
  b <- if (is.null(f)) rep(NA, length(g)) else unname(coef(f)); b[is.na(b)] <- 0; setNames(b, g) }
jobs <- CJ(arm = "A", src = names(pub), kind = "null", rep = seq_len(RA))
jobs[, id := .I]; jobs[, chunk := (id - 1L) %/% 100L + 1L]
fwrite(jobs, file.path(O, "jobs.tsv"), sep = "\t")
jobs7 <- CJ(arm = "A", src = names(pub), kind = "null", rep = seq_len(RA7))
jobs7[, id := .I]; jobs7[, chunk := (id - 1L) %/% 20L + 1L]
cat("jobs M0-M6:", nrow(jobs), " chunks:", max(jobs$chunk), " | jobs M7:", nrow(jobs7), " chunks:", max(jobs7$chunk), "\n")

# 每个任务的签名：确定性种子 -> 可复现、可续跑
getsig <- function(j) {
  seed <- 5e6 + 1e6 * POOL + 1e6 * (j$arm == "B") + 1e3 * match(j$src, if (j$arm == "A") names(pub) else names(sets)) + j$rep
  set.seed(seed)
  if (j$arm == "A") { k <- length(pub[[j$src]])
    if (j$kind == "real") return(list(b = pub[[j$src]], genes = names(pub[[j$src]]), seed = seed))
    g <- sample(U, if (POOL) 100L else k)
    b <- if (NM == "fitsig_pool100_k") fitsig_cap(g, cohorts$TCGA$z, cohorts$TCGA$y, k) else fitsig(g, cohorts$TCGA$z, cohorts$TCGA$y)
    return(list(b = b, genes = if (POOL) names(b) else g, seed = seed)) }
  k <- length(sets[[j$src]]); g <- if (j$kind == "real") sets[[j$src]] else sample(U, k)
  list(b = fitsig(g, cohorts$TCGA$z, cohorts$TCGA$y), genes = g, seed = seed)
}
run1 <- function(i) { j <- jobs[i]; s <- getsig(j); b <- s$b
  r <- tryCatch(c(m0_external(b), m1_indep(b, P), m2_immune(b, P), m4_gsea(b, P, seed = s$seed), m5_tmb(b, P), m6_nomo(b)),
                error = function(e) c(err = 1))
  data.table(j, k = length(s$genes), n_model = length(b), t(r)) }
run7 <- function(i) { j <- jobs7[i]; s <- getsig(j)
  o <- tryCatch(m7_multi(s$genes, seed = s$seed), error = function(e) NULL)
  if (is.null(o)) return(data.table(j, k = length(s$genes), err = 1))
  data.table(j, k = length(s$genes), t(o$summary), best = o$best) }

options(rf.cores = 1, mc.cores = 1)
cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(limma); library(fgsea); library(survcomp) }); options(rf.cores = 1, mc.cores = 1) }))
clusterExport(cl, c("jobs", "jobs7", "getsig", "run1", "run7", "coxfit", "pub", "sets", "U", "Uall", "cohorts", "P", "fitsig", "score", "cidx",
                    "NM", "POOL", "fitsig_cap", "m0_external", "m1_indep", "m2_immune", "m4_gsea", "m5_tmb", "m6_nomo", "m7_multi", "grp"))
t0 <- Sys.time()
for (c in sort(unique(jobs$chunk))) { f <- file.path(CH, sprintf("m06_%04d.tsv", c)); if (file.exists(f)) next
  fwrite(rbindlist(parLapply(cl, which(jobs$chunk == c), run1), fill = TRUE), paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  if (c %% 10 == 0) cat("m06 chunk", c, "/", max(jobs$chunk), round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n") }
for (c in sort(unique(jobs7$chunk))) { f <- file.path(CH, sprintf("m7_%04d.tsv", c)); if (file.exists(f)) next
  fwrite(rbindlist(parLapply(cl, which(jobs7$chunk == c), run7), fill = TRUE), paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  if (c %% 10 == 0) cat("m7 chunk", c, "/", max(jobs7$chunk), round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n") }
stopCluster(cl)
R06 <- rbindlist(lapply(list.files(CH, "^m06_.*\\.tsv$", full.names = TRUE), fread), fill = TRUE)
R7 <- rbindlist(lapply(list.files(CH, "^m7_.*\\.tsv$", full.names = TRUE), fread), fill = TRUE)
fwrite(R06, file.path(O, "results_m06.tsv"), sep = "\t"); fwrite(R7, file.path(O, "results_m7.tsv"), sep = "\t")
cat("rows m06", nrow(R06), "of", nrow(jobs), " errors", sum(!is.na(R06$err)), " | rows m7", nrow(R7), "of", nrow(jobs7), " errors", sum(!is.na(R7$err)), "\nDONE\n")
