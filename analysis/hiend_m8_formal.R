# 高端模块 M8 单细胞定位 —— 正式查表（2026-10-15 夜）
# 逐基因结果来自 hiend_m8_sc_pilot.R（课程 265 12.Seurat 写法，判据照抄被审计模板）：
#   M8a：FindAllMarkers |avg_log2FC|>1 & p_val_adj<0.05，签名 ≥1 个基因命中即阳性
#   M8b：各细胞类型 Tumor vs Normal，|avg_log2FC|>1 & p_val_adj<0.05，≥1 个基因、≥1 个类型即阳性
# 签名清单与第 2 环完全相同（同 getsig、同种子）：A 真实 / A 随机 fitsig_k / Ak 尺寸匹配 / B 真实 / B 随机
#   签名基因 = 最终模型中系数非零的基因（已发表签名 = 原文全部基因）
# 判读（第 3 环已定规则）：随机阳性率 Clopper–Pearson 95% CI 下界 > 0.075 → I 类错误失控（Bradley 1978 宽松标准）；真实 vs 随机 Fisher 双侧
# 另落盘全部签名基因清单（后续 M3 / DSigDB / 对接复用）
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/hiend_m8"); M8 <- O
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
O <- M8   # load_cohorts.R 会改写 O

## ---- 签名（照抄 ring2_proxy.R）----
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
L <- L[sig_id %in% I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id & !is.na(coef)]
sdT <- apply(cohorts$TCGA$x[Uall, colnames(cohorts$TCGA$z)], 1, sd)
pub <- lapply(split(L, L$sig_id), function(d) { d <- d[gene_hgnc %in% Uall]; d <- d[!duplicated(gene_hgnc)]; setNames(d$coef * sdT[d$gene_hgnc], d$gene_hgnc) })
pub <- pub[vapply(pub, length, 1L) >= 2]
jobs <- fread(file.path(B, "N11_out/ring2/jobs.tsv"))
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
getsig <- function(j) {
  idx <- match(j$src, if (j$arm %in% c("A", "Ak")) names(pub) else names(sets))
  seed <- c(A = 5e6, B = 6e6, Ak = 7e6)[[j$arm]] + 1e3 * idx + j$rep; set.seed(seed)
  if (j$arm == "Ak") return(fitsig_cap(sample(U, 100L), cohorts$TCGA$z, cohorts$TCGA$y, length(pub[[j$src]])))
  if (j$arm == "A" && j$kind == "real") return(pub[[j$src]])
  k <- if (j$arm == "A") length(pub[[j$src]]) else length(sets[[j$src]])
  g <- if (j$arm == "B" && j$kind == "real") sets[[j$src]] else sample(U, k)
  fitsig(g, cohorts$TCGA$z, cohorts$TCGA$y)
}
cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })))
clusterExport(cl, c("jobs", "getsig", "fitsig_cap", "pub", "sets", "U", "cohorts", "fitsig"))
G <- rbindlist(parLapply(cl, seq_len(nrow(jobs)), function(i) { j <- jobs[i]; b <- getsig(j)
  data.table(arm = j$arm, src = j$src, kind = j$kind, rep = j$rep, gene = names(b), coef = unname(b)) }))
stopCluster(cl)
fwrite(G, file.path(O, "signature_genes_all_arms.tsv"), sep = "\t")

## ---- 查表 ----
mk <- fread(file.path(O, "m8a_findallmarkers_all.tsv")); gd <- fread(file.path(O, "m8b_tumor_vs_normal_all.tsv"))
pos_a <- unique(mk[abs(avg_log2FC) > 1 & p_val_adj < 0.05]$gene); pos_b <- unique(gd[abs(avg_log2FC) > 1 & p_val_adj < 0.05]$gene)
sc_genes <- unique(c(mk$gene, gd$gene))
S <- G[, .(k = .N, n_in_sc_tested = sum(gene %in% sc_genes), n_pos_a = sum(gene %in% pos_a), n_pos_b = sum(gene %in% pos_b)), by = .(arm, src, kind, rep)]
S[, `:=`(M8a_pos = as.integer(n_pos_a >= 1), M8b_pos = as.integer(n_pos_b >= 1))]
fwrite(S, file.path(O, "m8_per_signature.tsv"), sep = "\t")

full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
comps <- list(A_vs_fitsig_k = list(S[arm == "A" & kind == "real"], S[arm == "A" & kind == "null"]),
              A_vs_Ak_sizematched = list(S[arm == "A" & kind == "real"], S[arm == "Ak"]),
              B_vs_fitsig_k = list(S[arm == "B" & kind == "real"], S[arm == "B" & kind == "null"]),
              Afull_vs_Ak_sizematched = list(S[arm == "A" & kind == "real" & src %in% full_ids], S[arm == "Ak" & src %in% full_ids]))
D <- rbindlist(lapply(names(comps), function(nm) rbindlist(lapply(c("M8a_pos", "M8b_pos"), function(m) {
  re <- comps[[nm]][[1]][[m]]; nu <- comps[[nm]][[2]][[m]]
  xr <- sum(re); nr <- length(re); xn <- sum(nu); nn <- length(nu); ci <- binom.test(xn, nn)$conf.int
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(comparison = nm, module = m, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn,
             FPR_CP_lo = ci[1], FPR_CP_hi = ci[2], typeI_inflated = ci[1] > 0.075, fisher_OR = unname(ft$estimate), fisher_P = ft$p.value,
             k_median_real = median(comps[[nm]][[1]]$k), k_median_null = median(comps[[nm]][[2]]$k)) }))))
fwrite(D, file.path(O, "m8_decision.tsv"), sep = "\t")
fwrite(data.table(p_gene_M8a_in_U = mean(U %in% pos_a), p_gene_M8b_in_U = mean(U %in% pos_b), U_n = length(U), U_in_sc = sum(U %in% sc_genes)),
       file.path(O, "m8_gene_level_rates.tsv"), sep = "\t")
print(D); cat("DONE\n")
