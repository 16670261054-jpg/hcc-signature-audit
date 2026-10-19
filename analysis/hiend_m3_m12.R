# 高端模块 M3 药敏 + M12 DSigDB 药物富集 —— 正式（2026-10-15 夜）
# 被审计对象 = 课程模板写法，判据照抄模板（不是我们的判据）：
#   M3（课程 161/209/225 40.oncoPredict + 41.boxplot；已逐行对照 161Anoikis/40.oncoPredict/Anoikis40.oncoPredict.R）：
#       GDSC2_Res=exp(GDSC2_Res)；测试矩阵 rowMeans>0.5；calcPhenotype(batchCorrect='eb', powerTransformPhenotype=TRUE,
#       removeLowVaryingGenes=0.2, minNumSamples=10, removeLowVaringGenesFrom='rawData')，set.seed(12345)；高/低风险组逐药 Wilcoxon。
#       阳性判据按第 3 环预注册 N14：≥1 个药 BH-FDR<0.05（模板原文 p<0.001 未校正，另报）。
#       修坑145：测试矩阵若非 log 尺度先 log2(x+1)；修坑146：保留 DRUG_ID 列名，不去后缀；数值是预测 IC50（越高越耐药）。
#   M12（课程 236 25.drugEnrich / 228 34.drugEnrich）：enricher(genes, TERM2GENE=DSigDB[,1:2], minGSSize=10, maxGSSize=500,
#       pvalueCutoff=1, qvalueCutoff=1)，阳性 = ≥1 个药 pvalue<0.05 & p.adjust<0.05（模板原文）。
# 签名 = hiend_m8_formal.R 落盘的 signature_genes_all_arms.tsv（与第 2 环同一批）
# 判读：随机阳性率 Clopper–Pearson 95% CI 下界 > 0.075 → I 类错误失控；真实 vs 随机 Fisher 双侧
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/hiend_m3_m12"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
DL <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "external")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))
G[, sid := paste(arm, src, kind, rep, sep = "|")]
sigs <- split(G, G$sid)
meta <- unique(G[, .(sid, arm, src, kind, rep)])

## ================= M3 =================
pred_f <- file.path(O, "DrugPredictions_TCGA.tsv")
if (!file.exists(pred_f)) {
  ex_dir <- file.path(DL, "oncoPredict_osf/extracted")
  if (!file.exists(file.path(ex_dir, "GDSC2_Res.rds"))) {
    dir.create(ex_dir, showWarnings = FALSE)
    for (n in c("GDSC2_Expr (RMA Normalized and Log Transformed).rds", "GDSC2_Res.rds")) {
      unzip(file.path(DL, "oncoPredict_osf/DataFiles.zip"), files = paste0("DataFiles/DataFiles/Training Data/", n), exdir = ex_dir, junkpaths = TRUE)
      stopifnot(file.exists(file.path(ex_dir, n))) }
  }
  trainExpr <- readRDS(file.path(ex_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
  trainRes <- readRDS(file.path(ex_dir, "GDSC2_Res.rds")); trainRes <- exp(trainRes)   # 模板写法：GDSC2_Res 为 ln(IC50)，模板先 exp 再交给 powerTransform
  X <- cohorts$TCGA$x[, colnames(cohorts$TCGA$z)]
  scale_note <- sprintf("TCGA max=%.2f", max(X, na.rm = TRUE))
  if (max(X, na.rm = TRUE) > 50) X <- log2(X + 1)                                       # 坑145
  X <- X[rowMeans(X) > 0.5, ]                                                            # 模板 L(data[rowMeans(data)>0.5,])
  fwrite(data.table(note = scale_note, logged = max(cohorts$TCGA$x, na.rm = TRUE) > 50, train_max = max(trainExpr)), file.path(O, "m3_scale_check.tsv"), sep = "\t")
  # oncoPredict 1.3.1：默认 folder=FALSE 不写 DrugPredictions.csv（课程旧版会写），结果只在返回值里 → 直接接返回值
  set.seed(12345)
  pm <- oncoPredict::calcPhenotype(trainingExprData = trainExpr, trainingPtype = trainRes, testExprData = X,
                             batchCorrect = "eb", powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2,
                             minNumSamples = 10, printOutput = FALSE, removeLowVaringGenesFrom = "rawData")
  if (is.list(pm)) pm <- pm$DrugPredictions
  stopifnot(is.matrix(pm), nrow(pm) == ncol(X), !anyDuplicated(colnames(pm)))           # 坑146
  pr <- data.table(sample = rownames(pm), pm)
  fwrite(pr, pred_f, sep = "\t"); tm(sprintf("calcPhenotype done: %d samples x %d drugs", nrow(pr), ncol(pr) - 1))
}
pr <- fread(pred_f); smp <- pr[[1]]; DM <- as.matrix(pr[, -1]); rownames(DM) <- smp
Z <- cohorts$TCGA$z; stopifnot(all(colnames(Z) %in% smp)); DM <- DM[colnames(Z), ]
m3_one <- function(d) { b <- setNames(d$coef, d$gene); s <- drop(t(Z[names(b), , drop = FALSE]) %*% b)
  hi <- s > median(s); p <- apply(DM, 2, function(v) wilcox.test(v[hi], v[!hi], exact = FALSE)$p.value)
  q <- p.adjust(p, "BH")
  data.table(sid = d$sid[1], M3_n_fdr05 = sum(q < 0.05), M3_pos = as.integer(any(q < 0.05)), M3_n_p001_template = sum(p < 0.001)) }
cl <- makeCluster(min(20, detectCores() - 2)); clusterExport(cl, c("Z", "DM", "m3_one")); invisible(clusterEvalQ(cl, library(data.table)))
R3 <- rbindlist(parLapply(cl, sigs, m3_one)); stopCluster(cl)
tm("M3 per-signature done")

## ================= M12 DSigDB =================
suppressPackageStartupMessages(library(clusterProfiler))
dr <- fread(file.path(DL, "DSigDB/DSigDB_All_detailed.txt"), sep = "\t", quote = "", select = 1:2, col.names = c("Drug", "Gene"))
dr <- unique(dr[!is.na(Gene) & Gene != ""])
m12_one <- function(d) { kk <- suppressMessages(enricher(unique(d$gene), pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, TERM2GENE = dr))
  r <- if (is.null(kk)) data.frame() else as.data.frame(kk)
  n <- if (nrow(r)) sum(r$pvalue < 0.05 & r$p.adjust < 0.05) else 0L
  data.table(sid = d$sid[1], M12_n_sig = n, M12_pos = as.integer(n > 0)) }
cl <- makeCluster(min(16, detectCores() - 2)); clusterExport(cl, c("dr", "m12_one"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(clusterProfiler) })))
R12 <- rbindlist(parLapply(cl, sigs, m12_one)); stopCluster(cl)
tm("M12 per-signature done")

## ================= 判读 =================
S <- Reduce(function(a, b) merge(a, b, by = "sid"), list(meta, R3, R12))
fwrite(S, file.path(O, "m3_m12_per_signature.tsv"), sep = "\t")
full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
comps <- list(A_vs_fitsig_k = list(S[arm == "A" & kind == "real"], S[arm == "A" & kind == "null"]),
              A_vs_Ak_sizematched = list(S[arm == "A" & kind == "real"], S[arm == "Ak"]),
              B_vs_fitsig_k = list(S[arm == "B" & kind == "real"], S[arm == "B" & kind == "null"]),
              Afull_vs_Ak_sizematched = list(S[arm == "A" & kind == "real" & src %in% full_ids], S[arm == "Ak" & src %in% full_ids]))
D <- rbindlist(lapply(names(comps), function(nm) rbindlist(lapply(c("M3_pos", "M12_pos"), function(m) {
  re <- comps[[nm]][[1]][[m]]; nu <- comps[[nm]][[2]][[m]]
  xr <- sum(re); nr <- length(re); xn <- sum(nu); nn <- length(nu); ci <- binom.test(xn, nn)$conf.int
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(comparison = nm, module = m, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn,
             FPR_CP_lo = ci[1], FPR_CP_hi = ci[2], typeI_inflated = ci[1] > 0.075, fisher_OR = unname(ft$estimate), fisher_P = ft$p.value) }))))
fwrite(D, file.path(O, "m3_m12_decision.tsv"), sep = "\t")
print(D); tm("DONE")
