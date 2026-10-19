# 第 3 环模块函数库：每个模块对一个签名（系数向量 b）返回"模板阳性条件"是否达标 + 关键统计量
# 依赖：load_cohorts.R 已 source（cohorts, U, sets, fitsig, score, cidx）
suppressPackageStartupMessages({ library(survival); library(limma); library(fgsea); library(survcomp) })

## ---- 一次性预计算（与签名无关）----
prep_ring3 <- function(B, TG) {
  P <- list()
  ph <- data.table::fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
  ids <- colnames(cohorts$TCGA$z)
  P$age <- as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids, ph[[1]])])
  P$sex <- factor(ph$gender[match(ids, ph[[1]])])
  # ssGSEA 29 个免疫集（课程 126 同款 immune.gmt）+ ESTIMATE 三分数：只算一次
  rg <- function(f) { s <- strsplit(readLines(f, warn = FALSE), "\t"); setNames(lapply(s, function(x) x[-(1:2)]), vapply(s, `[`, "", 1)) }
  imm <- rg(file.path(B, "N11_geneset/immune29_course126.gmt"))
  X <- cohorts$TCGA$x[, ids]; X <- X[rowSums(X > 0) > 0, ]
  par <- GSVA::ssgseaParam(X, imm, normalize = TRUE)
  P$imm <- GSVA::gsva(par, verbose = FALSE)
  tf <- tempfile(fileext = ".txt"); tg <- tempfile(fileext = ".gct"); te <- tempfile(fileext = ".gct")
  write.table(data.frame(NAME = rownames(X), X, check.names = FALSE), tf, sep = "\t", quote = FALSE, row.names = FALSE)
  suppressPackageStartupMessages(library(estimate)); utils::data("common_genes", package = "estimate", envir = environment())
  filterCommonGenes(tf, tg, id = "GeneSymbol")
  estimateScore(tg, te, platform = "illumina")
  es <- read.table(te, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  em <- as.matrix(es[, -(1:2)]); rownames(em) <- es$NAME; colnames(em) <- colnames(X)
  P$imm <- rbind(P$imm, em[c("StromalScore", "ImmuneScore", "ESTIMATEScore"), , drop = FALSE])
  # TMB：mc3 非同义突变数 / 38 Mb
  mf <- data.table::fread(file.path(TG, "LIHC_mc3.txt.gz"), data.table = FALSE)
  ns <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site")
  cnt <- table(mf$sample[mf$effect %in% ns]); sampled <- unique(mf$sample)
  P$tmb <- ifelse(ids %in% sampled, as.numeric(cnt[ids]) / 38, NA); P$tmb[ids %in% sampled & is.na(P$tmb)] <- 0
  # Hallmark
  P$H <- lapply(rg(file.path(B, "N11_geneset/h.all.v2024.1.Hs.symbols.gmt")), function(g) intersect(g, U))
  P
}

## ---- 各模块 ----
# M0 外部验证：外部队列 C-index 显著 > 0.5（survcomp::concordance.index 单侧检验）
m0_external <- function(b) {
  r <- lapply(c("GSE14520", "LIRI"), function(cn) { C <- cohorts[[cn]]; s <- score(b, C$z)
    ci <- concordance.index(s, C$time, C$event, method = "noether")
    c(c = ci$c.index, p = ci$p.value) })
  c(M0_c_GSE = r[[1]][["c"]], M0_p_GSE = r[[1]][["p"]], M0_c_LIRI = r[[2]][["c"]], M0_p_LIRI = r[[2]][["p"]],
    M0_pos = as.numeric(r[[1]][["p"]] < 0.05 & r[[2]][["p"]] < 0.05))
}
# M1 独立预后（TCGA，模板做法：风险分 + 年龄 + 性别 + 分期 多因素 Cox）
m1_indep <- function(b, P) {
  C <- cohorts$TCGA; d <- data.frame(s = as.numeric(scale(score(b, C$z))), age = P$age, sex = P$sex, stage = factor(C$stage))
  m <- summary(coxph(C$y ~ s + age + sex + stage, data = d))$coefficients
  c(M1_HR = m["s", 2], M1_p = m["s", 5], M1_pos = as.numeric(m["s", 5] < 0.05))
}
grp <- function(b) { s <- score(b, cohorts$TCGA$z); factor(ifelse(s > median(s), "high", "low"), levels = c("low", "high")) }
# M2 免疫浸润差异：32 个分数 Wilcoxon + BH，任一 FDR<0.05
m2_immune <- function(b, P) {
  g <- grp(b); p <- apply(P$imm, 1, function(v) wilcox.test(v[g == "high"], v[g == "low"], exact = FALSE)$p.value)
  q <- p.adjust(p, "BH"); c(M2_nsig = sum(q < 0.05), M2_pos = as.numeric(any(q < 0.05)))
}
# M4 GSEA：limma t 排序 -> fgsea Hallmark，任一 padj<0.05
m4_gsea <- function(b, P, seed = 1) {
  g <- grp(b); des <- model.matrix(~ g); fit <- eBayes(lmFit(cohorts$TCGA$z, des))
  t <- fit$t[, 2]; set.seed(seed); r <- suppressWarnings(fgsea(P$H, t, minSize = 15, maxSize = 500, eps = 1e-10, nproc = 1))
  c(M4_nsig = sum(r$padj < 0.05, na.rm = TRUE), M4_pos = as.numeric(any(r$padj < 0.05, na.rm = TRUE)))
}
# M5 TMB：a) 高低风险组 TMB Wilcoxon；b) TMB 中位 × 风险 四组 log-rank
m5_tmb <- function(b, P) {
  g <- grp(b); ok <- !is.na(P$tmb); C <- cohorts$TCGA
  pa <- wilcox.test(P$tmb[ok & g == "high"], P$tmb[ok & g == "low"], exact = FALSE)$p.value
  tg <- ifelse(P$tmb > median(P$tmb, na.rm = TRUE), "TMBhi", "TMBlo"); four <- interaction(g, tg)
  pb <- 1 - pchisq(survdiff(C$y[ok] ~ four[ok])$chisq, df = nlevels(droplevels(four[ok])) - 1)
  c(M5a_p = pa, M5a_pos = as.numeric(pa < 0.05), M5b_p = pb, M5b_pos = as.numeric(pb < 0.05))
}
# M6 列线图增量：外部队列 Cox(分期) vs Cox(分期 + 风险分) 嵌套似然比检验
m6_nomo <- function(b) {
  r <- vapply(c("GSE14520", "LIRI"), function(cn) { C <- cohorts[[cn]]; s <- score(b, C$z); ok <- !is.na(C$stage)
    d <- data.frame(s = s[ok], st = factor(C$stage[ok])); y <- C$y[ok]
    m0 <- coxph(y ~ st, data = d); m1 <- coxph(y ~ st + s, data = d)
    anova(m0, m1)[["Pr(>|Chi|)"]][2] }, numeric(1))
  c(M6_p_GSE = r[[1]], M6_p_LIRI = r[[2]], M6_pos = as.numeric(all(r < 0.05)))
}

# M7 多算法挑最优（113 组合精简版）：同一基因集喂 7 种生存算法，按两个外部队列平均 C 选第一名（模板做法）
#   阳性条件：选出模型在两个外部队列 C-index 都显著 > 0.5（survcomp）
m7_multi <- function(g, seed = 1) {
  set.seed(seed); g <- intersect(g, U); Xt <- t(cohorts$TCGA$z[g, , drop = FALSE]); yt <- cohorts$TCGA$y
  ext <- lapply(c("GSE14520", "LIRI"), function(cn) t(cohorts[[cn]]$z[g, , drop = FALSE]))
  pred <- list()
  for (a in c(0, 0.5, 1)) { f <- tryCatch(glmnet::cv.glmnet(Xt, yt, family = "cox", alpha = a, nfolds = 5), error = function(e) NULL)
    if (!is.null(f)) pred[[paste0("glmnet_a", a)]] <- local({ ff <- f; function(X) drop(predict(ff, X, s = "lambda.min")) }) }
  dt <- data.frame(Xt, check.names = TRUE); nm <- colnames(dt)
  sc <- tryCatch(MASS::stepAIC(coxph(yt ~ ., data = dt), direction = "both", trace = 0), error = function(e) NULL)
  if (!is.null(sc)) pred$stepCox <- function(X) predict(sc, newdata = setNames(data.frame(X, check.names = TRUE), nm), type = "lp")
  rf <- tryCatch(randomForestSRC::rfsrc(Surv(time, status) ~ ., data = cbind(dt, time = yt[, 1], status = yt[, 2]), ntree = 500, nodesize = 15, seed = -seed), error = function(e) NULL)
  if (!is.null(rf)) pred$RSF <- function(X) randomForestSRC::predict.rfsrc(rf, newdata = setNames(data.frame(X, check.names = TRUE), nm))$predicted
  gb <- tryCatch(suppressWarnings(gbm::gbm(yt ~ ., data = dt, distribution = "coxph", n.trees = 1000, interaction.depth = 3, shrinkage = 0.01, n.minobsinnode = 10, cv.folds = 5, n.cores = 1)), error = function(e) NULL)
  if (!is.null(gb)) { bi <- suppressWarnings(gbm::gbm.perf(gb, plot.it = FALSE, method = "cv"))
    pred$GBM <- function(X) predict(gb, newdata = setNames(data.frame(X, check.names = TRUE), nm), n.trees = bi, type = "link") }
  sp <- tryCatch({ d0 <- list(x = t(Xt), y = yt[, 1], censoring.status = yt[, 2], featurenames = colnames(Xt))
    fs <- superpc::superpc.train(d0, type = "survival")
    cv <- superpc::superpc.cv(fs, d0, n.threshold = 10, n.fold = 5, n.components = 1, min.features = 2, max.features = ncol(Xt))
    list(f = fs, d0 = d0, th = cv$thresholds[which.max(cv[["scor"]][1, ])]) }, error = function(e) NULL)
  if (!is.null(sp)) pred$superpc <- function(X) { d1 <- list(x = t(X), y = rep(1, nrow(X)), censoring.status = rep(1, nrow(X)), featurenames = colnames(X))
    drop(superpc::superpc.predict(sp$f, sp$d0, d1, threshold = sp$th, n.components = 1)$v.pred) }
  res <- rbindlist(lapply(names(pred), function(k) {
    cs <- vapply(1:2, function(j) { s <- tryCatch(pred[[k]](ext[[j]]), error = function(e) rep(NA, nrow(ext[[j]])))
      if (all(is.na(s))) return(c(NA, NA)); C <- cohorts[[c("GSE14520", "LIRI")[j]]]
      ci <- concordance.index(s, C$time, C$event, method = "noether"); c(ci$c.index, ci$p.value) }, numeric(2))
    data.table(algo = k, c_GSE = cs[1, 1], p_GSE = cs[2, 1], c_LIRI = cs[1, 2], p_LIRI = cs[2, 2]) }))
  res[, mean_c := (c_GSE + c_LIRI) / 2]; best <- res[which.max(mean_c)]
  lasso <- res[algo == "glmnet_a1"]
  list(table = res, summary = c(M7_n_algo = nrow(res), M7_best_mean_c = best$mean_c, M7_lasso_mean_c = if (nrow(lasso)) lasso$mean_c else NA,
       M7_pos = as.numeric(best$p_GSE < 0.05 & best$p_LIRI < 0.05)), best = best$algo)
}
