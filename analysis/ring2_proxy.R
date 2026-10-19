# 第 2 环：签名在测什么（预注册 N17_ring2_preregistration.md，结果出来后不改）
# 每个签名 × 每个队列：Spearman ρ（增殖/Stromal/Immune/OXPHOS）、代理 R²、外部队列嵌套 Cox LRT（基础 = 分期+增殖+Stromal+Immune+年龄+性别）
# A 组已发表签名（原文系数 × TCGA sd）/ B 组 8 个机制集（fitsig）vs 同尺寸随机 fitsig_k；分块落盘、可续跑
# R2_SMOKE=1 -> 只跑少量任务做冒烟测试
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("R2_SMOKE", "0") == "1"
O <- file.path(B, if (SMOKE) "N11_out/ring2_smoke" else "N11_out/ring2"); CH <- file.path(O, "chunks"); dir.create(CH, showWarnings = FALSE, recursive = TRUE)
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
RA <- if (SMOKE) 2L else 25L; RB <- if (SMOKE) 3L else 200L

## ---- 年龄 / 性别 / ESTIMATE（照抄 n16_gate_pilot.R）----
ids <- lapply(cohorts, function(c) colnames(c$z))
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
cohorts$TCGA$age <- as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids$TCGA, ph[[1]])])
cohorts$TCGA$sex <- toupper(substr(ph$gender[match(ids$TCGA, ph[[1]])], 1, 1))
su <- fread(text = readLines(gzfile(file.path(GE, "GSE14520/GSE14520_Extra_Supplement.txt.gz"))), data.table = FALSE)
su <- su[su[["Tissue Type"]] == "Tumor", ]; su <- su[!duplicated(su$Affy_GSM), ]
cohorts$GSE14520$age <- as.numeric(su$Age[match(ids$GSE14520, su$Affy_GSM)])
cohorts$GSE14520$sex <- toupper(substr(su$Gender[match(ids$GSE14520, su$Affy_GSM)], 1, 1))
ls2 <- rd(file.path(GE, "HCCDB18_LIRI-JP/HCCDB18.sample.txt")); lp2 <- rd(file.path(GE, "HCCDB18_LIRI-JP/HCCDB18.patient.txt"))
pid <- ls2$PATIENT_ID[match(ids$LIRI, ls2$SAMPLE_ID)]
cohorts$LIRI$age <- as.numeric(lp2$AGE[match(pid, lp2$PATIENT_ID)])
cohorts$LIRI$sex <- toupper(substr(lp2$GENDER[match(pid, lp2$PATIENT_ID)], 1, 1))
suppressPackageStartupMessages(library(estimate)); utils::data("common_genes", package = "estimate", envir = environment())
plat <- c(TCGA = "illumina", GSE14520 = "affymetrix", LIRI = "illumina")
for (k in names(cohorts)) {
  X <- cohorts[[k]]$x[, ids[[k]]]; X <- X[rowSums(X > 0) > 0, ]
  tf <- tempfile(fileext = ".txt"); tg <- tempfile(fileext = ".gct"); te <- tempfile(fileext = ".gct")
  write.table(data.frame(NAME = rownames(X), X, check.names = FALSE), tf, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(tf, tg, id = "GeneSymbol"); estimateScore(tg, te, platform = plat[[k]])
  es <- read.table(te, skip = 2, header = TRUE, sep = "\t", check.names = FALSE); em <- as.matrix(es[, -(1:2)]); rownames(em) <- es$NAME
  stopifnot(ncol(em) == length(ids[[k]]))
  cohorts[[k]]$stromal <- as.numeric(scale(em["StromalScore", ])); cohorts[[k]]$immune <- as.numeric(scale(em["ImmuneScore", ]))
}
D <- lapply(cohorts, function(c) { sx <- ifelse(c$sex %in% c("M", "F"), c$sex, NA)
  data.frame(stage = factor(c$stage, levels = c("I-II", "III-IV")), prolif = as.numeric(scale(c$prolif)), oxphos = as.numeric(scale(c$oxphos)),
             stromal = c$stromal, immune = c$immune, age = c$age, sex = factor(sx, levels = c("F", "M"))) })
fwrite(rbindlist(lapply(names(D), function(k) data.table(cohort = k, n = nrow(D[[k]]), events = sum(cohorts[[k]]$event),
  complete_base = sum(complete.cases(D[[k]]))))), file.path(O, "covariate_audit.tsv"), sep = "\t")

## ---- 签名 ----
L <- fread(file.path(B, "N12_signatures/signatures_long.tsv"), encoding = "UTF-8"); L[, coef := suppressWarnings(as.numeric(coef))]
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
L <- L[sig_id %in% I[group_primary == "是" & 是否有完整系数 == "是"]$sig_id & !is.na(coef)]
sdT <- apply(cohorts$TCGA$x[Uall, colnames(cohorts$TCGA$z)], 1, sd)
pub <- lapply(split(L, L$sig_id), function(d) { d <- d[gene_hgnc %in% Uall]; d <- d[!duplicated(gene_hgnc)]; setNames(d$coef * sdT[d$gene_hgnc], d$gene_hgnc) })
pub <- pub[vapply(pub, length, 1L) >= 2]
if (SMOKE) pub <- pub[1:3]

jobs <- rbind(CJ(arm = "A", src = names(pub), kind = c("real", "null"), rep = seq_len(RA))[!(kind == "real" & rep > 1)],
              CJ(arm = "Ak", src = names(pub), kind = "null", rep = seq_len(RA)),
              CJ(arm = "B", src = names(sets), kind = c("real", "null"), rep = seq_len(RB))[!(kind == "real" & rep > 1)])
jobs[kind == "real", rep := 0L]
jobs[, id := .I]; jobs[, chunk := (id - 1L) %/% 100L + 1L]
fwrite(jobs, file.path(O, "jobs.tsv"), sep = "\t")
cat("jobs:", nrow(jobs), " chunks:", max(jobs$chunk), "\n")

# Ak（预注册变更 10-15）：尺寸匹配对照，照抄 ring3_armA_null.R 的 fitsig_cap
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
BASE <- c("stage", "prolif", "stromal", "immune", "age", "sex")
lrt <- function(s, k, base) { d <- cbind(D[[k]], s = as.numeric(scale(s)))
  ok <- complete.cases(d[, c(base, "s")]); d <- d[ok, ]; y <- cohorts[[k]]$y[ok]
  m0 <- coxph(as.formula(paste("y ~", paste(base, collapse = "+"))), data = d)
  m1 <- coxph(as.formula(paste("y ~", paste(c(base, "s"), collapse = "+"))), data = d)
  c(p = anova(m0, m1)[["Pr(>|Chi|)"]][2], hr_adj = exp(coef(m1)[["s"]]), hr_unadj = exp(coef(coxph(y ~ s, data = d))[["s"]]), n = nrow(d)) }
run1 <- function(i) { j <- jobs[i]; b <- getsig(j)
  out <- lapply(names(cohorts), function(k) { s <- score(b, cohorts[[k]]$z); d <- D[[k]]
    r <- c(rho_prolif = cor(s, d$prolif, method = "spearman"), rho_stromal = cor(s, d$stromal, method = "spearman"),
           rho_immune = cor(s, d$immune, method = "spearman"), rho_oxphos = cor(s, d$oxphos, method = "spearman"),
           R2_proxy = summary(lm(s ~ prolif + stromal + immune + oxphos, data = d))$r.squared)
    if (k != "TCGA") { a <- lrt(s, k, BASE); a2 <- lrt(s, k, c(BASE, "oxphos"))
      r <- c(r, lrt_p = a[["p"]], hr_adj = a[["hr_adj"]], hr_unadj = a[["hr_unadj"]], n_lrt = a[["n"]], lrt_p_oxphos = a2[["p"]]) }
    data.table(j, cohort = k, k_model = length(b), t(r)) })
  rbindlist(out, fill = TRUE) }

cl <- makeCluster(min(20, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })))
clusterExport(cl, c("jobs", "getsig", "fitsig_cap", "run1", "lrt", "BASE", "D", "pub", "sets", "U", "cohorts", "fitsig", "score"))
t0 <- Sys.time()
for (c in sort(unique(jobs$chunk))) { f <- file.path(CH, sprintf("r2_%04d.tsv", c)); if (file.exists(f)) next
  r <- rbindlist(parLapply(cl, which(jobs$chunk == c), function(i) tryCatch(run1(i), error = function(e) data.table(jobs[i], err = conditionMessage(e)))), fill = TRUE)
  fwrite(r, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  if (c %% 5 == 0) cat("chunk", c, "/", max(jobs$chunk), round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n") }
stopCluster(cl)
R <- rbindlist(lapply(list.files(CH, "^r2_.*\\.tsv$", full.names = TRUE), fread), fill = TRUE)
fwrite(R, file.path(O, "ring2_rows.tsv"), sep = "\t")
cat("rows", nrow(R), " errors", if ("err" %in% names(R)) sum(!is.na(R$err)) else 0, "\nDONE\n")
