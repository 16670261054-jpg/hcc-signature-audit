# 补充模块 M13 免疫检查点 / M14 共识聚类分型 / M15 校准(ICI) + DCA（预注册 N24_extra_modules_prereg.md，结果出来后不改）
# 环境变量：M1315_SMOKE=1 → 每组少量签名，计时；M1315_MAXREP=N → A/Ak 随机组每签名只取前 N 个 rep（预注册允许的缩减）
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("M1315_SMOKE", "0") == "1"; MAXREP <- as.integer(Sys.getenv("M1315_MAXREP", "0"))
OUTD <- file.path(B, if (SMOKE) "N11_out/hiend_m13_m15_smoke" else "N11_out/hiend_m13_m15"); dir.create(file.path(OUTD, "chunks"), showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %6.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
O <- OUTD; source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); O <- OUTD

## ---- TCGA 临床（年龄 / 性别 / 分期）----
C <- cohorts$TCGA; ids <- colnames(C$z)
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
age <- as.numeric(ph$age_at_initial_pathologic_diagnosis[match(ids, ph[[1]])])
sex <- toupper(substr(ph$gender[match(ids, ph[[1]])], 1, 1)); sex[!sex %in% c("M", "F")] <- NA
CL <- data.frame(time = C$time, event = C$event, age = age, sex = factor(sex, levels = c("F", "M")), stage = factor(C$stage, levels = c("I-II", "III-IV")))
ok_cl <- complete.cases(CL); tm(sprintf("TCGA n=%d, complete clinical n=%d", nrow(CL), sum(ok_cl)))

## ---- 检查点基因（模板 gene.txt）----
gf <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "immune_checkpoint_genes.txt")
cp <- trimws(gsub("\r", "", readLines(gf))); cp <- cp[nzchar(cp)]
cp[cp == "C10orf54"] <- "VSIR"
cp_use <- intersect(cp, rownames(C$x)); fwrite(data.table(gene = cp, in_data = cp %in% rownames(C$x)), file.path(O, "checkpoint_genes.tsv"), sep = "\t")
tm(sprintf("checkpoint genes %d listed, %d in data", length(cp), length(cp_use)))
XCP <- C$x[cp_use, ids]

## ---- 签名 ----
G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"), na.strings = "")
stopifnot(all(G$kind %in% c("real", "null")), any(G$kind == "null"))
if (MAXREP > 0) G <- G[!(arm %in% c("A", "Ak") & kind == "null" & rep > MAXREP)]
G[, sid := paste(arm, src, kind, rep, sep = "|")]
if (SMOKE) { keep <- unique(rbind(G[arm == "A" & kind == "real", .(sid)][1:6], G[arm == "Ak", .(sid)][1:6], G[arm == "B" & kind == "real", .(sid)][1:2])$sid); G <- G[sid %in% keep] }
sigs <- split(G, G$sid); meta <- unique(G[, .(sid, arm, src, kind, rep)])
tm(sprintf("signatures %d", length(sigs)))

one <- function(d) {
  b <- setNames(d$coef, d$gene); s <- drop(t(C$z[names(b), ids, drop = FALSE]) %*% b)
  hi <- s > median(s); grpf <- factor(ifelse(hi, "high", "low"), levels = c("low", "high"))
  # M13
  pv <- apply(XCP, 1, function(v) if (sd(v) < 0.001) NA_real_ else wilcox.test(v ~ grpf, exact = FALSE)$p.value)
  m13 <- c(M13_n_sig = sum(pv < 0.05, na.rm = TRUE), M13_pos = as.integer(any(pv < 0.05, na.rm = TRUE)),
           M13_pos_BH = as.integer(any(p.adjust(pv, "BH") < 0.05, na.rm = TRUE)))
  # M14
  m14 <- c(M14_p = NA_real_, M14_pos = NA_real_, M14_skipped = 1)
  if (length(b) >= 2) {
    X <- C$x[names(b), ids, drop = FALSE]
    rr <- tryCatch(suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(X, maxK = 9, reps = 50, pItem = 0.8, pFeature = 1,
                     title = tempfile(), clusterAlg = "km", distance = "euclidean", seed = 123, plot = NULL, verbose = FALSE)), error = function(e) NULL)
    if (!is.null(rr)) { cl <- rr[[2]][["consensusClass"]][ids]
      if (length(unique(cl)) == 2) { p <- 1 - pchisq(survdiff(C$y ~ cl)$chisq, df = 1); m14 <- c(M14_p = p, M14_pos = as.integer(p < 0.05), M14_skipped = 0) } }
  }
  # M15
  m15 <- c(ICI_3y = NA_real_, DCA_pos = NA_real_, DCA_frac_thr = NA_real_)
  dd <- cbind(CL, risk = grpf)[ok_cl, ]
  r15 <- tryCatch({
    y <- Surv(dd$time, dd$event); u <- 36
    pr <- function(f) 1 - summary(survfit(f, newdata = dd), times = u, extend = TRUE)$surv[1, ]
    fn <- coxph(y ~ risk + age + sex + stage, data = dd)
    p_n <- pr(fn); p_r <- pr(coxph(y ~ risk, data = dd)); p_a <- pr(coxph(y ~ age, data = dd))
    p_s <- pr(coxph(y ~ sex, data = dd)); p_t <- pr(coxph(y ~ stage, data = dd))
    # ICI（Austin 2020）
    cll <- log(-log(1 - pmin(pmax(p_n, 1e-8), 1 - 1e-8)))
    fc <- rms::cph(y ~ rms::rcs(cll, 3), x = TRUE, y = TRUE, surv = TRUE)
    obs <- 1 - rms::survest(fc, times = u, newdata = data.frame(cll = cll))$surv
    ici <- mean(abs(p_n - obs))
    df <- data.frame(time = dd$time, event = dd$event, nomo = p_n, risk = p_r, age = p_a, sex = p_s, stage = p_t)
    dc <- dcurves::dca(Surv(time, event) ~ nomo + risk + age + sex + stage, data = df, time = u, thresholds = seq(0.01, 0.99, 0.01))
    nb <- as.data.table(dc$dca)[, .(variable, threshold, net_benefit)]
    w <- dcast(nb, threshold ~ variable, value.var = "net_benefit")
    best_other <- pmax(w$all, 0, w$age, w$sex, w$stage, na.rm = TRUE)
    win <- w$nomo > best_other
    c(ICI_3y = ici, DCA_pos = as.integer(any(win, na.rm = TRUE)), DCA_frac_thr = mean(win, na.rm = TRUE))
  }, error = function(e) m15)
  data.table(sid = d$sid[1], k = length(b), t(m13), t(m14), t(r15))
}

cl <- makeCluster(min(14, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table) })))
clusterExport(cl, c("C", "ids", "CL", "ok_cl", "XCP", "one"))
names(sigs) <- NULL; ch <- split(seq_along(sigs), ceiling(seq_along(sigs) / 200))
for (i in seq_along(ch)) { f <- file.path(OUTD, "chunks", sprintf("c_%04d.tsv", i)); if (file.exists(f)) next
  r <- rbindlist(parLapply(cl, sigs[ch[[i]]], function(d) tryCatch(one(d), error = function(e) data.table(sid = d$sid[1], err = conditionMessage(e)))), fill = TRUE)
  fwrite(r, paste0(f, ".tmp"), sep = "\t"); file.rename(paste0(f, ".tmp"), f)
  tm(sprintf("chunk %d/%d", i, length(ch))) }
stopCluster(cl)
R <- rbindlist(lapply(list.files(file.path(OUTD, "chunks"), "^c_.*\\.tsv$", full.names = TRUE), fread, na.strings = ""), fill = TRUE)
S <- merge(meta, R, by = "sid"); fwrite(S, file.path(OUTD, "m13_m15_per_signature.tsv"), sep = "\t")
tm(sprintf("rows %d, errors %d", nrow(S), if ("err" %in% names(S)) sum(!is.na(S$err)) else 0))

## ---- 判读 ----
full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
comps <- list(A_vs_fitsig_k = list(S[arm == "A" & kind == "real"], S[arm == "A" & kind == "null"]),
              A_vs_Ak_sizematched = list(S[arm == "A" & kind == "real"], S[arm == "Ak"]),
              B_vs_fitsig_k = list(S[arm == "B" & kind == "real"], S[arm == "B" & kind == "null"]),
              Afull_vs_Ak_sizematched = list(S[arm == "A" & kind == "real" & src %in% full_ids], S[arm == "Ak" & src %in% full_ids]))
Dd <- rbindlist(lapply(names(comps), function(nm) rbindlist(lapply(c("M13_pos", "M13_pos_BH", "M14_pos", "DCA_pos"), function(m) {
  re <- na.omit(comps[[nm]][[1]][[m]]); nu <- na.omit(comps[[nm]][[2]][[m]])
  if (!length(re) || !length(nu)) return(NULL)
  xr <- sum(re); nr <- length(re); xn <- sum(nu); nn <- length(nu); ci <- binom.test(xn, nn)$conf.int
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(comparison = nm, module = m, n_real = nr, rate_real = xr / nr, n_null = nn, FPR_null = xn / nn, FPR_CP_lo = ci[1], FPR_CP_hi = ci[2],
             typeI_inflated = ci[1] > 0.075, fisher_OR = unname(ft$estimate), fisher_P = ft$p.value) }))))
Di <- rbindlist(lapply(names(comps), function(nm) { a <- na.omit(comps[[nm]][[1]]$ICI_3y); b <- na.omit(comps[[nm]][[2]]$ICI_3y)
  if (length(a) < 2 || length(b) < 2) return(NULL)
  data.table(comparison = nm, n_real = length(a), n_null = length(b), ICI_median_real = median(a), ICI_median_null = median(b), MWU_P = wilcox.test(a, b, exact = FALSE)$p.value) }))
fwrite(Dd, file.path(OUTD, "m13_m15_decision.tsv"), sep = "\t"); fwrite(Di, file.path(OUTD, "m15_ici_compare.tsv"), sep = "\t")
fwrite(S[, .(n = .N, M14_skipped = sum(M14_skipped == 1, na.rm = TRUE), ICI_NA = sum(is.na(ICI_3y)), DCA_NA = sum(is.na(DCA_pos))), by = .(arm, kind)],
       file.path(OUTD, "m13_m15_missing.tsv"), sep = "\t")
print(Dd); print(Di); tm("DONE")
