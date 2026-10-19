# M14 解释检查（预注册 N25_m14_prolif_check_prereg.md，结果出来后不改）
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
OUTD <- file.path(B, "N11_out/hiend_m14_prolif"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %5.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
O <- OUTD; source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8"); O <- OUTD
C <- cohorts$TCGA; ids <- colnames(C$z); prolif <- as.numeric(scale(C$prolif))

G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"), na.strings = "")
stopifnot(all(G$kind %in% c("real", "null")), any(G$kind == "null"))
G <- G[(arm == "A" & kind == "real") | (arm == "Ak" & rep <= 5)]
G[, sid := paste(arm, src, kind, rep, sep = "|")]
M14 <- fread(file.path(B, "N11_out/hiend_m13_m15/m13_m15_per_signature.tsv"), na.strings = "")[, .(sid, M14_p_old = M14_p)]
sigs <- split(G, G$sid); names(sigs) <- NULL; tm(sprintf("signatures %d", length(sigs)))

one <- function(d) {
  b <- setNames(d$coef, d$gene); out <- data.table(sid = d$sid[1], k = length(b), M14_p = NA_real_, prolif_sep = NA_real_, adj_p = NA_real_)
  if (length(b) < 2) return(out)
  X <- C$x[names(b), ids, drop = FALSE]
  rr <- tryCatch(suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(X, maxK = 9, reps = 50, pItem = 0.8, pFeature = 1,
          title = tempfile(), clusterAlg = "km", distance = "euclidean", seed = 123, plot = NULL, verbose = FALSE)), error = function(e) NULL)
  if (is.null(rr)) return(out)
  cl <- rr[[2]][["consensusClass"]][ids]; if (length(unique(cl)) != 2) return(out)
  g2 <- as.integer(cl == 2)
  auc <- wilcox.test(prolif[g2 == 1], prolif[g2 == 0], exact = FALSE)$statistic / (sum(g2 == 1) * sum(g2 == 0))
  out$M14_p <- 1 - pchisq(survdiff(C$y ~ cl)$chisq, df = 1)
  out$prolif_sep <- abs(unname(auc) - 0.5) * 2
  out$adj_p <- summary(coxph(C$y ~ factor(cl) + prolif))$coefficients[1, 5]
  out }
cl <- makeCluster(min(14, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(survival); library(data.table) })))
clusterExport(cl, c("C", "ids", "prolif", "one"))
R <- rbindlist(parLapply(cl, sigs, function(d) tryCatch(one(d), error = function(e) data.table(sid = d$sid[1], err = conditionMessage(e)))), fill = TRUE)
stopCluster(cl)
R <- merge(R, M14, by = "sid", all.x = TRUE); R[, arm := sub("\\|.*", "", sid)]
fwrite(R, file.path(OUTD, "m14_prolif_rows.tsv"), sep = "\t")
chk <- R[!is.na(M14_p) & !is.na(M14_p_old), .(n = .N, max_abs_diff_vs_M14 = max(abs(M14_p - M14_p_old)))]
fwrite(chk, file.path(OUTD, "m14_reproduce_check.tsv"), sep = "\t"); print(chk)

V <- R[!is.na(prolif_sep)]; a <- V[arm == "A"]; k <- V[arm == "Ak"]
mw <- wilcox.test(k$prolif_sep, a$prolif_sep, exact = FALSE)
cmp <- function(x, y) { xr <- sum(x < 0.05); nr <- length(x); xn <- sum(y < 0.05); nn <- length(y)
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2)); c(rate_pub = xr / nr, rate_Ak = xn / nn, OR = unname(ft$estimate), P = ft$p.value) }
raw <- cmp(a$M14_p, k$M14_p); adj <- cmp(a$adj_p, k$adj_p)
ca <- mw$p.value < 0.05 && median(k$prolif_sep) > median(a$prolif_sep); cb <- adj[["P"]] >= 0.05
D <- data.table(n_pub = nrow(a), n_Ak = nrow(k), prolif_sep_median_pub = median(a$prolif_sep), prolif_sep_median_Ak = median(k$prolif_sep), MWU_P = mw$p.value,
                raw_rate_pub = raw[["rate_pub"]], raw_rate_Ak = raw[["rate_Ak"]], raw_fisher_P = raw[["P"]],
                adj_rate_pub = adj[["rate_pub"]], adj_rate_Ak = adj[["rate_Ak"]], adj_OR = adj[["OR"]], adj_fisher_P = adj[["P"]],
                cond_a = ca, cond_b = cb, verdict = if (ca && cb) "supported" else if (ca) "partially_supported" else "not_supported")
fwrite(D, file.path(OUTD, "m14_prolif_decision.tsv"), sep = "\t"); print(D)
cat("errors", if ("err" %in% names(R)) sum(!is.na(R$err)) else 0, "\n"); tm("DONE")
