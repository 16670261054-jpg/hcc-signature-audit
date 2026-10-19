# Proxies recomputed after removing signature genes (decision rule written before running)
# Candidate signatures and gate-1 references are produced by the template pipeline on the simulated TCGA outcome; gates pooling the two external cohorts and a "missing proxy" scenario are added
# R4_SMOKE=1 -> 3 runs per scenario, 10 references
suppressPackageStartupMessages({ library(survival); library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
SMOKE <- Sys.getenv("R4_SMOKE", "0") == "1"
OUTD <- file.path(B, if (SMOKE) "N11_out/n4_loo_smoke" else "N11_out/n4_loo"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
NSIM <- if (SMOKE) 3L else 200L; NNULL <- if (SMOKE) 10L else 100L
BETAS <- c(0, 0.2, 0.3, 0.4, 0.6); K <- 5L; NPOOL <- 100L; GAMMA_OX <- 0.3
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%s | %5.1f min] %s\n", format(Sys.time(), "%H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
O <- OUTD   # load_cohorts.R writes cohort_meta to O
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
O <- OUTD

## ---- Covariates (as in ring2_proxy.R) ----
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

G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))
G <- G[(arm == "A") | (arm == "Ak")]; G[, sid := paste(arm, src, kind, rep, sep = "|")]
if (SMOKE) { u <- unique(G[, .(sid, arm, kind)]); G <- G[sid %in% c(u[arm == "A" & kind == "real"]$sid[1:8], u[arm == "A" & kind == "null"]$sid[1:8], u[arm == "Ak"]$sid[1:8])] }
PG <- prox_genes$prolif; OG <- prox_genes$oxphos
one <- function(s) { d <- G[sid == s]; b <- setNames(d$coef, d$gene)
  rbindlist(lapply(names(cohorts), function(k) { z <- cohorts[[k]]$z; sc <- drop(t(z[names(b), , drop = FALSE]) %*% b)
    pl <- colMeans(z[setdiff(PG, names(b)), , drop = FALSE]); ox <- colMeans(z[setdiff(OG, names(b)), , drop = FALSE])
    dd <- data.frame(s = sc, prolif = pl, oxphos = ox, stromal = cohorts[[k]]$stromal, immune = cohorts[[k]]$immune)
    df <- data.frame(s = sc, prolif = cohorts[[k]]$prolif, oxphos = cohorts[[k]]$oxphos, stromal = cohorts[[k]]$stromal, immune = cohorts[[k]]$immune)
    data.table(sid = s, arm = d$arm[1], kind = d$kind[1], cohort = k, n_prolif_removed = sum(PG %in% names(b)),
               rho_prolif_loo = cor(sc, pl, method = "spearman"), R2_proxy_loo = summary(lm(s ~ prolif + stromal + immune + oxphos, data = dd))$r.squared,
               R2_proxy_full = summary(lm(s ~ prolif + stromal + immune + oxphos, data = df))$r.squared) })) }
cl <- makeCluster(as.integer(Sys.getenv("NCORE", "2")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table))))
clusterExport(cl, c("G", "PG", "OG", "cohorts", "one"))
R <- rbindlist(parLapply(cl, unique(G$sid), one)); stopCluster(cl)
R[, grp := fcase(arm == "A" & kind == "real", "Published", arm == "Ak", "Size-matched random", arm == "A" & kind == "null", "Input-matched random")]
fwrite(R, file.path(O, "loo_rows.tsv"), sep = "\t")
r2 <- fread(file.path(B, "N11_out/ring2/ring2_rows.tsv"))[arm %in% c("A", "Ak")][, sid := paste(arm, src, kind, rep, sep = "|")]
chk <- merge(R[, .(sid, cohort, R2_proxy_full)], r2[, .(sid, cohort, R2_proxy)], by = c("sid", "cohort"))
fwrite(data.table(n_matched = nrow(chk), n_rows = nrow(R), max_abs_diff = max(abs(chk$R2_proxy_full - chk$R2_proxy))), file.path(O, "loo_consistency_check.tsv"), sep = "	")
S <- rbindlist(lapply(c("TCGA", "GSE14520", "LIRI"), function(k) { d <- R[cohort == k]
  rbindlist(lapply(c("R2_proxy_loo", "rho_prolif_loo"), function(m) { med <- d[, .(m = median(get(m), na.rm = TRUE)), by = grp]
    data.table(cohort = k, metric = m, pub = med[grp == "Published"]$m, sizematched = med[grp == "Size-matched random"]$m, inputmatched = med[grp == "Input-matched random"]$m,
               P_pub_vs_size = wilcox.test(d[grp == "Published"][[m]], d[grp == "Size-matched random"][[m]])$p.value,
               P_pub_vs_input = wilcox.test(d[grp == "Published"][[m]], d[grp == "Input-matched random"][[m]])$p.value) })) }))
S[, order_holds := sizematched > pub & pub > inputmatched][, robust := order_holds & P_pub_vs_size < 0.05]
fwrite(S, file.path(O, "loo_summary.tsv"), sep = "\t")
fwrite(R[cohort == "TCGA", .(median_removed = median(n_prolif_removed), mean_removed = mean(n_prolif_removed), max_removed = max(n_prolif_removed)), by = grp], file.path(O, "loo_overlap.tsv"), sep = "\t")
print(S); print(fread(file.path(O, "loo_overlap.tsv"))); tm("DONE")
