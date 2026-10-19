# M12 DSigDB annotation-matched null (analysis plan fixed before running)
#   Real signatures (A published / B mechanism sets; genes = rows with kind "real" in signature_genes_all_arms.tsv)
#   vs N_annot: size-matched random genes matched stratum by stratum on DSigDB annotation degree (one stratum for 0 + 10 deciles within pool U)
#   vs N_plain: size-matched plain random genes (no annotation matching)
#   enricher parameters and positivity criterion as in M12 (template workflow)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
O <- file.path(B, "N11_out/hiend_m12_annot"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
DL <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "external")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
RA <- 25L; RB <- 200L

U <- readLines(file.path(B, "N11_out/hiend_m8/U_pool.txt"))
dr <- fread(file.path(DL, "DSigDB/DSigDB_All_detailed.txt"), sep = "\t", quote = "", select = 1:2, col.names = c("Drug", "Gene"))
dr <- unique(dr[!is.na(Gene) & Gene != ""])
ann <- dr[, .(n_drug = uniqueN(Drug)), by = Gene]
acnt <- function(g) { v <- ann$n_drug[match(g, ann$Gene)]; v[is.na(v)] <- 0L; v }
uc <- acnt(U)
brk <- unique(quantile(uc[uc > 0], 0:10 / 10))
layer <- function(v) ifelse(v == 0, 0L, as.integer(cut(v, brk, include.lowest = TRUE, labels = FALSE)))
# Annotation degrees outside the range of pool U go to the end strata
layer2 <- function(v) { l <- layer(pmin(pmax(v, min(brk)), max(brk))); l[v == 0] <- 0L; l }
uL <- layer2(uc); pool <- split(U, uL)
fwrite(data.table(layer = names(pool), n_genes = lengths(pool), min_ann = sapply(split(uc, uL), min), max_ann = sapply(split(uc, uL), max)),
       file.path(O, "annotation_layers.tsv"), sep = "\t")

G <- fread(file.path(B, "N11_out/hiend_m8/signature_genes_all_arms.tsv"))
real <- G[kind == "real" & arm %in% c("A", "B")]
real_sets <- split(real$gene, paste(real$arm, real$src, sep = "|"))
draw_annot <- function(g, seed) { set.seed(seed); L <- layer2(acnt(g)); tb <- table(L)
  unlist(lapply(names(tb), function(l) sample(pool[[l]], min(tb[[l]], length(pool[[l]]))))) }
draw_plain <- function(g, seed) { set.seed(seed); sample(U, length(g)) }

jobs <- rbindlist(lapply(seq_along(real_sets), function(i) { nm <- names(real_sets)[i]; arm <- sub("\\|.*", "", nm)
  R <- if (arm == "A") RA else RB
  rbind(data.table(sid = nm, arm = arm, type = "real", rep = 0L, idx = i),
        data.table(sid = nm, arm = arm, type = "N_annot", rep = seq_len(R), idx = i),
        data.table(sid = nm, arm = arm, type = "N_plain", rep = seq_len(R), idx = i)) }))
tm(sprintf("jobs %d", nrow(jobs)))

one <- function(i) { j <- jobs[i]; g0 <- real_sets[[j$sid]]
  g <- switch(j$type, real = g0, N_annot = draw_annot(g0, 8e6 + 1e3 * j$idx + j$rep), N_plain = draw_plain(g0, 9e6 + 1e3 * j$idx + j$rep))
  kk <- suppressMessages(clusterProfiler::enricher(unique(g), pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, TERM2GENE = dr))
  r <- if (is.null(kk)) data.frame() else as.data.frame(kk)
  n <- if (nrow(r)) sum(r$pvalue < 0.05 & r$p.adjust < 0.05) else 0L
  data.table(j, k = length(g), mean_ann = mean(acnt(g)), n_sig = n, pos = as.integer(n > 0)) }
cl <- makeCluster(min(16, detectCores() - 2))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(clusterProfiler) })))
clusterExport(cl, c("jobs", "real_sets", "draw_annot", "draw_plain", "layer2", "layer", "brk", "acnt", "ann", "pool", "U", "dr", "one"))
R <- rbindlist(parLapply(cl, seq_len(nrow(jobs)), one)); stopCluster(cl)
fwrite(R, file.path(O, "m12_annot_rows.tsv"), sep = "\t"); tm("enrichment done")

full_ids <- fread(file.path(B, "N11_out/ring3_full/published_coverage.tsv"))[full == TRUE]$sig_id
R[, src := sub("^[^|]*\\|", "", sid)]
cmp <- function(d, lab) rbindlist(lapply(c("N_annot", "N_plain"), function(nt) {
  re <- d[type == "real"]$pos; nu <- d[type == nt]$pos
  xr <- sum(re); nr <- length(re); xn <- sum(nu); nn <- length(nu); ci <- binom.test(xn, nn)$conf.int
  ft <- fisher.test(matrix(c(xr, nr - xr, xn, nn - xn), 2))
  data.table(group = lab, null = nt, n_real = nr, rate_real = xr / nr, n_null = nn, rate_null = xn / nn, null_CP_lo = ci[1], null_CP_hi = ci[2],
             fisher_OR = unname(ft$estimate), fisher_P = ft$p.value,
             mean_ann_real_median = median(d[type == "real"]$mean_ann), mean_ann_null_median = median(d[type == nt]$mean_ann),
             n_sig_real_median = median(d[type == "real"]$n_sig), n_sig_null_median = median(d[type == nt]$n_sig)) }))
D <- rbind(cmp(R[arm == "A"], "A_published"), cmp(R[arm == "B"], "B_mechanism"), cmp(R[arm == "A" & src %in% full_ids], "A_fullcoverage"))
D[, verdict := ifelse(null == "N_annot", ifelse(fisher_P >= 0.05, "explained_by_annotation", "difference_remains"), "reference")]
fwrite(D, file.path(O, "m12_annot_decision.tsv"), sep = "\t")
mw <- wilcox.test(R[type == "real" & arm == "A"]$mean_ann, uc, exact = FALSE)
fwrite(data.table(group = "A_published_vs_Upool", median_real = median(R[type == "real" & arm == "A"]$mean_ann), median_U_gene = median(uc),
                  mean_U_gene = mean(uc), MWU_P = mw$p.value), file.path(O, "m12_annotation_level.tsv"), sep = "\t")
print(D); tm("DONE")
