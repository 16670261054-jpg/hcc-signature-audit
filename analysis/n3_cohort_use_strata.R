# External validation stratified by whether the source article used an external cohort (descriptive only)
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out"); OD <- file.path(O, "n3_cohort_use"); dir.create(OD, showWarnings = FALSE)
I <- fread(file.path(B, "N12_signatures/signatures_index.tsv"), encoding = "UTF-8")
body <- function(pm) { f <- file.path(B, "N12_signatures/raw/txt", paste0(pm, ".txt")); if (!file.exists(f)) return(NA_character_)
  L <- readLines(f, warn = FALSE, encoding = "UTF-8"); cut <- grep("^\\S*\\s*references\\s*$", L, ignore.case = TRUE)[1]
  if (!is.na(cut)) L <- L[seq_len(cut - 1)]; paste(L, collapse = " ") }
I[, txt := vapply(as.character(PMID), body, "")]
I[, `:=`(used_GSE = grepl("GSE14520", txt), used_LIRI = grepl("LIRI|ICGC", txt))]
LS <- fread(file.path(B, "N28_fig/source_data/Fig4c_landscape_per_signature_modules.tsv"))
P <- merge(LS[, .(src, M0_pos)], I[, .(src = sig_id, used_GSE, used_LIRI)], by = "src"); stopifnot(nrow(P) == 148, !anyNA(P$used_GSE))
N <- rbindlist(lapply(list.files(file.path(O, "ring3_armA_fitsig_pool100_k/chunks"), "^m06_", full.names = TRUE), fread))
N[, M0_pos := as.integer(M0_p_GSE < 0.05 & M0_p_LIRI < 0.05 & M0_c_GSE > 0.5 & M0_c_LIRI > 0.5)]
N <- N[arm == "A" & kind == "null", .(nn = .N, xn = sum(M0_pos)), by = src]
D <- merge(P, N, by = "src"); stopifnot(nrow(D) == 148)
cmh <- function(d) { a <- array(0, c(2, 2, nrow(d))); for (i in seq_len(nrow(d))) a[, , i] <- matrix(c(d$M0_pos[i], 1 - d$M0_pos[i], d$xn[i], d$nn[i] - d$xn[i]), 2)
  k <- apply(a, 3, function(x) all(rowSums(x) > 0) & all(colSums(x) > 0)); if (sum(k) < 2) return(NA_real_); mantelhaen.test(a[, , k], exact = FALSE)$p.value }
D[, used_any := used_GSE | used_LIRI]
out <- rbindlist(lapply(c("used_GSE", "used_LIRI", "used_any"), function(v) rbindlist(lapply(c(TRUE, FALSE), function(l) { d <- D[get(v) == l]
  data.table(stratifier = v, level = l, n_pub = nrow(d), pub_M0_rate = mean(d$M0_pos), random_M0_rate = sum(d$xn) / sum(d$nn), CMH_P_vs_random = cmh(d)) }))))
btw <- rbindlist(lapply(c("used_GSE", "used_LIRI", "used_any"), function(v) data.table(stratifier = v,
  fisher_P_between_levels = fisher.test(table(D[[v]], D$M0_pos))$p.value)))
fwrite(out, file.path(OD, "m0_by_cohort_use.tsv"), sep = "\t"); fwrite(btw, file.path(OD, "m0_by_cohort_use_between.tsv"), sep = "\t")
fwrite(D[, .(src, used_GSE, used_LIRI, M0_pos)], file.path(OD, "per_signature_cohort_use.tsv"), sep = "\t")
print(out); print(btw)
