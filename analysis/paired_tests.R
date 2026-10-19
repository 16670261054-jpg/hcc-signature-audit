# Published vs matched random signatures: between-group tests that respect the matching structure (decision rules fixed before running)
#   Binary outcomes: CMH stratified by src (mantelhaen.test, exact = FALSE, informative strata only; same construction as analysis/cmh_modules.R)
#   Continuous outcomes: each published signature vs the median of its matched random signatures, paired Wilcoxon signed-rank (wilcox.test(paired = TRUE))
#   No additional correction; descriptive statistics unchanged (the script first asserts that recomputed descriptive values equal the earlier source data)
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); O <- file.path(B, "N11_out"); N35 <- file.path(B, "N35_paired_tests")
OLD <- file.path(N35, "old_source_data")
rd_old <- function(f) fread(file.path(OLD, f), sep = "\t", quote = "")
stars <- function(p) fcase(is.na(p), "", p < 1e-4, "****", p < 1e-3, "***", p < 1e-2, "**", p < 0.05, "*", default = "ns")
T_CMH <- "Cochran-Mantel-Haenszel test stratified by matched set (published signature + its matched random signatures)"
T_WSR <- "Paired Wilcoxon signed-rank test (each published signature vs median of its matched random signatures)"
near <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) < tol))

## Matching check: every random signature's src belongs to the published set, and every published src has random signatures
chk_match <- function(pub_src, nul_src, what) {
  ok <- all(nul_src %in% pub_src) && all(pub_src %in% nul_src)
  if (!ok) stop("not one-to-one matched: ", what)
  invisible(TRUE) }

## CMH (array construction and informative-stratum rule as in cmh_modules.R)
cmh <- function(P, N, what) {   # P: data.table(src, x 0/1); N: data.table(src, x 0/1)
  P <- P[!is.na(x)]; N <- N[!is.na(x)]
  chk_match(unique(P$src), unique(N$src), what)
  stopifnot(!anyDuplicated(P$src))
  Ns <- N[, .(nn = .N, xn = sum(x)), by = src]
  d <- merge(P[, .(src, xr = x)], Ns, by = "src")
  a <- array(0, c(2, 2, nrow(d))); for (i in seq_len(nrow(d))) a[, , i] <- matrix(c(d$xr[i], 1 - d$xr[i], d$xn[i], d$nn[i] - d$xn[i]), 2)
  inf <- apply(a, 3, function(z) all(rowSums(z) > 0) & all(colSums(z) > 0))
  ct <- if (sum(inf) >= 1) tryCatch(mantelhaen.test(a[, , inf, drop = FALSE], exact = FALSE), error = function(e) NULL) else NULL
  list(P = if (is.null(ct)) 1 else ct$p.value, OR = if (is.null(ct)) NA_real_ else unname(ct$estimate),
       n_pub = nrow(d), n_null = sum(d$nn), n_unit = sum(inf), unit = "informative strata",
       rate_pub = mean(d$xr), rate_null = sum(d$xn) / sum(d$nn)) }

## Paired Wilcoxon
wsr <- function(P, N, what) {   # P: data.table(src, v); N: data.table(src, v)
  P <- P[!is.na(v)]; N <- N[!is.na(v)]
  chk_match(unique(P$src), unique(N$src), what)
  stopifnot(!anyDuplicated(P$src))
  d <- merge(P[, .(src, vp = v)], N[, .(vn = median(v)), by = src], by = "src")
  w <- wilcox.test(d$vp, d$vn, paired = TRUE)
  list(P = w$p.value, OR = NA_real_, n_pub = nrow(d), n_null = nrow(N), n_unit = nrow(d), unit = "pairs",
       median_pub = median(P$v), median_null = median(N$v), median_paired_diff = median(d$vp - d$vn)) }

res <- list()
put <- function(key, panel, comparison, old_test, old_P, r, new_test) {
  res[[length(res) + 1]] <<- data.table(key, panel, comparison, old_test, old_P, new_test, new_P = r$P,
    old_star = stars(old_P), new_star = stars(r$P), flip_0.05 = (old_P < 0.05) != (r$P < 0.05),
    n_pub = r$n_pub, n_null = r$n_null, n_strata_or_pairs = r$n_unit, unit = r$unit, CMH_OR = r$OR) }

## ================= Fig 3a-c (proxy analysis, ring2) =================
R <- fread(file.path(O, "ring2/ring2_rows.tsv"))
rA <- R[arm == "A" & kind == "real"]
NULLS <- list("Size-matched random" = R[arm == "Ak"], "Input-matched random" = R[arm == "A" & kind == "null"])
ST_old <- rd_old("Fig3ab_stats.tsv")
for (g in names(NULLS)) for (mt in c("R2_proxy", "rho_prolif")) for (cn in c("TCGA", "GSE14520", "LIRI")) {
  r <- wsr(rA[cohort == cn, .(src, v = get(mt))], NULLS[[g]][cohort == cn, .(src, v = get(mt))], paste("Fig3ab", mt, cn, g))
  o <- ST_old[metric == mt & cohort == cn & grp2 == g]; stopifnot(nrow(o) == 1)
  stopifnot(near(o$median_pub, r$median_pub), near(o$median_null, r$median_null))   # descriptive values equal those of the earlier figure
  put(paste("Fig3ab", mt, cn, g, sep = "|"), ifelse(mt == "R2_proxy", "Fig3a", "Fig3b"),
      sprintf("%s, %s: published vs %s", mt, cn, tolower(g)), "Mann-Whitney U (independent samples)", o$P, r, T_WSR) }
both <- function(d, col) { w <- dcast(d[cohort != "TCGA"], arm + src + kind + rep ~ cohort, value.var = col)
  w[, pos := as.integer(GSE14520 < 0.05 & LIRI < 0.05)]; w }
C3_old <- rd_old("Fig3c_residual.tsv")
for (g in names(NULLS)) for (col in c("lrt_p", "lrt_p_oxphos")) {
  a <- both(rA, col); b <- both(NULLS[[g]], col)
  r <- cmh(a[, .(src, x = pos)], b[, .(src, x = pos)], paste("Fig3c", col, g))
  ml <- if (col == "lrt_p") "Base" else "+ OXPHOS"
  o <- C3_old[grp == g & grepl(if (col == "lrt_p") "^Base: " else "OXPHOS$", model)]; stopifnot(nrow(o) == 1)
  op <- C3_old[grp == "Published" & model == o$model]
  stopifnot(o$n == r$n_null, near(o$rate, r$rate_null), op$n == r$n_pub, near(op$rate, r$rate_pub))
  put(paste("Fig3c", col, g, sep = "|"), "Fig3c", sprintf("Residual signal (%s model): published vs %s", ml, tolower(g)),
      "Fisher exact (independent samples)", o$fisher_P_vs_published, r, T_CMH) }

## ================= Fig 4d DSigDB (hiend_m12_annot) =================
AR <- fread(file.path(O, "hiend_m12_annot/m12_annot_rows.tsv"))[arm == "A"]
AR[, src := sub("^A[|]", "", sid)]
D4_old <- rd_old("Fig4d_dsigdb_rates.tsv")
D4_old[, code := c("real", "N_plain", "N_annot")[match(gsub("[ ]+", " ", grp), c("Published", "Random", "Annotation- matched"))]]
stopifnot(!anyNA(D4_old$code))
for (nl in c("N_plain", "N_annot")) {
  r <- cmh(AR[type == "real", .(src, x = pos)], AR[type == nl, .(src, x = pos)], paste("Fig4d", nl))
  o <- D4_old[code == nl]; op <- D4_old[code == "real"]
  stopifnot(o$n == r$n_null, o$pos == round(r$rate_null * r$n_null), op$n == r$n_pub, near(op$rate, r$rate_pub))
  put(paste("Fig4d", nl, sep = "|"), "Fig4d", paste("DSigDB positive: published vs", c(N_plain = "random", N_annot = "annotation-matched random")[[nl]]),
      "Fisher exact (independent samples)", o$P, r, T_CMH) }

## ================= Fig 6f M8 input-matched (hiend_m8; the size-matched comparison is already CMH + Holm and is listed unchanged) =================
M8 <- fread(file.path(O, "hiend_m8/m8_per_signature.tsv"))
F6_old <- rd_old("Fig6f_M8_rates.tsv")
CMHmod <- fread(file.path(O, "cmh_modules.tsv"))
for (m in c("M8a_pos", "M8b_pos")) {
  r <- cmh(M8[arm == "A" & kind == "real", .(src, x = get(m))], M8[arm == "A" & kind == "null", .(src, x = get(m))], paste("Fig6f", m))
  o <- F6_old[module == m & grp == "Input-matched random"]; op <- F6_old[module == m & grp == "Published"]
  stopifnot(o$n == r$n_null, near(o$rate, r$rate_null), near(op$rate, r$rate_pub))
  put(paste("Fig6f", m, "Input-matched random", sep = "|"), "Fig6f", paste(sub("_pos", "", m), "positive: published vs input-matched random"),
      "Fisher exact (independent samples)", o$P, r, T_CMH)
  os <- F6_old[module == m & grp == "Size-matched random"]; cm <- CMHmod[module == sub("_pos", "", m)]
  stopifnot(near(os$P, cm$CMH_P_holm))
  put(paste("Fig6f", m, "Size-matched random", sep = "|"), "Fig6f", paste(sub("_pos", "", m), "positive: published vs size-matched random (unchanged)"),
      "CMH stratified by matched set, Holm across 16 modules (Fig 4a family)", os$P,
      list(P = cm$CMH_P_holm, OR = cm$CMH_OR, n_pub = cm$n_pub, n_null = cm$n_rand, n_unit = cm$n_informative_strata, unit = "informative strata"),
      "CMH stratified by matched set, Holm across 16 modules (Fig 4a family; unchanged)") }

## ================= Fig 6g right, Geneformer abundance-matched (hiend_m10) =================
G <- fread(file.path(O, "hiend_m10/m10_abund_matched_rows.tsv"))
r <- cmh(G[type == "real", .(src, x = pos)], G[type == "abund_matched", .(src, x = pos)], "Fig6g")
G_old <- rd_old("Fig6g_geneformer_abundance_matched.tsv")
stopifnot(G_old[grp != "Published"]$n == r$n_null, near(G_old[grp != "Published"]$rate, r$rate_null), near(G_old[grp == "Published"]$rate, r$rate_pub))
put("Fig6g|abund_matched", "Fig6g", "Geneformer positive: published vs abundance-matched random", "Fisher exact (independent samples)", G_old[grp != "Published"]$P, r, T_CMH)

## ================= Fig S2 (hiend_m14_prolif) =================
P14 <- fread(file.path(O, "hiend_m14_prolif/m14_prolif_rows.tsv"))
P14[, src := vapply(strsplit(sid, "|", fixed = TRUE), `[`, "", 2)]
V <- P14[!is.na(prolif_sep)]; a <- V[arm == "A"]; k <- V[arm == "Ak"]   # same analysis set as hiend_m14_prolif_check.R
D14_old <- rd_old("FigS2b_m14_prolif_decision.tsv")
r <- wsr(a[, .(src, v = prolif_sep)], k[, .(src, v = prolif_sep)], "FigS2a")
stopifnot(near(D14_old$prolif_sep_median_pub, r$median_pub), near(D14_old$prolif_sep_median_Ak, r$median_null), D14_old$n_Ak == nrow(k), D14_old$n_pub == nrow(a))
put("FigS2a|sep", "FigS2a", "Subtype separation by proliferation: published vs size-matched random", "Mann-Whitney U (independent samples)", D14_old$MWU_P, r, T_WSR)
for (st in c("raw", "adj")) { col <- c(raw = "M14_p", adj = "adj_p")[[st]]
  r <- cmh(a[, .(src, x = as.integer(get(col) < 0.05))], k[, .(src, x = as.integer(get(col) < 0.05))], paste("FigS2b", st))
  stopifnot(near(D14_old[[paste0(st, "_rate_pub")]], r$rate_pub), near(D14_old[[paste0(st, "_rate_Ak")]], r$rate_null))
  put(paste("FigS2b", st, sep = "|"), "FigS2b", paste0("Consensus subtypes positive (", c(raw = "unadjusted", adj = "proliferation-adjusted")[[st]], "): published vs size-matched random"),
      "Fisher exact (independent samples)", D14_old[[paste0(st, "_fisher_P")]], r, T_CMH) }

RES <- rbindlist(res)
f <- file.path(N35, "paired_tests_results.tsv")
fwrite(RES, f, sep = "\t", quote = FALSE)
back <- fread(f, sep = "\t", quote = "")
stopifnot(nrow(back) == nrow(RES), nrow(RES) == 26)
print(RES[, .(panel, comparison, old_P = signif(old_P, 3), new_P = signif(new_P, 3), old_star, new_star, flip_0.05, n_pub, n_strata_or_pairs)])
cat("PAIRED TESTS OK:", nrow(RES), "rows; flips =", sum(RES$flip_0.05), "\n")
