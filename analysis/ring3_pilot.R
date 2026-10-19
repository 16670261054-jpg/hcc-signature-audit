# 第 3 环可行性试跑：每个模块都在 1 个真实集 + 3 个随机集上跑通，记录耗时，检查输出不是 NA；不做任何判读
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table) })
set.seed(20261017)
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); TG <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "TCGA")
O <- file.path(B, "N11_out/ring3_pilot"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
source(file.path(B, "analysis/load_cohorts.R"), encoding = "UTF-8")
source(file.path(B, "analysis/ring3_lib.R"), encoding = "UTF-8")
tm <- list()
t0 <- Sys.time(); P <- prep_ring3(B, TG); tm$prep <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("prep secs", tm$prep, " imm dims", dim(P$imm), " tmb nonNA", sum(!is.na(P$tmb)), " H", length(P$H), "\n")
fwrite(data.table(score = rownames(P$imm), n_NA = rowSums(is.na(P$imm)), sd = apply(P$imm, 1, sd)), file.path(O, "prep_immune_check.tsv"), sep = "\t")

g_real <- sets$cuproptosis_Tsvetkov19; k <- length(g_real)
glist <- c(list(real = g_real), setNames(lapply(1:3, function(i) sample(U, k)), paste0("rand", 1:3)))
rows <- list(); mt <- list()
for (nm in names(glist)) {
  b <- fitsig(glist[[nm]], cohorts$TCGA$z, cohorts$TCGA$y); r <- c()
  for (m in c("m0", "m1", "m2", "m4", "m5", "m6")) {
    t0 <- Sys.time()
    v <- switch(m, m0 = m0_external(b), m1 = m1_indep(b, P), m2 = m2_immune(b, P), m4 = m4_gsea(b, P),
                m5 = m5_tmb(b, P), m6 = m6_nomo(b))
    mt[[length(mt) + 1]] <- data.table(set = nm, module = m, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    r <- c(r, v)
  }
  rows[[nm]] <- data.table(set = nm, n_model = length(b), t(r))
}
res <- rbindlist(rows, fill = TRUE); fwrite(res, file.path(O, "pilot_modules.tsv"), sep = "\t"); print(res)
m7t <- list()
for (nm in c("real", "rand1")) { t0 <- Sys.time(); o <- m7_multi(glist[[nm]])
  m7t[[nm]] <- data.table(set = nm, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), t(o$summary), best = o$best)
  fwrite(o$table, file.path(O, paste0("pilot_m7_", nm, ".tsv")), sep = "\t"); print(o$table) }
m7t <- rbindlist(m7t); fwrite(m7t, file.path(O, "pilot_m7_summary.tsv"), sep = "\t"); print(m7t)
mt <- rbindlist(mt); fwrite(mt, file.path(O, "pilot_timing.tsv"), sep = "\t")
print(mt[, .(mean_secs = mean(secs)), by = module])
na_cols <- names(res)[vapply(res, function(x) any(is.na(x)), logical(1))]
cat("NA_COLUMNS:", if (length(na_cols)) paste(na_cols, collapse = ",") else "none", "\nDONE\n")
