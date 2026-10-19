# Shared: loads the three cohorts + gene sets + proxy scores + template pipeline functions (sourced by the analysis scripts)
GE <- file.path(Sys.getenv("HCC_AUDIT_DATA", "data/raw"), "GEO_bulk")

stage2 <- function(x) ifelse(grepl("IV|III", x), "III-IV", ifelse(grepl("^I$|^II$|Stage I[A-C]?$|Stage II[A-C]?$|^I[A-C]?$|^II[A-C]?$", x), "I-II", NA))

## ---- TCGA-LIHC ----
ex <- fread(file.path(TG, "TCGA-LIHC_gene_expression.tsv.gz"), data.table = FALSE)
hg <- fread(file.path(B, "N12_signatures/ref/hgnc.txt"), select = c("symbol", "entrez_id", "status", "prev_symbol"), data.table = FALSE, quote = "")
hg <- hg[hg$status == "Approved", ]
sp <- strsplit(hg$prev_symbol, "|", fixed = TRUE); pv <- data.frame(old = trimws(unlist(sp)), new = rep(hg$symbol, lengths(sp))); pv <- pv[nzchar(pv$old), ]
pv <- pv[!duplicated(pv$old) & !(pv$old %in% pv$old[duplicated(pv$old)]) & !(pv$old %in% hg$symbol), ]
hg <- hg[!is.na(hg$entrez_id), ]
upd <- function(x) { m <- match(x, pv$old); ifelse(is.na(m), x, pv$new[m]) }
ex[[1]] <- upd(ex[[1]]); ex <- ex[!duplicated(ex[[1]]), ]
rownames(ex) <- ex[[1]]; ex <- as.matrix(ex[, -1]); ex <- ex[, substr(colnames(ex), 14, 15) == "01"]
sv <- fread(file.path(TG, "TCGA-LIHC_survival.tsv"), data.table = FALSE)
ph <- fread(file.path(TG, "TCGA-LIHC_phenotype.tsv"), data.table = FALSE)
sv <- sv[sv$sample %in% colnames(ex) & !is.na(sv$OS.time) & sv$OS.time > 0 & !is.na(sv$OS), ]
stcol <- grep("pathologic_stage", colnames(ph), value = TRUE)[1]
tc <- list(x = ex[, sv$sample], time = sv$OS.time / 30.44, event = sv$OS,
           stage = stage2(gsub("Stage ", "", ph[[stcol]][match(sv$sample, ph[[1]])])))

## ---- GSE14520 (GPL3921) ----
gzr <- function(f) readLines(gzfile(f))
r0 <- gzr(file.path(GE, "GSE14520/GSE14520-GPL3921_series_matrix.txt.gz")); r0 <- r0[!startsWith(r0, "!") & nzchar(r0)]
mx <- fread(text = r0, data.table = FALSE)
pr <- mx[[1]]; mx <- as.matrix(mx[, -1]); rownames(mx) <- pr
a0 <- gzr(file.path(GE, "GSE14520/GPL3921.annot.gz")); a0 <- a0[grep("^ID\t", a0)[1]:length(a0)]; a0 <- a0[!startsWith(a0, "!")]
an <- fread(text = a0, data.table = FALSE, sep = "\t", quote = "", fill = TRUE)
e2s <- function(e) hg$symbol[match(suppressWarnings(as.integer(e)), hg$entrez_id)]
gid <- an[["Gene ID"]][match(rownames(mx), an$ID)]
sym <- ifelse(grepl("///", gid), NA, e2s(gid))
ok <- !is.na(sym) & sym != ""
mx <- mx[ok, ]; sym <- sym[ok]
if (max(mx, na.rm = TRUE) > 50) mx <- log2(mx + 1)
mm <- rowMeans(mx); o <- order(sym, -mm); keep <- o[!duplicated(sym[o])]
mx <- mx[keep, ]; rownames(mx) <- sym[keep]
su <- fread(text = gzr(file.path(GE, "GSE14520/GSE14520_Extra_Supplement.txt.gz")), data.table = FALSE)
su <- su[su[["Tissue Type"]] == "Tumor" & su$Affy_GSM %in% colnames(mx) & !is.na(su[["Survival months"]]) & su[["Survival months"]] > 0, ]
g1 <- list(x = mx[, su$Affy_GSM], time = su[["Survival months"]], event = as.integer(su[["Survival status"]]),
           stage = ifelse(su[["TNM staging"]] %in% c("III", "IIIA", "IIIB", "IIIC", "IV"), "III-IV",
                          ifelse(su[["TNM staging"]] %in% c("I", "II"), "I-II", NA)))

## ---- LIRI-JP (HCCDB18) ----
L <- file.path(GE, "HCCDB18_LIRI-JP")
lx <- fread(file.path(L, "HCCDB18_mRNA_level3.txt"), data.table = FALSE)
lx$Symbol <- e2s(lx$Entrez_ID); lx <- lx[!is.na(lx$Symbol) & lx$Symbol != "", ]
lm <- as.matrix(lx[, -(1:2)]); lmm <- rowMeans(lm); o <- order(lx$Symbol, -lmm); keep <- o[!duplicated(lx$Symbol[o])]
lm <- lm[keep, ]; rownames(lm) <- lx$Symbol[keep]
rd <- function(f) { r <- readLines(f); s <- strsplit(r, "\t"); setNames(lapply(s, `[`, -1), vapply(s, `[`, "", 1)) }
ls_ <- rd(file.path(L, "HCCDB18.sample.txt")); lp <- rd(file.path(L, "HCCDB18.patient.txt"))
smp <- data.frame(sid = ls_$SAMPLE_ID, type = ls_$TYPE, pid = ls_$PATIENT_ID)
smp <- smp[smp$type == "HCC" & smp$sid %in% colnames(lm), ]; smp <- smp[!duplicated(smp$pid), ]
pi <- match(smp$pid, lp$PATIENT_ID)
lt <- as.numeric(lp$SUR[pi]); le <- as.integer(lp$STATUS[pi] == "Dead")
okl <- !is.na(lt) & lt > 0 & !is.na(le)
smp <- smp[okl, ]; pi <- pi[okl]
l1 <- list(x = lm[, smp$sid], time = lt[okl], event = le[okl], stage = stage2(lp$TNM_STAGE_T1[pi]))

## ---- Common gene universe + expression filter + within-cohort z ----
cohorts <- list(TCGA = tc, GSE14520 = g1, LIRI = l1)
Uall <- Reduce(intersect, lapply(cohorts, function(c) rownames(c$x)))   # for scoring published signatures: measured in all three cohorts
U <- Uall[rowMeans(tc$x[Uall, ] > 1) > 0.5]                          # random sampling pool: TCGA expression filter
for (k in names(cohorts)) {
  x <- cohorts[[k]]$x[Uall, ]; x <- x[, colSums(is.na(x)) == 0]
  sdv <- apply(x, 1, sd); x[sdv == 0, ] <- 0
  cohorts[[k]]$z <- t(scale(t(x))); cohorts[[k]]$z[!is.finite(cohorts[[k]]$z)] <- 0
  cohorts[[k]]$y <- Surv(cohorts[[k]]$time, cohorts[[k]]$event)
}
meta <- rbindlist(lapply(names(cohorts), function(k) data.table(cohort = k, n = ncol(cohorts[[k]]$z),
  events = sum(cohorts[[k]]$event), stage_NA = sum(is.na(cohorts[[k]]$stage)), n_genes_universe = length(U), n_genes_all = length(Uall))))
fwrite(meta, file.path(O, "cohort_meta.tsv"), sep = "\t"); print(meta)

## ---- Gene sets ----
rg <- function(f) { s <- strsplit(readLines(file.path(B, "N11_geneset", f)), "\t"); setNames(lapply(s, function(x) x[-(1:2)]), vapply(s, `[`, "", 1)) }
H <- rg("h.all.v2024.1.Hs.symbols.gmt")
C2 <- c(rg("c2.cp.reactome.v2024.1.Hs.symbols.gmt"), rg("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt"))
C5 <- rg("c5.go.bp.v2024.1.Hs.symbols.gmt")
sets <- list(
  cuproptosis_Tsvetkov19 = c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1","GLS","CDKN2A",
                             "NFE2L2","NLRP3","ATP7B","ATP7A","SLC31A1","LIPT2","DBT","GCSH","DLST"),
  disulfidptosis_10 = c("GYS1","NDUFS1","OXSM","LRPPRC","NDUFA11","NUBPL","NCKAP1","RPN1","SLC3A2","SLC7A11"),
  WP_FERROPTOSIS = C2$WP_FERROPTOSIS, GOBP_FERROPTOSIS = C5$GOBP_FERROPTOSIS,
  REACTOME_PYROPTOSIS = C2$REACTOME_PYROPTOSIS, GOBP_NECROPTOTIC_PROCESS = C5$GOBP_NECROPTOTIC_PROCESS,
  REACTOME_REGULATED_NECROSIS = C2$REACTOME_REGULATED_NECROSIS, GOBP_ANOIKIS = C5$GOBP_ANOIKIS)
sets <- lapply(sets, function(g) intersect(g, U))
prox_genes <- list(prolif = intersect(union(H$HALLMARK_E2F_TARGETS, H$HALLMARK_G2M_CHECKPOINT), U),
                   oxphos = intersect(H$HALLMARK_OXIDATIVE_PHOSPHORYLATION, U))
for (k in names(cohorts)) for (p in names(prox_genes)) cohorts[[k]][[p]] <- colMeans(cohorts[[k]]$z[prox_genes[[p]], ])

## ---- Template pipeline: univariable Cox P<0.05 -> LASSO-Cox -> coefficients; if none selected, mean of the full set signed by the training direction ----
fitsig <- function(g, Z, y) {
  X <- t(Z[g, , drop = FALSE])
  p <- apply(X, 2, function(v) { s <- tryCatch(summary(coxph(y ~ v))$coefficients[5], error = function(e) 1); s })
  kp <- names(p)[p < 0.05]
  if (length(kp) >= 2) {
    f <- tryCatch(suppressWarnings(cv.glmnet(X[, kp], y, family = "cox", nfolds = 5)), error = function(e) NULL)
    if (!is.null(f)) { b <- as.matrix(coef(f, s = "lambda.min"))[, 1]; b <- b[b != 0]; if (length(b)) return(b) }
  }
  if (length(kp) == 1) return(setNames(sign(coef(coxph(y ~ X[, kp]))), kp))
  sg <- sign(coef(coxph(y ~ rowMeans(X)))); setNames(rep(sg / length(g), length(g)), g)
}
score <- function(b, Z) drop(t(Z[names(b), , drop = FALSE]) %*% b)
cidx <- function(s, y) unname(concordance(y ~ s, reverse = TRUE)$concordance)
