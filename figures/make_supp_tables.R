# Supplementary Tables S1-S3 (xlsx): reads files on disk only
#  S1 149 signatures (source + per-gene coefficients + source sentences); per-article module calls / gate results go to repo_source_data/ (code repository only)
#  S2 search queries and hit counts, corpus module keyword regexes, criteria and parameters of the 16 modules, 8 cell-death gene sets
#  S3 = TableS3_numbers.tsv (all numbers)
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
B <- Sys.getenv("HCC_AUDIT_ROOT", "."); S12 <- file.path(B, "N12_signatures"); OUT <- file.path(B, "N11_out"); FD <- file.path(B, "N28_fig")
SUP <- file.path(FD, "supplementary"); dir.create(SUP, showWarnings = FALSE)
I <- fread(file.path(S12, "signatures_index.tsv"), encoding = "UTF-8")[group_primary == "是" & 是否有完整系数 == "是"]
stopifnot(nrow(I) == 149)
L <- fread(file.path(S12, "signatures_long.tsv"), encoding = "UTF-8")[sig_id %in% I$sig_id]
## S1a signature level
S1a <- I[, .(signature_id = sig_id, PMID, PMCID, first_author, year, journal, cell_death_type = 机制, outcome, n_genes = 基因数, source_of_coefficients = source)]
LS <- fread(file.path(FD, "source_data/Fig4c_landscape_per_signature_modules.tsv"))
G2 <- fread(file.path(OUT, "n2_gates/published_gates.tsv"))
S1a <- merge(S1a, LS[, !c("mech", "year")], by.x = "signature_id", by.y = "src", all.x = TRUE)
S1a <- merge(S1a, G2[, .(signature_id = src, template_criterion = T0, gate1_pooled = T1p, gate2_pooled = T2p, both_gates = T1p & T2p,
                         C_GSE14520 = C_GSE, C_LIRI = C_LIRI, empirical_P_pooled = empP_pooled, LRT_P_pooled = lrt_pooled)], by = "signature_id", all.x = TRUE)
setnames(S1a, "k", "genes_in_scoring_set", skip_absent = TRUE)
## Per-article module calls and gate results are not in the supplementary tables; they are provided as Source Data in the code repository
REPO <- file.path(FD, "repo_source_data"); dir.create(REPO, showWarnings = FALSE)
src_cols <- c("signature_id", "PMID", "PMCID", "first_author", "year", "journal", "cell_death_type", "outcome", "n_genes", "source_of_coefficients")
## English display labels (display text only, content unchanged)
EN_SRC <- c("正文公式(Results/Methods段)" = "Risk-score formula in main text (Results/Methods)", "正文公式(人工抄录)" = "Risk-score formula in main text (manually transcribed)")
en_src <- function(x) { y <- unname(EN_SRC[x]); ifelse(is.na(y), x, y) }
S1a[, source_of_coefficients := en_src(source_of_coefficients)]
fwrite(S1a, file.path(REPO, "per_signature_module_calls_and_gates.tsv"), sep = "	")
stopifnot(nrow(fread(file.path(REPO, "per_signature_module_calls_and_gates.tsv"))) == 149)
S1a <- S1a[, ..src_cols]
## S1b gene level
S1b <- L[, .(signature_id = sig_id, PMID, gene_as_reported = gene_original, gene_HGNC = gene_hgnc, HGNC_mapping = hgnc_status, coefficient = coef, coefficient_source = coef_source, source_sentence = evidence_quote)]
S1b[, coefficient_source := en_src(coefficient_source)]
## S2
SL <- fread(file.path(S12, "search_log.tsv"), encoding = "UTF-8")
EN_SL <- "Manual review of title-screening misses: 13 articles using spellings such as hepatic cancer, hepatocarcinoma, cuprotosis, RCD, entosis or cell-in-cell were added to the candidates (PMIDs 41611937, 38737687, 37608927, 36110946, 36248857, 36909185, 38799554, 37014321, 38197309, 38662757, 39528670, 39937284, 40794333)"
for (j in names(SL)) if (is.character(SL[[j]])) set(SL, which(startsWith(SL[[j]], "人工复核标题过滤漏检")), j, EN_SL)
MF <- fread(file.path(S12, "module_frequency_fulltext.tsv"))[, .(analysis = module, articles = n_papers, of = n_fulltext, regex)]
MOD <- fread(file.path(FD, "source_data/Fig4a_module_rates.tsv"))[, .(module, label, category)]
CRIT <- data.table(module = c("M0", "M1", "M2", "M3", "M4", "M5a", "M5b", "M6", "M7", "M8a", "M8b", "M10", "M12", "M13", "M14", "M15"),
  criterion = c("C > 0.5 and two-sided Noether P < 0.05 (survcomp) in both GSE14520 and LIRI-JP",
                "Score P < 0.05 in Cox(score + age + sex + stage), TCGA-LIHC",
                "Any of 29 ssGSEA immune sets or 3 ESTIMATE scores: BH-adjusted Wilcoxon P < 0.05, high vs low risk",
                "Any of 198 GDSC2 compounds (oncoPredict, batchCorrect = eb, powerTransformPhenotype = TRUE): BH-adjusted Wilcoxon P < 0.05",
                "Any Hallmark set: fgsea padj < 0.05, genes ranked by limma t (high vs low), minSize 15, maxSize 500",
                "TMB (MC3 non-synonymous / 38 Mb) Wilcoxon P < 0.05, high vs low risk",
                "Log-rank P < 0.05 across four TMB-median x risk-median groups",
                "Nested Cox LRT P < 0.05 for adding score to stage, in both external cohorts",
                "Best of glmnet (alpha 0, 0.5, 1), stepwise Cox, RSF, GBM, superpc by mean external C; Noether P < 0.05 in both external cohorts",
                "Any signature gene among FindAllMarkers markers (|log2FC| > 1, adjusted P < 0.05), GSE149614",
                "Any signature gene differentially expressed tumour vs normal within any cell type (FindMarkers, |log2FC| > 1, adjusted P < 0.05, >= 5 cells per group)",
                "Any signature gene with Geneformer V1-10M in silico deletion goal-state-shift FDR < 0.05 (1,000 tumour hepatocytes)",
                "Any DSigDB drug set with adjusted P < 0.05 (clusterProfiler enricher, minGSSize 10, maxGSSize 500)",
                "Any of 47 immune checkpoint genes Wilcoxon P < 0.05, high vs low risk",
                "Log-rank P < 0.05 between k = 2 ConsensusClusterPlus subtypes (maxK 9, reps 50, pItem 0.8, k-means, Euclidean)",
                "Nomogram (risk group + age + sex + stage) net benefit above all alternatives at >= 1 threshold, 36 months (dcurves)"))
S2c <- merge(MOD, CRIT, by = "module"); stopifnot(nrow(S2c) == 16)
gs <- c("cuproptosis_Tsvetkov19", "disulfidptosis_10", "WP_FERROPTOSIS", "GOBP_FERROPTOSIS", "REACTOME_PYROPTOSIS", "GOBP_NECROPTOTIC_PROCESS", "REACTOME_REGULATED_NECROSIS", "GOBP_ANOIKIS")
EX <- fread(file.path(OUT, "robust_v2/external.tsv"))[cohort == "GSE14520", .(gene_set = set, genes_in_pool = k)]
S2d <- EX[match(gs, gene_set)][, source := c("Tsvetkov et al. 2022 (19 genes)", "Liu et al. 2023 (10 genes)", rep("MSigDB v2024.1", 6))]
## S3
S3 <- fread(file.path(FD, "TableS3_numbers.tsv"))
wb <- createWorkbook()
add <- function(nm, d, note) { addWorksheet(wb, nm); writeData(wb, nm, note, startRow = 1); writeData(wb, nm, d, startRow = 3, headerStyle = createStyle(textDecoration = "bold")); freezePane(wb, nm, firstActiveRow = 4) }
add("S1a_signatures", S1a, "Supplementary Table S1a. Published HCC cell-death signatures analysed (n = 149).")
add("S1b_genes", S1b, "Supplementary Table S1b. Per-gene coefficients and the source sentence from which each coefficient was extracted.")
saveWorkbook(wb, file.path(SUP, "Supplementary_Table_S1.xlsx"), overwrite = TRUE)
wb <- createWorkbook()
add("S2a_search", SL, "Supplementary Table S2a. PubMed queries (searched 15 October 2026), hits and records retrieved.")
add("S2b_corpus_keywords", MF, "Supplementary Table S2b. Keyword patterns used to record downstream analyses in the main text of the 149 source articles (case-insensitive; references excluded).")
add("S2c_modules", S2c, "Supplementary Table S2c. Template downstream modules and positivity criteria.")
add("S2d_gene_sets", S2d, "Supplementary Table S2d. Cell-death gene sets (genes in the sampling pool).")
CK <- fread(file.path(OUT, "hiend_m13_m15/checkpoint_genes.tsv"))[, .(gene)]   # the 47 checkpoint genes actually used by M13
stopifnot(nrow(CK) == 47, !anyDuplicated(CK$gene))
add("S2e_checkpoint_genes", CK, "Supplementary Table S2e. Immune checkpoint genes tested in the immune-checkpoint module (M13; n = 47).")
saveWorkbook(wb, file.path(SUP, "Supplementary_Table_S2.xlsx"), overwrite = TRUE)
wb <- createWorkbook(); add("S3_numbers", S3, "Supplementary Table S3. Every quantity reported in the text and figures, with source file and extraction rule.")
saveWorkbook(wb, file.path(SUP, "Supplementary_Table_S3.xlsx"), overwrite = TRUE)
for (f in c("S1", "S2", "S3")) stopifnot(file.exists(file.path(SUP, paste0("Supplementary_Table_", f, ".xlsx"))))
cat("SUPP TABLES OK", nrow(S1a), nrow(S1b), nrow(S2c), nrow(S3), "\n")
