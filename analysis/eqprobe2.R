if (nzchar(Sys.getenv("HCC_AUDIT_RLIB"))) .libPaths(Sys.getenv("HCC_AUDIT_RLIB")); suppressPackageStartupMessages({library(Matrix); library(data.table)})
O <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N11_out/hiend_m9_hvg")
WTs <- readRDS(file.path(O, "WT_network.rds"))
DRs <- as.data.table(readRDS(file.path(O, "ko_speedtest/ALAD.rds"))); DRo <- as.data.table(readRDS(file.path(O, "official_ALAD_diffRegulation.rds")))
TIN <- as.data.table(readRDS(file.path(O, "ko_speedtest/TINAGL1.rds")))
top <- function(d, n = 20) d[order(p.value)][gene != d$gene[1]][1:n]$gene
out <- data.table(
  check = c("sig genes shared (excl. target)", "sig genes official (excl. target)", "top20 overlap shared vs official (ALAD)",
            "Spearman FC shared vs official (ALAD)", "top20 overlap ALAD vs TINAGL1 (shared)"),
  value = c(paste(setdiff(DRs[p.adj < 0.05]$gene, "ALAD"), collapse = ","), paste(setdiff(DRo[p.adj < 0.05]$gene, "ALAD"), collapse = ","),
            length(intersect(top(DRs), top(DRo))), round(cor(DRs[order(gene)]$FC, DRo[order(gene)]$FC, method = "spearman"), 4),
            length(intersect(top(DRs), top(TIN)))))
print(out); fwrite(out, file.path(O, "equivalence_aligned_summary_ALAD.tsv"), sep = "\t")
