# Use of downstream modules in the corpus (fixed version; replaces an earlier table without a generating script)
# The earlier table module_frequency_fulltext.tsv (denominator 151) had no generating script and could not be reproduced from the full texts -> discarded and renamed *_UNREPRODUCIBLE_20261017.tsv
# Definitions (fixed):
#   Corpus = PMIDs of the primary signatures in signatures_index.tsv (group_primary == yes & complete coefficients == yes), each article once
#   Text = main text in raw/txt/<PMID>.txt, truncated before the first line "References" (words in reference titles not counted)
#   Criterion = regex column of the earlier table unchanged, case-insensitive
suppressPackageStartupMessages(library(data.table))
B <- file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N12_signatures")
old_f <- file.path(B, "module_frequency_fulltext.tsv")
arch <- file.path(B, "module_frequency_fulltext_UNREPRODUCIBLE_20261017.tsv")
if (!file.exists(arch)) file.copy(old_f, arch)
rx <- fread(arch, sep = "\t", quote = "")[, .(module, regex)]
IX <- fread(file.path(B, "signatures_index.tsv"), encoding = "UTF-8")
pm <- unique(IX[group_primary == "是" & 是否有完整系数 == "是"]$PMID)
body <- function(p) { f <- file.path(B, "raw/txt", paste0(p, ".txt")); if (!file.exists(f)) return(NA_character_)
  L <- readLines(f, warn = FALSE, encoding = "UTF-8"); cut <- grep("^\\S*\\s*references\\s*$", L, ignore.case = TRUE)[1]
  if (!is.na(cut)) L <- L[seq_len(cut - 1)]; tolower(paste(L, collapse = " ")) }
tx <- vapply(pm, body, ""); miss <- sum(is.na(tx)); tx <- tx[!is.na(tx)]
out <- rx[, .(n_papers = vapply(regex, function(r) sum(grepl(r, tx, perl = TRUE)), 1L)), by = .(module, regex)]
out[, `:=`(n_fulltext = length(tx), frac = n_papers / length(tx))]
setcolorder(out, c("module", "n_papers", "frac", "n_fulltext", "regex")); setorder(out, -n_papers)
fwrite(out, old_f, sep = "\t")
fwrite(data.table(n_pmid = length(pm), n_txt_found = length(tx), n_missing = miss,
                  rule = "primary signatures; body text before 'References'"), file.path(B, "module_frequency_fulltext_meta.tsv"), sep = "\t")
print(out[, .(module, n_papers, n_fulltext)])
