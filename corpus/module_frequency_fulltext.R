# 语料下游模块使用频率（2026-10-18 固化版，替代 10-17 无脚本的旧表）
# 旧表 module_frequency_fulltext.tsv（分母 151）未存生成脚本，且按现有全文无法复现（单模块差值最大 37 篇）→ 弃用，改名 *_UNREPRODUCIBLE_20261017.tsv
# 本版口径（写死）：
#   语料 = signatures_index.tsv 中主分析签名（group_primary == "是" & 是否有完整系数 == "是"）对应的 PMID，每篇一次
#   文本 = raw/txt/<PMID>.txt 正文，截到第一行 "References" 之前（不计参考文献标题里的词）
#   判据 = 旧表 regex 列原样、大小写不敏感
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
