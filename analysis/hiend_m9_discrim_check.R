# M9 判别性检查（预注册 N20 变更记录 2026-10-17；用户要求"结果都一样或毫不相关就停"）—— 运行停止规则，不是论文判据
# 输入：KNK_OUT 目录下 WT_network.rds + KNK_KODIR 里的逐基因 diffRegulation（.rds 全精度）
# D1 "都是同一个结果"：在有 ≥1 个显著基因的 KO 中，任一非靶基因出现在 ≥50% 的 KO 里 → 触发（无判别性）；另报两两 Jaccard 中位数
# D2 "毫不相关"：每个有显著基因的 KO，显著基因与靶基因的网络边权 |WT[g,s]|+|WT[s,g]| vs 其余基因，Wilcoxon 单侧(greater)；
#     P<0.05 的 KO 占比 ≥ 50% → 相关性成立；有显著基因的 KO 少于 3 个时 D2 不可评估（照实标注）
# PASS = D1 未触发 且（D2 成立 或 D2 不可评估）
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("HCC_AUDIT_ROOT", ".")
O <- file.path(B, "N11_out", Sys.getenv("KNK_OUT", "hiend_m9_hvg")); KD <- file.path(O, Sys.getenv("KNK_KODIR", "ko_speedtest"))
WT <- readRDS(file.path(O, "WT_network.rds"))
fs <- list.files(KD, "\\.rds$", full.names = TRUE)
R <- rbindlist(lapply(fs, function(f) { g <- sub("\\.rds$", "", basename(f)); d <- as.data.table(readRDS(f))
  sig <- d[gene != g & p.adj < 0.05]$gene
  w <- abs(WT[g, ]) + abs(WT[, g]); w <- w[names(w) != g]
  wp <- if (length(sig) >= 1 && length(sig) < length(w)) wilcox.test(w[sig], w[setdiff(names(w), sig)], alternative = "greater", exact = FALSE)$p.value else NA_real_
  data.table(ko = g, n_sig = length(sig), sig_genes = paste(sig, collapse = ","), frac_sig_with_edge = if (length(sig)) mean(w[sig] > 0) else NA_real_,
             median_w_sig = if (length(sig)) median(w[sig]) else NA_real_, median_w_other = median(w), wilcox_P = wp) }))
fwrite(R, file.path(O, sprintf("discrim_per_ko_%s.tsv", basename(KD))), sep = "\t")
H <- R[n_sig >= 1]
sets <- lapply(strsplit(H$sig_genes, ","), unique)
freq <- if (nrow(H)) sort(table(unlist(sets)), decreasing = TRUE) / nrow(H) else numeric(0)
jac <- if (nrow(H) >= 2) { cb <- combn(nrow(H), 2); apply(cb, 2, function(ij) { a <- sets[[ij[1]]]; b <- sets[[ij[2]]]; length(intersect(a, b)) / length(union(a, b)) }) } else NA_real_
d1 <- length(freq) > 0 && max(freq) >= 0.5 && nrow(H) >= 2
d2_eval <- nrow(H) >= 3
d2 <- if (d2_eval) mean(H$wilcox_P < 0.05, na.rm = TRUE) >= 0.5 else NA
out <- data.table(kodir = basename(KD), n_ko = nrow(R), n_ko_with_sig = nrow(H),
  top_gene = if (length(freq)) names(freq)[1] else NA_character_, top_gene_freq = if (length(freq)) as.numeric(freq[1]) else NA_real_,
  top5 = if (length(freq)) paste(sprintf("%s:%.2f", names(freq)[1:min(5, length(freq))], freq[1:min(5, length(freq))]), collapse = ";") else NA_character_,
  jaccard_median = median(jac), D1_same_result_triggered = d1,
  D2_evaluable = d2_eval, frac_ko_wilcoxP_lt0.05 = if (d2_eval) mean(H$wilcox_P < 0.05, na.rm = TRUE) else NA_real_, D2_related = d2,
  PASS = !d1 && (isTRUE(d2) || !d2_eval))
fwrite(out, file.path(O, sprintf("discrim_check_%s.tsv", basename(KD))), sep = "\t"); print(out); print(R)
cat("DONE\n")
