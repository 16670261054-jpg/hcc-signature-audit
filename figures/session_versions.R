# Writes software versions (the Methods "Software" section takes version numbers from here)
pk <- c("survival", "glmnet", "survcomp", "randomForestSRC", "gbm", "superpc", "MASS", "rms", "dcurves", "GSVA", "estimate",
        "fgsea", "limma", "clusterProfiler", "ConsensusClusterPlus", "oncoPredict", "Seurat", "SeuratObject", "scTenifoldKnk",
        "scTenifoldNet", "harmony", "data.table", "ggplot2", "patchwork", "ggtree", "ggtreeExtra", "circlize", "ggdist", "timeROC")
v <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
out <- data.frame(software = c("R", pk), version = c(paste(R.version$major, R.version$minor, sep = "."), v))
write.table(out, file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N28_fig/software_versions_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("HCC_AUDIT_ROOT", "."), "N28_fig/sessionInfo_R.txt"))
print(out)
