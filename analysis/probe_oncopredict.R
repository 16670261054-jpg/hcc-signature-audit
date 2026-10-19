if (nzchar(Sys.getenv("HCC_AUDIT_RLIB"))) .libPaths(Sys.getenv("HCC_AUDIT_RLIB"))
x <- deparse(oncoPredict::calcPhenotype)
cat(paste(seq_along(x), x)[c(36:70, 300:377)], sep = "\n")
