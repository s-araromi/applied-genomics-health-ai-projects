#!/usr/bin/env Rscript
# Install a Bioconductor release compatible with the running R version.
# Reproduction was verified on R 4.3.3 / DESeq2 1.42.0; newer versions can differ.
if (getRversion() < "4.3.0") stop("R 4.3 or newer is required")
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "SummarizedExperiment"), ask = FALSE, update = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
cat("\nInstalled versions:\n")
for (p in c("DESeq2", "SummarizedExperiment", "digest")) cat(p, as.character(packageVersion(p)), "\n")
