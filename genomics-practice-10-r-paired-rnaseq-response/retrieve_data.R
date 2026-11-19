#!/usr/bin/env Rscript
# Rebuild the three analytical input tables from the pinned, genuine archive.
# The temporary package is not installed; its data object is read directly.
source("R/export_airway.R")
if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest first")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install SummarizedExperiment first")
args <- commandArgs(trailingOnly = TRUE)
dest <- if (length(args)) args[1] else "retrieved_data"
if (dir.exists(dest) && length(list.files(dest, all.files = TRUE, no.. = TRUE)))
  stop("Choose an empty output directory; existing inputs are never overwritten")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
expected_hash <- "829f654fc8075ac9d810574e8a8d4672e69e847c45fa2ce598d1ffbc30c7fe19"
urls <- c("https://bioconductor.org/packages/3.23/data/experiment/src/contrib/airway_1.32.0.tar.gz",
          "https://bioconductor.posit.co/packages/3.23/data/experiment/src/contrib/airway_1.32.0.tar.gz")
work <- tempfile("airway-source-"); dir.create(work)
archive <- file.path(work, "airway_1.32.0.tar.gz")
options(timeout = 300)
ok <- FALSE
for (u in urls) {
  status <- tryCatch({download.file(u, archive, mode = "wb", method = "libcurl", quiet = TRUE); TRUE},
                     error = function(e) FALSE)
  if (status && identical(digest::digest(file = archive, algo = "sha256"), expected_hash)) {
    ok <- TRUE; message("Verified pinned archive from ", u); break
  }
}
if (!ok) stop("Pinned archive unavailable or identity changed; no substitute is permitted")
utils::untar(archive, files = c("airway/data/airway.RData", "airway/DESCRIPTION"), exdir = work)
desc <- read.dcf(file.path(work,"airway/DESCRIPTION"))
stopifnot(desc[1,"Version"] == "1.32.0", desc[1,"License"] == "LGPL")
export_airway(file.path(work, "airway/data/airway.RData"), dest)
m <- read.csv("data/source_manifest.csv")
for (f in c("counts.csv", "samples.csv", "gene_annotation.csv")) {
  expected <- m$sha256[match(paste0("data/", f), m$path)]
  actual <- digest::digest(file = file.path(dest, f), algo = "sha256")
  if (!identical(actual, expected)) stop("Export identity mismatch: ", f)
}
unlink(work, recursive = TRUE)
cat("All three exported input tables reproduce the committed SHA-256 values.\n")
