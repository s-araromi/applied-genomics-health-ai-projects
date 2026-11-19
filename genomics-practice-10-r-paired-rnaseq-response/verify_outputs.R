#!/usr/bin/env Rscript
# Compare independently generated CSV tables against the delivered reference.
# Usage: Rscript --vanilla verify_outputs.R outputs rerun_outputs [--strict]
args <- commandArgs(trailingOnly = TRUE)
reference <- if (length(args)>=1) args[1] else "outputs"
candidate <- if (length(args)>=2) args[2] else "rerun_outputs"
strict <- "--strict" %in% args
files <- sort(list.files(reference, pattern = "[.]csv$"))
if (!length(files) || !dir.exists(candidate)) stop("Both output directories are required")
other <- sort(list.files(candidate, pattern = "[.]csv$"))
if (!identical(files, other)) stop("CSV file coverage differs")
for (f in files) {
  a <- read.csv(file.path(reference,f),check.names=FALSE,stringsAsFactors=FALSE)
  b <- read.csv(file.path(candidate,f),check.names=FALSE,stringsAsFactors=FALSE)
  equal <- isTRUE(all.equal(a,b,tolerance=1e-7,check.attributes=TRUE))
  if (!equal) stop("Numerical or structural mismatch: ",f,
                   "; inspect package versions and scientific changes before replacing outputs")
  byte_equal <- identical(readBin(file.path(reference,f),"raw",n=file.info(file.path(reference,f))$size),
                          readBin(file.path(candidate,f),"raw",n=file.info(file.path(candidate,f))$size))
  if (strict && !byte_equal) stop("Byte mismatch: ",f)
  cat("OK",f,if (byte_equal) "(identical bytes)" else "(numeric tolerance 1e-7)","\n")
}
cat(length(files),"CSV tables verified. Plots and session metadata require separate review.\n")
