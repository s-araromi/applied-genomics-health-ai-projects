# Export every gene and all eight dexamethasone/untreated libraries from the
# versioned airway object. No FPKM values or example BAM subsets are used.
export_airway <- function(rdata_path, destination) {
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  x <- env$airway
  counts <- SummarizedExperiment::assay(x, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(x))
  genes <- as.data.frame(SummarizedExperiment::rowData(x))
  expected <- c("SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
                "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521")
  stopifnot(identical(dim(counts), c(63677L, 8L)),
            identical(colnames(counts), expected), all(meta$albut == "untrt"),
            identical(rownames(counts), genes$gene_id),
            all(counts >= 0), all(counts == floor(counts)))
  samples <- data.frame(run = meta$Run, geo_sample = meta$SampleName,
                        donor = meta$cell, dex = meta$dex,
                        albuterol = meta$albut, biosample = meta$BioSample,
                        experiment = meta$Experiment)
  # Gene symbols are historical labels; stable gene IDs remain the primary key.
  annotation <- genes[, c("gene_id", "gene_name", "gene_biotype", "seq_name",
                           "gene_seq_start", "gene_seq_end", "seq_strand")]
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            file.path(destination, "counts.csv"), row.names = FALSE, na = "")
  write.csv(samples, file.path(destination, "samples.csv"), row.names = FALSE, na = "")
  write.csv(annotation, file.path(destination, "gene_annotation.csv"),
            row.names = FALSE, na = "")
  invisible(samples)
}
