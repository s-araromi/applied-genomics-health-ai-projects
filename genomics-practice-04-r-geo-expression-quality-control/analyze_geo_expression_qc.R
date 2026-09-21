#!/usr/bin/env Rscript

DATA_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz"
EXPECTED_BYTES <- 21530998
EXPECTED_MD5 <- "c817eb7b8b74d1cbdd4db7580fb9ea08"

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  result <- list(input = file.path("data", "GSE19804_series_matrix.txt.gz"),
                 output_dir = "outputs")
  index <- 1L
  while (index <= length(args)) {
    if (args[[index]] == "--input") {
      index <- index + 1L
      if (index > length(args)) stop("--input requires a path")
      result$input <- args[[index]]
    } else if (args[[index]] == "--output-dir") {
      index <- index + 1L
      if (index > length(args)) stop("--output-dir requires a path")
      result$output_dir <- args[[index]]
    } else {
      stop(sprintf("Unknown argument: %s", args[[index]]))
    }
    index <- index + 1L
  }
  result
}

download_matrix <- function(path, url = DATA_URL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial")
  on.exit(unlink(temporary), add = TRUE)
  message("Downloading the complete GSE19804 processed expression matrix...")
  download.file(url, temporary, mode = "wb", quiet = FALSE)
  file.rename(temporary, path)
}

validate_source_file <- function(path, expected_bytes = EXPECTED_BYTES,
                                 expected_md5 = EXPECTED_MD5) {
  if (!file.exists(path)) stop(sprintf("Input file not found: %s", path))
  observed_bytes <- file.info(path)$size
  if (!identical(as.numeric(observed_bytes), as.numeric(expected_bytes))) {
    stop(sprintf("Source byte count differs: observed %s; expected %s",
                 observed_bytes, expected_bytes))
  }
  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(tolower(observed_md5), tolower(expected_md5))) {
    stop(sprintf("Source MD5 differs: observed %s; expected %s",
                 observed_md5, expected_md5))
  }
  invisible(TRUE)
}

strip_quotes <- function(values) sub('^"(.*)"$', "\\1", values)

read_sample_metadata <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("Series matrix table marker was not found")
    if (identical(line, "!series_matrix_table_begin")) break
    header <- c(header, line)
  }
  extract <- function(prefix) {
    matches <- header[startsWith(header, paste0(prefix, "\t"))]
    if (length(matches) != 1L) stop(sprintf("Expected one %s metadata row", prefix))
    strip_quotes(strsplit(matches, "\t", fixed = TRUE)[[1L]][-1L])
  }
  titles <- extract("!Sample_title")
  accessions <- extract("!Sample_geo_accession")
  sources <- extract("!Sample_source_name_ch1")
  if (length(unique(c(length(titles), length(accessions), length(sources)))) != 1L) {
    stop("Sample metadata fields have inconsistent lengths")
  }
  tissue <- ifelse(grepl("primary tumor", sources, ignore.case = TRUE),
                   "Tumor",
                   ifelse(grepl("adjacent normal", sources, ignore.case = TRUE),
                          "Normal", NA_character_))
  patient_id <- sub(".* ([0-9]+)[TN]$", "\\1", titles)
  if (anyNA(tissue) || any(!grepl("^[0-9]+$", patient_id))) {
    stop("Could not derive tissue group or patient identifier")
  }
  data.frame(sample_accession = accessions, sample_title = titles,
             patient_id = patient_id, tissue = tissue,
             stringsAsFactors = FALSE)
}

read_expression_matrix <- function(path) {
  table <- read.delim(gzfile(path), comment.char = "!", quote = "\"",
                      check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(table) < 2L) stop("Expression table contains no sample columns")
  probe_id <- table[[1L]]
  expression <- as.matrix(table[-1L])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe_id
  expression
}

validate_expression <- function(expression, metadata) {
  if (ncol(expression) != nrow(metadata)) {
    stop("Expression columns do not match sample metadata")
  }
  if (!identical(colnames(expression), metadata$sample_accession)) {
    stop("Expression columns are not aligned with GEO accessions")
  }
  if (anyDuplicated(rownames(expression))) stop("Duplicate probe identifiers detected")
  if (anyDuplicated(colnames(expression))) stop("Duplicate sample identifiers detected")
  if (any(!is.finite(expression))) stop("Missing or non-finite expression values detected")
  invisible(TRUE)
}

review_threshold <- function(values, multiplier = 2) {
  quartiles <- quantile(values, c(0.25, 0.75), names = FALSE)
  quartiles[[1L]] - multiplier * (quartiles[[2L]] - quartiles[[1L]])
}

compute_sample_metrics <- function(expression, metadata) {
  result <- metadata
  result$sample_median <- apply(expression, 2L, median)
  result$sample_iqr <- apply(expression, 2L, IQR)
  result$sample_min <- apply(expression, 2L, min)
  result$sample_max <- apply(expression, 2L, max)
  result$mean_within_tissue_correlation <- NA_real_
  result$review_threshold <- NA_real_
  result$review_flag <- FALSE
  for (group in unique(metadata$tissue)) {
    indices <- which(metadata$tissue == group)
    correlations <- cor(expression[, indices, drop = FALSE], method = "pearson")
    diag(correlations) <- NA_real_
    group_means <- colMeans(correlations, na.rm = TRUE)
    threshold <- review_threshold(group_means)
    result$mean_within_tissue_correlation[indices] <- group_means
    result$review_threshold[indices] <- threshold
    result$review_flag[indices] <- group_means < threshold
  }
  result
}

compute_pair_metrics <- function(expression, metadata) {
  patients <- unique(metadata$patient_id)
  rows <- lapply(patients, function(patient) {
    tumor <- which(metadata$patient_id == patient & metadata$tissue == "Tumor")
    normal <- which(metadata$patient_id == patient & metadata$tissue == "Normal")
    if (length(tumor) != 1L || length(normal) != 1L) {
      stop(sprintf("Patient %s does not have exactly one tumor-normal pair", patient))
    }
    data.frame(patient_id = patient,
               tumor_accession = metadata$sample_accession[tumor],
               normal_accession = metadata$sample_accession[normal],
               pair_correlation = cor(expression[, tumor], expression[, normal]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

build_summary <- function(expression, metadata, sample_metrics, pair_metrics) {
  data.frame(
    metric = c("probe_sets", "samples", "tumor_samples", "normal_samples",
               "complete_pairs", "missing_or_nonfinite_values",
               "duplicate_probe_ids", "duplicate_sample_ids", "expression_min",
               "expression_max", "sample_median_min", "sample_median_max",
               "median_pair_correlation", "samples_marked_for_review",
               "samples_automatically_excluded"),
    value = c(nrow(expression), ncol(expression), sum(metadata$tissue == "Tumor"),
              sum(metadata$tissue == "Normal"), nrow(pair_metrics),
              sum(!is.finite(expression)), anyDuplicated(rownames(expression)),
              anyDuplicated(colnames(expression)), min(expression), max(expression),
              min(sample_metrics$sample_median), max(sample_metrics$sample_median),
              median(pair_metrics$pair_correlation), sum(sample_metrics$review_flag), 0),
    stringsAsFactors = FALSE
  )
}

write_report <- function(path, summary, sample_metrics, pair_metrics) {
  value <- function(name) summary$value[summary$metric == name]
  flagged <- sample_metrics[sample_metrics$review_flag,
                            c("sample_accession", "sample_title", "tissue",
                              "mean_within_tissue_correlation")]
  lines <- c(
    "# Paired GEO Expression-Matrix Quality Assessment",
    "",
    "## Dataset",
    "",
    "- GEO series: `GSE19804`",
    "- Platform: `GPL570`",
    "- Design: paired lung tumor and adjacent-normal tissues",
    "- Processed values: authors' quantile-normalized log2 expression matrix",
    "",
    "## Results",
    "",
    sprintf("- Probe sets: %s", format(value("probe_sets"), big.mark = ",", scientific = FALSE)),
    sprintf("- Samples: %s", value("samples")),
    sprintf("- Complete tumor-normal pairs: %s", value("complete_pairs")),
    sprintf("- Expression range: %.4f to %.4f", value("expression_min"), value("expression_max")),
    sprintf("- Sample-median range: %.5f to %.5f", value("sample_median_min"), value("sample_median_max")),
    sprintf("- Median within-pair correlation: %.4f", value("median_pair_correlation")),
    sprintf("- Samples marked for review: %s", value("samples_marked_for_review")),
    "- Samples automatically excluded: 0",
    "",
    "## Samples marked for review",
    "",
    "| GEO accession | Sample | Tissue | Mean within-tissue correlation |",
    "| --- | --- | --- | ---: |",
    apply(flagged, 1L, function(row) sprintf("| %s | %s | %s | %.4f |",
                                             row[[1L]], row[[2L]], row[[3L]],
                                             as.numeric(row[[4L]]))),
    "",
    "## Interpretation",
    "",
    "The narrow sample-median range is consistent with the normalization reported",
    "for the processed GEO matrix. Review flags identify unusually low mean",
    "correlation within a tissue group using a Q1 minus 2 x IQR rule. Flagged",
    "samples remain in the dataset because the differences may reflect technical",
    "variation, participant-specific biology, or tumor heterogeneity.",
    "",
    "Quality-control patterns do not establish differential expression, disease",
    "causation, diagnostic performance, or clinical validity."
  )
  writeLines(lines, path, useBytes = TRUE)
}

plot_qc <- function(path_png, path_svg, expression, sample_metrics, pair_metrics) {
  draw <- function() {
    old <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.7, 1), oma = c(0, 0, 2, 0))
    on.exit(par(old), add = TRUE)
    colors <- ifelse(sample_metrics$tissue == "Tumor", "#B64926", "#2A6F97")
    plot(sample_metrics$sample_median, pch = 19, col = colors,
         xlab = "Sample order", ylab = "Median log2 expression",
         main = "Sample medians")
    legend("topright", c("Tumor", "Normal"), pch = 19,
           col = c("#B64926", "#2A6F97"), bty = "n")
    boxplot(split(sample_metrics$sample_iqr, sample_metrics$tissue),
            col = c("#9EC5E0", "#E6A184"), ylab = "Expression IQR",
            main = "Within-sample spread")
    plot(sample_metrics$mean_within_tissue_correlation, pch = 19, col = colors,
         xlab = "Sample order", ylab = "Mean correlation",
         main = "Within-tissue similarity")
    flagged <- which(sample_metrics$review_flag)
    points(flagged, sample_metrics$mean_within_tissue_correlation[flagged],
           pch = 1, cex = 1.7, lwd = 2, col = "#7F1D1D")
    hist(pair_metrics$pair_correlation, breaks = 12, col = "#6A8D73",
         border = "white", xlab = "Pearson correlation",
         main = "Tumor-normal pair correlation")
    abline(v = median(pair_metrics$pair_correlation), lty = 2, lwd = 2)
    mtext("GSE19804 expression-matrix quality assessment", outer = TRUE,
          cex = 1.25, font = 2)
  }
  png(path_png, width = 1800, height = 1250, res = 180)
  draw()
  dev.off()
  svg(path_svg, width = 10, height = 7)
  draw()
  dev.off()
}

main <- function() {
  args <- parse_args()
  if (!file.exists(args$input)) download_matrix(args$input)
  validate_source_file(args$input)
  message("Reading the complete processed expression matrix...")
  metadata <- read_sample_metadata(args$input)
  expression <- read_expression_matrix(args$input)
  validate_expression(expression, metadata)
  sample_metrics <- compute_sample_metrics(expression, metadata)
  pair_metrics <- compute_pair_metrics(expression, metadata)
  summary <- build_summary(expression, metadata, sample_metrics, pair_metrics)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(sample_metrics, file.path(args$output_dir, "sample_quality_metrics.csv"),
            row.names = FALSE, quote = TRUE)
  write.csv(pair_metrics, file.path(args$output_dir, "pair_correlations.csv"),
            row.names = FALSE, quote = TRUE)
  write.csv(summary, file.path(args$output_dir, "expression_qc_summary.csv"),
            row.names = FALSE, quote = TRUE)
  write_report(file.path(args$output_dir, "expression_qc_report.md"),
               summary, sample_metrics, pair_metrics)
  plot_qc(file.path(args$output_dir, "expression_qc_overview.png"),
          file.path(args$output_dir, "expression_qc_overview.svg"),
          expression, sample_metrics, pair_metrics)
  cat(sprintf("Probe sets: %s\n", format(nrow(expression), big.mark = ",")))
  cat(sprintf("Samples: %s (%s complete pairs)\n", ncol(expression), nrow(pair_metrics)))
  cat(sprintf("Expression range: %.4f to %.4f\n", min(expression), max(expression)))
  cat(sprintf("Median within-pair correlation: %.4f\n", median(pair_metrics$pair_correlation)))
  cat(sprintf("Samples marked for review: %s; automatically excluded: 0\n",
              sum(sample_metrics$review_flag)))
  cat(sprintf("Outputs written to %s\n", normalizePath(args$output_dir)))
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(error) {
    message("Analysis stopped: ", conditionMessage(error))
    quit(status = 1L)
  })
}
