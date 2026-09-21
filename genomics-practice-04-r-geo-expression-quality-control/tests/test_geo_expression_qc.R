script_path <- "analyze_geo_expression_qc.R"
if (!file.exists(script_path)) script_path <- file.path("..", "analyze_geo_expression_qc.R")
source(script_path)

passed <- 0L
failed <- 0L

check <- function(name, expression) {
  outcome <- tryCatch(isTRUE(force(expression)), error = function(error) FALSE)
  if (outcome) {
    cat("PASS:", name, "\n")
    passed <<- passed + 1L
  } else {
    cat("FAIL:", name, "\n")
    failed <<- failed + 1L
  }
}

expect_error <- function(expression) {
  inherits(tryCatch({ force(expression); NULL }, error = identity), "error")
}

metadata <- data.frame(
  sample_accession = c("GSM1", "GSM2", "GSM3", "GSM4"),
  sample_title = c("Lung Cancer 1T", "Lung Cancer 2T",
                   "Lung Normal 1N", "Lung Normal 2N"),
  patient_id = c("1", "2", "1", "2"),
  tissue = c("Tumor", "Tumor", "Normal", "Normal"),
  stringsAsFactors = FALSE
)
expression <- matrix(c(4, 5, 6, 7, 8,
                       4.1, 5.2, 6.1, 7.2, 8.1,
                       4.2, 5.1, 6.2, 7.1, 8.2,
                       4.3, 5.3, 6.3, 7.3, 8.3), nrow = 5)
rownames(expression) <- paste0("probe", 1:5)
colnames(expression) <- metadata$sample_accession

check("Argument defaults identify the real GEO matrix",
      grepl("GSE19804_series_matrix", parse_args(character())$input))
check("Pinned source size is 21,530,998 bytes", EXPECTED_BYTES == 21530998)
check("Pinned source MD5 has 32 hexadecimal characters",
      grepl("^[0-9a-f]{32}$", EXPECTED_MD5))
check("Valid synthetic expression data pass structural validation",
      identical(validate_expression(expression, metadata), invisible(TRUE)))

bad_columns <- expression
colnames(bad_columns)[1] <- "WRONG"
check("Misaligned GEO accessions are rejected",
      expect_error(validate_expression(bad_columns, metadata)))

bad_probe <- expression
rownames(bad_probe)[2] <- rownames(bad_probe)[1]
check("Duplicate probe identifiers are rejected",
      expect_error(validate_expression(bad_probe, metadata)))

bad_value <- expression
bad_value[1, 1] <- NA_real_
check("Missing expression values are rejected",
      expect_error(validate_expression(bad_value, metadata)))

check("Review threshold uses Q1 minus two IQRs",
      isTRUE(all.equal(review_threshold(c(1, 2, 3, 4, 5)), -2)))

sample_metrics <- compute_sample_metrics(expression, metadata)
check("Sample metrics retain all four samples", nrow(sample_metrics) == 4L)
check("Sample medians are calculated correctly",
      isTRUE(all.equal(unname(sample_metrics$sample_median),
                       unname(apply(expression, 2, median)))))

pair_metrics <- compute_pair_metrics(expression, metadata)
check("Two complete synthetic pairs are recovered", nrow(pair_metrics) == 2L)
check("Pair correlations are finite", all(is.finite(pair_metrics$pair_correlation)))

cat(sprintf("\nResults: %s passed, %s failed\n", passed, failed))
if (failed > 0L) quit(status = 1L)
