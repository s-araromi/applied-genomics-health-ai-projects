#!/usr/bin/env Rscript

# Validate coordinate conversion, assembly selection, source-backed fixtures,
# and integration with the genuine Project 01 RefSeq results. No simulated
# nucleotide sequences or fabricated genomic records are used.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_argument) != 1L) {
  stop("Run this test file with: Rscript tests/test_gene_annotations.R")
}

test_directory <- dirname(
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
)
project_directory <- dirname(test_directory)

source(file.path(project_directory, "analyze_gene_annotations.R"))

passed <- 0L
failed <- 0L

run_test <- function(label, expression) {
  outcome <- tryCatch(
    {
      force(expression)
      TRUE
    },
    error = function(error) {
      cat(sprintf("FAIL: %s\n      %s\n", label, conditionMessage(error)))
      FALSE
    }
  )

  if (outcome) {
    passed <<- passed + 1L
    cat(sprintf("PASS: %s\n", label))
  } else {
    failed <<- failed + 1L
  }
}

run_test("Documented NCBI Gene records are authentic and ordered", {
  annotations <- load_verified_fixture(project_directory)
  stopifnot(identical(as.character(annotations$gene), names(GENE_IDS)))
  stopifnot(identical(as.character(annotations$gene_id), unname(GENE_IDS)))
  stopifnot(all(startsWith(annotations$source_url, "https://www.ncbi.nlm.nih.gov/gene/")))
  stopifnot(all(annotations$assembly == "GRCh38.p14"))
})

run_test("Real MTHFR minus-strand coordinates are normalized correctly", {
  coordinates <- normalize_gene_coordinates(11805963, 11785722)
  stopifnot(coordinates$start_1_based == 11785723)
  stopifnot(coordinates$end_1_based == 11805964)
  stopifnot(identical(coordinates$strand, "-"))
  stopifnot(coordinates$genomic_span_bp == 20242)
})

run_test("Real ACE plus-strand coordinates are normalized correctly", {
  coordinates <- normalize_gene_coordinates(63477060, 63498372)
  stopifnot(coordinates$start_1_based == 63477061)
  stopifnot(coordinates$end_1_based == 63498373)
  stopifnot(identical(coordinates$strand, "+"))
  stopifnot(coordinates$genomic_span_bp == 21313)
})

run_test("Real AGT and NOS3 genomic spans match their NCBI Gene pages", {
  annotations <- load_verified_fixture(project_directory)
  agt <- annotations[annotations$gene == "AGT", , drop = FALSE]
  nos3 <- annotations[annotations$gene == "NOS3", , drop = FALSE]
  stopifnot(agt$genomic_span_bp == 43061, identical(agt$strand, "-"))
  stopifnot(nos3$genomic_span_bp == 23572, identical(nos3$strand, "+"))
})

run_test("Invalid chromosome positions are rejected", {
  negative_rejected <- inherits(
    try(normalize_gene_coordinates(-1, 63498372), silent = TRUE),
    "try-error"
  )
  fractional_rejected <- inherits(
    try(normalize_gene_coordinates(63477060.5, 63498372), silent = TRUE),
    "try-error"
  )
  stopifnot(negative_rejected, fractional_rejected)
})

run_test("The NCBI request contains all four genuine Gene IDs", {
  url <- build_ncbi_url("researcher@example.org")
  stopifnot(grepl("db=gene", url, fixed = TRUE))
  stopifnot(grepl("4524,1636,183,4846", url, fixed = TRUE))
  stopifnot(grepl("researcher%40example.org", url, fixed = TRUE))
})

run_test("Project 01 and Project 02 real datasets integrate correctly", {
  annotations <- load_verified_fixture(project_directory)
  transcripts <- load_project_01_results(project_directory)
  combined <- combine_gene_and_transcript_data(annotations, transcripts)

  stopifnot(nrow(combined) == 4L)
  stopifnot(identical(as.character(combined$gene), names(GENE_IDS)))
  stopifnot(combined$refseq_accession[combined$gene == "AGT"] == "NM_001382817.3")
  stopifnot(combined$genomic_to_transcript_ratio[combined$gene == "AGT"] == 20.05)
  stopifnot(combined$transcript_fraction_of_genomic_span_percent[combined$gene == "MTHFR"] == 34.67)
  stopifnot(combined$exons_per_10kb[combined$gene == "NOS3"] == 11.88)
})

run_test("A missing real Project 01 gene is detected", {
  annotations <- load_verified_fixture(project_directory)
  transcripts <- load_project_01_results(project_directory)
  incomplete <- transcripts[transcripts$gene != "ACE", , drop = FALSE]
  rejected <- inherits(
    try(combine_gene_and_transcript_data(annotations, incomplete), silent = TRUE),
    "try-error"
  )
  stopifnot(rejected)
})

cat(sprintf("\nResults: %d passed, %d failed\n", passed, failed))

if (failed > 0L) {
  quit(status = 1L)
}
