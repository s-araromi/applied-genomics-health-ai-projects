#!/usr/bin/env Rscript

# Automated validation for Project 06. The compact fixtures are exact excerpts
# from the authoritative files used by the full analysis.

arguments_all <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments_all, value = TRUE)
test_path <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_dir <- dirname(dirname(test_path))
fixture_dir <- file.path(dirname(test_path), "fixtures")
source(file.path(project_dir, "R", "enrichment_functions.R"))

tests_run <- 0L

check <- function(name, expression) {
  tests_run <<- tests_run + 1L
  tryCatch(
    {
      force(expression)
      message(sprintf("PASS: %s", name))
    },
    error = function(error) {
      stop(sprintf("FAIL: %s\n%s", name, conditionMessage(error)), call. = FALSE)
    }
  )
}

check("pinned source releases and checksums are unchanged", {
  stopifnot(SOURCE_FILES$go_annotations$bytes == 15041996)
  stopifnot(SOURCE_FILES$go_ontology$sha256 == "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc")
  stopifnot(SOURCE_FILES$reactome_mapping$md5 == "8ec3d8afe228764cf3111dc2eac6e838")
})

check("Project 05-derived inputs retain their verified dimensions", {
  inputs <- read_project05_inputs(file.path(project_dir, "data"))
  stopifnot(nrow(inputs$universe) == 20848L)
  stopifnot(nrow(inputs$genes) == 1328L)
  stopifnot(sum(inputs$genes$direction == "Higher in tumor") == 430L)
  stopifnot(sum(inputs$genes$direction == "Lower in tumor") == 898L)
})

terms <- read_go_terms(file.path(fixture_dir, "verified_go_terms.obo"))

check("GO OBO parsing preserves names and hierarchy relations", {
  mitosis <- terms[terms$term_id == "GO:1903047", ]
  stopifnot(mitosis$term_name == "mitotic cell cycle process")
  stopifnot(grepl("GO:0022402", mitosis$parent_ids, fixed = TRUE))
  stopifnot(grepl("GO:0000278", mitosis$parent_ids, fixed = TRUE))
})

check("GOA parsing keeps human biological-process annotations and removes NOT", {
  temporary_gaf <- tempfile(fileext = ".gaf.gz")
  input <- file(file.path(fixture_dir, "verified_go_annotations.gaf"), "rb")
  output <- gzfile(temporary_gaf, "wb")
  writeBin(readBin(input, "raw", n = 100000L), output)
  close(input)
  close(output)
  annotations <- read_go_annotations(temporary_gaf, terms)
  unlink(temporary_gaf)
  stopifnot(nrow(annotations) == 3L)
  stopifnot(!"GO:0008283" %in% annotations$term_id)
  stopifnot(all(annotations$gene_symbol %in% c("BUB1", "MKI67")))
})

check("GO annotations propagate through is_a and part_of ancestors", {
  direct <- data.frame(
    gene_symbol = "GENE",
    term_id = "GO:0140014",
    evidence_code = "IDA",
    stringsAsFactors = FALSE
  )
  propagated <- propagate_go_annotations(direct, terms)
  expected <- c("GO:0140014", "GO:0000280", "GO:1903047", "GO:0022402", "GO:0000278", "GO:0007049")
  stopifnot(all(expected %in% propagated$term_id))
})

check("Reactome parser preserves official human mappings", {
  mappings <- read_reactome_annotations(file.path(fixture_dir, "verified_reactome_mappings.tsv"))
  stopifnot(nrow(mappings) == 7L)
  stopifnot(any(mappings$gene_id == "1300" & mappings$term_id == "R-HSA-1442490"))
  stopifnot(any(mappings$gene_id == "699" & mappings$term_name == "Cell Cycle"))
})

check("hypergeometric over-representation uses the annotated universe", {
  mappings <- data.frame(
    gene_symbol = c("A", "B", "C", "D", "E", "F", "G", "H"),
    term_id = c(rep("GO:1", 4), rep("GO:2", 4)),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(term_id = c("GO:1", "GO:2"), term_name = c("one", "two"), namespace = "biological_process")
  result <- run_over_representation(
    input_genes = c("A", "B", "C"),
    universe_genes = LETTERS[1:10],
    mappings = mappings,
    term_metadata = metadata,
    gene_column = "gene_symbol",
    database = "GO Biological Process",
    direction = "Higher in tumor",
    label_lookup = setNames(LETTERS[1:10], LETTERS[1:10]),
    minimum_term_size = 1L,
    maximum_term_size = 10L,
    minimum_overlap = 1L
  )
  first <- result[result$term_id == "GO:1", ]
  stopifnot(first$background_size == 8L)
  stopifnot(first$input_size == 3L)
  stopifnot(first$overlap_count == 3L)
  stopifnot(all(result$fdr_bh >= 0 & result$fdr_bh <= 1))
})

check("GO experimental-evidence sensitivity is attached by term", {
  primary <- data.frame(term_id = c("GO:1", "GO:2"), fdr_bh = c(0.01, 0.02), enrichment_ratio = c(3, 2))
  experimental <- data.frame(
    term_id = "GO:1", input_size = 5, term_size = 10, overlap_count = 3,
    enrichment_ratio = 2, p_value = 0.001, fdr_bh = 0.01
  )
  joined <- attach_go_sensitivity(primary, experimental)
  stopifnot(joined$experimental_supported[joined$term_id == "GO:1"])
  stopifnot(!joined$experimental_supported[joined$term_id == "GO:2"])
})

check("Jaccard similarity is calculated from overlapping genes", {
  stopifnot(all.equal(jaccard_similarity(c("A", "B"), c("B", "C")), 1 / 3))
  stopifnot(jaccard_similarity(character(), character()) == 0)
})

check("representative selection removes highly redundant terms", {
  candidates <- data.frame(
    term_id = c("T1", "T2", "T3"),
    fdr_bh = c(0.001, 0.002, 0.003),
    overlapping_gene_symbols = c("A;B;C;D", "A;B;C;E", "X;Y;Z"),
    stringsAsFactors = FALSE
  )
  selected <- select_representative_terms(candidates, maximum_terms = 3L, threshold = 0.5)
  stopifnot(identical(selected$term_id, c("T1", "T3")))
})

check("verified output dimensions and summary counts remain reproducible", {
  results <- read.csv(file.path(project_dir, "outputs", "functional_enrichment_results.csv"), stringsAsFactors = FALSE)
  summary <- read.csv(file.path(project_dir, "outputs", "enrichment_summary.csv"), stringsAsFactors = FALSE)
  representatives <- read.csv(file.path(project_dir, "outputs", "representative_enriched_terms.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(results) == 4399L)
  stopifnot(nrow(representatives) == 48L)
  stopifnot(identical(summary$significant_terms, c(255L, 99L, 756L, 52L)))
  stopifnot(all(summary$representative_terms == 12L))
})

message(sprintf("All %s Project 06 tests passed.", tests_run))
