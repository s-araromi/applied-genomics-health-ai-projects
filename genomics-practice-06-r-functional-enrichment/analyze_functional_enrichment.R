#!/usr/bin/env Rscript

# Project 06 entry point: GO and Reactome over-representation analysis.

arguments_all <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments_all, value = TRUE)
script_path <- if (length(file_argument) == 1L) {
  normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
} else {
  normalizePath("analyze_functional_enrichment.R", mustWork = TRUE)
}
project_dir <- dirname(script_path)
source(file.path(project_dir, "R", "enrichment_functions.R"))

arguments <- commandArgs(trailingOnly = TRUE)
allowed <- "--refresh-source"
unknown <- setdiff(arguments, allowed)
if (length(unknown) > 0L) {
  stop(paste("Unknown argument(s):", paste(unknown, collapse = ", ")), call. = FALSE)
}

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Validating Project 05-derived inputs and authoritative annotation sources...")
input_data <- read_project05_inputs(data_dir)
source_paths <- ensure_sources(data_dir, refresh = "--refresh-source" %in% arguments)

message("Parsing Gene Ontology terms and human GOA annotations...")
go_terms <- read_go_terms(source_paths$go_ontology)
go_annotations <- read_go_annotations(source_paths$go_annotations, go_terms)
message("Propagating GO annotations through is_a and part_of ancestors...")
go_annotations <- propagate_go_annotations(go_annotations, go_terms)
go_term_metadata <- go_terms[go_terms$namespace == "biological_process", c("term_id", "term_name", "namespace")]
go_experimental <- go_annotations[
  go_annotations$evidence_code %in% EXPERIMENTAL_GO_EVIDENCE,
  ,
  drop = FALSE
]

message("Parsing human Reactome gene-to-pathway mappings...")
reactome <- read_reactome_annotations(source_paths$reactome_mapping)
reactome_terms <- unique(reactome[c("term_id", "term_name")])
reactome_terms$namespace <- "reactome_pathway"

symbol_lookup <- setNames(input_data$universe$gene_symbol, input_data$universe$gene_symbol)
id_to_symbol <- setNames(input_data$universe$gene_symbol, input_data$universe$gene_id)

result_sets <- list()
representative_sets <- list()
for (direction in c("Higher in tumor", "Lower in tumor")) {
  direction_genes <- input_data$genes[input_data$genes$direction == direction, ]

  go_primary <- run_over_representation(
    input_genes = direction_genes$gene_symbol,
    universe_genes = input_data$universe$gene_symbol,
    mappings = go_annotations,
    term_metadata = go_term_metadata,
    gene_column = "gene_symbol",
    database = "GO Biological Process",
    direction = direction,
    label_lookup = symbol_lookup
  )
  go_sensitivity <- run_over_representation(
    input_genes = direction_genes$gene_symbol,
    universe_genes = input_data$universe$gene_symbol,
    mappings = go_experimental,
    term_metadata = go_term_metadata,
    gene_column = "gene_symbol",
    database = "GO Biological Process",
    direction = direction,
    label_lookup = symbol_lookup
  )
  go_results <- attach_go_sensitivity(go_primary, go_sensitivity)

  reactome_results <- run_over_representation(
    input_genes = direction_genes$gene_id,
    universe_genes = input_data$universe$gene_id,
    mappings = reactome,
    term_metadata = reactome_terms,
    gene_column = "gene_id",
    database = "Reactome",
    direction = direction,
    label_lookup = id_to_symbol
  )
  reactome_results$experimental_supported <- NA

  result_sets[[paste("GO", direction)]] <- go_results
  result_sets[[paste("Reactome", direction)]] <- reactome_results
  representative_sets[[paste("GO", direction)]] <- select_representative_terms(go_results)
  representative_sets[[paste("Reactome", direction)]] <- select_representative_terms(reactome_results)
}

all_results <- bind_rows_by_name(result_sets)
representatives <- bind_rows_by_name(representative_sets)
rownames(all_results) <- NULL
rownames(representatives) <- NULL

write.csv(all_results, file.path(output_dir, "functional_enrichment_results.csv"), row.names = FALSE)
write.csv(representatives, file.path(output_dir, "representative_enriched_terms.csv"), row.names = FALSE)
summary <- summarize_results(
  all_results,
  representatives,
  input_data,
  file.path(output_dir, "enrichment_summary.csv")
)
write_source_manifest(
  source_paths,
  input_data$paths,
  file.path(data_dir, "source_manifest.csv")
)
write_report(summary, representatives, file.path(output_dir, "functional_enrichment_report.md"))
plot_enrichment(
  representatives,
  file.path(output_dir, "functional_enrichment_overview.png"),
  "png"
)
plot_enrichment(
  representatives,
  file.path(output_dir, "functional_enrichment_overview.svg"),
  "svg"
)

message(sprintf("Analyzed %s robust Project 05 genes against %s measured genes.", nrow(input_data$genes), nrow(input_data$universe)))
message(sprintf("Wrote %s full enrichment rows and %s representative terms.", nrow(all_results), nrow(representatives)))
message(paste("Outputs written to", output_dir))
