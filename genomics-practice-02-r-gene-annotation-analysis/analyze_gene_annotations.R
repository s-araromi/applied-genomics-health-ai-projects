#!/usr/bin/env Rscript

# Compare genuine NCBI Gene annotations with the RefSeq transcript measurements
# produced in Project 01. Live API retrieval requires jsonlite; the explicitly
# documented, source-verifiable offline snapshot requires only base R.

GENE_IDS <- c(MTHFR = "4524", ACE = "1636", AGT = "183", NOS3 = "4846")

# Pin the primary GRCh38 chromosome accessions so alternative assemblies or
# unlocalized scaffolds cannot silently replace the intended reference records.
GRCH38_CHROMOSOME_ACCESSIONS <- c(
  MTHFR = "NC_000001.11",
  ACE = "NC_000017.11",
  AGT = "NC_000001.11",
  NOS3 = "NC_000007.14"
)

# These real cross-references are displayed on the corresponding NCBI Gene pages.
ENSEMBL_GENE_IDS <- c(
  MTHFR = "ENSG00000177000",
  ACE = "ENSG00000159640",
  AGT = "ENSG00000135744",
  NOS3 = "ENSG00000164867"
)

NCBI_ESUMMARY_URL <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  "esummary.fcgi"
)

parse_arguments <- function(arguments) {
  options <- list(email = NULL, offline_fixture = FALSE, help = FALSE)
  position <- 1L

  while (position <= length(arguments)) {
    argument <- arguments[[position]]

    if (identical(argument, "--offline-fixture")) {
      options$offline_fixture <- TRUE
    } else if (identical(argument, "--help") || identical(argument, "-h")) {
      options$help <- TRUE
    } else if (identical(argument, "--email")) {
      if (position == length(arguments)) {
        stop("The --email option must be followed by an email address.")
      }

      position <- position + 1L
      options$email <- arguments[[position]]
    } else if (startsWith(argument, "--email=")) {
      options$email <- substring(argument, first = 9L)
    } else {
      stop(sprintf("Unrecognized option: %s", argument))
    }

    position <- position + 1L
  }

  if (!is.null(options$email) && !grepl("@", options$email, fixed = TRUE)) {
    stop("The --email value does not resemble a valid email address.")
  }

  options
}

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript analyze_gene_annotations.R --email you@example.com\n",
    "  Rscript analyze_gene_annotations.R --offline-fixture\n",
    "  Rscript tests/test_gene_annotations.R\n",
    sep = ""
  )
}

find_script_directory <- function() {
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) != 1L) {
    stop("Run the analysis with Rscript so its project directory can be found.")
  }

  script_path <- sub("^--file=", "", script_argument[[1L]])
  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}

build_ncbi_url <- function(email = NULL) {
  query <- paste0(
    "?db=gene",
    "&id=", paste(unname(GENE_IDS), collapse = ","),
    "&retmode=json",
    "&tool=applied_genomics_health_ai_projects"
  )

  if (!is.null(email)) {
    # Reserved characters are encoded so the address remains a single URL value.
    query <- paste0(
      query,
      "&email=",
      utils::URLencode(email, reserved = TRUE)
    )
  }

  paste0(NCBI_ESUMMARY_URL, query)
}

download_gene_summaries <- function(email = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Live NCBI retrieval requires jsonlite. Install it with:\n",
      "Rscript -e 'install.packages(\"jsonlite\", ",
      "repos=\"https://cloud.r-project.org\")'\n",
      "Alternatively, use the documented real-data --offline-fixture option."
    )
  }

  previous_timeout <- getOption("timeout")
  options(timeout = max(60, previous_timeout))
  on.exit(options(timeout = previous_timeout), add = TRUE)

  response_text <- paste(
    readLines(build_ncbi_url(email), warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )

  response <- jsonlite::fromJSON(response_text, simplifyVector = FALSE)

  if (is.null(response$result)) {
    stop("NCBI did not return a recognizable Gene ESummary result.")
  }

  list(result = response$result, raw_json = response_text)
}

normalize_gene_coordinates <- function(chrstart, chrstop) {
  raw_start <- suppressWarnings(as.numeric(chrstart))
  raw_stop <- suppressWarnings(as.numeric(chrstop))

  if (
    length(raw_start) != 1L || length(raw_stop) != 1L ||
      !is.finite(raw_start) || !is.finite(raw_stop) ||
      raw_start < 0 || raw_stop < 0 ||
      raw_start %% 1 != 0 || raw_stop %% 1 != 0
  ) {
    stop("NCBI chromosome coordinates must be non-negative whole numbers.")
  }

  # NCBI Gene ESummary ChrStart/ChrStop values are zero-based. For a minus-
  # strand gene the biological 5-prime ChrStart is larger than ChrStop, so
  # normalize position and infer strand before converting to one-based values.
  list(
    start_1_based = min(raw_start, raw_stop) + 1,
    end_1_based = max(raw_start, raw_stop) + 1,
    strand = if (raw_start <= raw_stop) "+" else "-",
    genomic_span_bp = abs(raw_stop - raw_start) + 1
  )
}

select_grch38_location <- function(record, gene) {
  genomic_locations <- record$genomicinfo

  if (is.null(genomic_locations) || length(genomic_locations) == 0L) {
    stop(sprintf("NCBI returned no genomic coordinates for %s.", gene))
  }

  expected_accession <- unname(GRCH38_CHROMOSOME_ACCESSIONS[[gene]])

  for (location in genomic_locations) {
    if (identical(as.character(location$chraccver), expected_accession)) {
      return(location)
    }
  }

  stop(
    sprintf(
      "%s did not include the expected GRCh38 chromosome accession %s.",
      gene,
      expected_accession
    )
  )
}

extract_gene_annotation <- function(record, expected_gene, expected_gene_id) {
  if (is.null(record)) {
    stop(sprintf("NCBI did not return Gene ID %s.", expected_gene_id))
  }

  if (!identical(as.character(record$name), expected_gene)) {
    stop(sprintf("Gene ID %s did not match %s.", expected_gene_id, expected_gene))
  }

  if (!identical(as.character(record$organism$scientificname), "Homo sapiens")) {
    stop(sprintf("Gene ID %s is not a human reference record.", expected_gene_id))
  }

  location <- select_grch38_location(record, expected_gene)
  coordinates <- normalize_gene_coordinates(location$chrstart, location$chrstop)
  exon_count <- suppressWarnings(as.integer(location$exoncount))

  if (length(exon_count) != 1L || is.na(exon_count) || exon_count <= 0L) {
    stop(sprintf("NCBI returned an invalid exon count for %s.", expected_gene))
  }

  data.frame(
    gene = expected_gene,
    gene_id = as.integer(expected_gene_id),
    full_name = as.character(record$description),
    chromosome = as.character(location$chrloc),
    cytoband = as.character(record$maplocation),
    chromosome_accession = as.character(location$chraccver),
    start_1_based = coordinates$start_1_based,
    end_1_based = coordinates$end_1_based,
    strand = coordinates$strand,
    exon_count = exon_count,
    genomic_span_bp = coordinates$genomic_span_bp,
    ensembl_gene_id = unname(ENSEMBL_GENE_IDS[[expected_gene]]),
    assembly = "GRCh38.p14",
    assembly_accession = "GCF_000001405.40",
    source_url = sprintf("https://www.ncbi.nlm.nih.gov/gene/%s", expected_gene_id),
    stringsAsFactors = FALSE
  )
}

parse_live_annotations <- function(summary_result) {
  rows <- vector("list", length(GENE_IDS))

  for (position in seq_along(GENE_IDS)) {
    gene <- names(GENE_IDS)[[position]]
    gene_id <- unname(GENE_IDS[[position]])
    rows[[position]] <- extract_gene_annotation(
      summary_result[[gene_id]],
      expected_gene = gene,
      expected_gene_id = gene_id
    )
  }

  do.call(rbind, rows)
}

load_verified_fixture <- function(project_directory) {
  fixture_path <- file.path(
    project_directory,
    "tests",
    "fixtures",
    "verified_ncbi_gene_annotations.csv"
  )

  if (!file.exists(fixture_path)) {
    stop(sprintf("The documented real-data snapshot is missing: %s", fixture_path))
  }

  snapshot <- read.csv(fixture_path, stringsAsFactors = FALSE)

  if (!identical(as.character(snapshot$gene), names(GENE_IDS))) {
    stop("The verified gene snapshot does not contain the four expected records.")
  }

  rows <- vector("list", nrow(snapshot))

  for (position in seq_len(nrow(snapshot))) {
    row <- snapshot[position, , drop = FALSE]

    # Reconstruct the documented zero-based API orientation from real,
    # independently verifiable one-based coordinates on the NCBI Gene page.
    raw_start <- if (row$strand == "+") {
      row$start_1_based - 1
    } else {
      row$end_1_based - 1
    }

    raw_stop <- if (row$strand == "+") {
      row$end_1_based - 1
    } else {
      row$start_1_based - 1
    }

    coordinates <- normalize_gene_coordinates(raw_start, raw_stop)

    if (
      coordinates$start_1_based != row$start_1_based ||
        coordinates$end_1_based != row$end_1_based ||
        !identical(coordinates$strand, row$strand)
    ) {
      stop(sprintf("Coordinate validation failed for verified snapshot row %s.", row$gene))
    }

    row$genomic_span_bp <- coordinates$genomic_span_bp
    row$source_checked_date <- NULL
    rows[[position]] <- row
  }

  annotations <- do.call(rbind, rows)
  annotations[, c(
    "gene", "gene_id", "full_name", "chromosome", "cytoband",
    "chromosome_accession", "start_1_based", "end_1_based", "strand",
    "exon_count", "genomic_span_bp", "ensembl_gene_id", "assembly",
    "assembly_accession", "source_url"
  )]
}

load_project_01_results <- function(project_directory) {
  results_path <- file.path(
    dirname(project_directory),
    "genomics-practice-01-python-gc-content-refseq",
    "outputs",
    "gc_content_results.csv"
  )

  if (!file.exists(results_path)) {
    stop(
      "Project 01 results were not found. Run its Python analysis first or ",
      "restore outputs/gc_content_results.csv."
    )
  }

  results <- read.csv(results_path, stringsAsFactors = FALSE)
  required_columns <- c("gene", "accession", "length_nt", "gc_percent")

  if (!all(required_columns %in% names(results))) {
    stop("The Project 01 CSV is missing required gene, accession, length, or GC columns.")
  }

  if (anyDuplicated(results$gene)) {
    stop("The Project 01 CSV contains duplicate gene records.")
  }

  results[, required_columns]
}

combine_gene_and_transcript_data <- function(annotations, transcript_results) {
  missing_genes <- setdiff(names(GENE_IDS), transcript_results$gene)

  if (length(missing_genes) > 0L) {
    stop(
      sprintf(
        "Project 01 results do not contain: %s",
        paste(missing_genes, collapse = ", ")
      )
    )
  }

  combined <- merge(annotations, transcript_results, by = "gene", sort = FALSE)
  combined <- combined[match(names(GENE_IDS), combined$gene), , drop = FALSE]
  rownames(combined) <- NULL

  names(combined)[names(combined) == "accession"] <- "refseq_accession"
  names(combined)[names(combined) == "length_nt"] <- "transcript_length_nt"
  names(combined)[names(combined) == "gc_percent"] <- "transcript_gc_percent"

  if (any(combined$transcript_length_nt <= 0)) {
    stop("Every RefSeq transcript length must be greater than zero.")
  }

  combined$genomic_to_transcript_ratio <- round(
    combined$genomic_span_bp / combined$transcript_length_nt,
    digits = 2
  )

  combined$transcript_fraction_of_genomic_span_percent <- round(
    (combined$transcript_length_nt / combined$genomic_span_bp) * 100,
    digits = 2
  )

  combined$exons_per_10kb <- round(
    (combined$exon_count / combined$genomic_span_bp) * 10000,
    digits = 2
  )

  combined
}

write_results_csv <- function(results, output_directory) {
  output_path <- file.path(output_directory, "gene_annotation_results.csv")
  write.csv(results, output_path, row.names = FALSE, quote = TRUE)
  output_path
}

write_structure_plot <- function(results, output_directory) {
  output_path <- file.path(output_directory, "gene_vs_transcript_lengths.svg")
  grDevices::svg(filename = output_path, width = 10.5, height = 6.3)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(
    mar = c(5.5, 5.2, 4.8, 1.8),
    family = "sans",
    bg = "#F8FAFC",
    fg = "#334155",
    col.axis = "#475569",
    col.lab = "#334155"
  )

  bar_values <- rbind(
    `GRCh38 genomic span` = results$genomic_span_bp,
    `Selected RefSeq transcript` = results$transcript_length_nt
  )

  graphics::barplot(
    bar_values,
    beside = TRUE,
    names.arg = results$gene,
    col = c("#2563EB", "#0F766E"),
    border = NA,
    ylim = c(0, max(bar_values) * 1.27),
    ylab = "Length (bases)",
    main = "Genomic gene span versus selected RefSeq transcript length",
    cex.main = 1.08,
    cex.names = 1.04
  )

  graphics::legend(
    "topright",
    legend = rownames(bar_values),
    fill = c("#2563EB", "#0F766E"),
    border = NA,
    bty = "n",
    cex = 0.87
  )

  graphics::mtext(
    "Real NCBI Gene GRCh38 annotations integrated with Project 01 RefSeq results",
    side = 1,
    line = 4.0,
    col = "#64748B",
    cex = 0.78
  )

  output_path
}

write_markdown_report <- function(results, output_directory, data_source_mode) {
  output_path <- file.path(output_directory, "gene_structure_report.md")
  largest <- results[which.max(results$genomic_span_bp), , drop = FALSE]
  most_exons <- results[which.max(results$exon_count), , drop = FALSE]
  largest_ratio <- results[which.max(results$genomic_to_transcript_ratio), , drop = FALSE]

  lines <- c(
    "# Comparative Genomic Annotation and Gene-Transcript Structure",
    "",
    sprintf("Data source: %s.", data_source_mode),
    "Reference assembly: GRCh38.p14 (GCF_000001405.40).",
    "",
    "| Gene | Chromosome | GRCh38 position | Strand | Exons | Genomic span (bp) | RefSeq transcript (nt) | Span/transcript |",
    "| --- | --- | --- | --- | ---: | ---: | ---: | ---: |"
  )

  for (position in seq_len(nrow(results))) {
    row <- results[position, , drop = FALSE]
    lines <- c(
      lines,
      sprintf(
        "| %s | %s | %s-%s | %s | %d | %s | %s | %.2f |",
        row$gene,
        row$chromosome,
        format(row$start_1_based, big.mark = ",", scientific = FALSE),
        format(row$end_1_based, big.mark = ",", scientific = FALSE),
        row$strand,
        row$exon_count,
        format(row$genomic_span_bp, big.mark = ",", scientific = FALSE),
        format(row$transcript_length_nt, big.mark = ",", scientific = FALSE),
        row$genomic_to_transcript_ratio
      )
    )
  }

  lines <- c(
    lines,
    "",
    sprintf(
      "- Largest genomic span: **%s (%s bp)**.",
      largest$gene,
      format(largest$genomic_span_bp, big.mark = ",", scientific = FALSE)
    ),
    sprintf("- Highest annotated exon count: **%s (%d exons)**.", most_exons$gene, most_exons$exon_count),
    sprintf(
      "- Highest genomic-to-transcript ratio: **%s (%.2f)**.",
      largest_ratio$gene,
      largest_ratio$genomic_to_transcript_ratio
    ),
    "",
    "## Interpretation",
    "",
    "Genomic spans include the complete annotated chromosome interval, whereas",
    "the selected RefSeq transcript is one processed RNA reference. Their ratio",
    "is descriptive and must not be interpreted as expression, disease risk,",
    "intron count, total exonic coverage, or a clinical diagnostic measurement.",
    "The NCBI Gene exon count refers to the gene-level annotation and does not",
    "necessarily equal the exon count of the single selected transcript.",
    "",
    "## Data sources",
    "",
    "- NCBI Gene: https://www.ncbi.nlm.nih.gov/gene/",
    "- NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/",
    "- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25499/",
    ""
  )

  writeLines(lines, con = output_path, useBytes = TRUE)
  output_path
}

print_console_summary <- function(results) {
  cat("\nComparative genomic annotation and transcript structure\n\n")

  for (position in seq_len(nrow(results))) {
    row <- results[position, , drop = FALSE]
    cat(
      sprintf(
        "%-5s | chr%-2s | strand %s | %2d exons | %6s bp | transcript %5s nt | ratio %.2f\n",
        row$gene,
        row$chromosome,
        row$strand,
        row$exon_count,
        format(row$genomic_span_bp, big.mark = ",", scientific = FALSE),
        format(row$transcript_length_nt, big.mark = ",", scientific = FALSE),
        row$genomic_to_transcript_ratio
      )
    )
  }
}

main <- function() {
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))

  if (arguments$help) {
    print_usage()
    return(invisible(NULL))
  }

  project_directory <- find_script_directory()
  output_directory <- file.path(project_directory, "outputs")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  if (arguments$offline_fixture) {
    cat("Loading the documented, genuine NCBI Gene annotation snapshot...\n")
    annotations <- load_verified_fixture(project_directory)
    source_mode <- "source-verified NCBI Gene annotation snapshot"
  } else {
    if (is.null(arguments$email)) {
      cat("Note: NCBI recommends supplying --email your-real-email@example.com\n")
    }

    cat("Downloading four genuine human NCBI Gene summaries...\n")
    response <- download_gene_summaries(arguments$email)
    annotations <- parse_live_annotations(response$result)
    writeLines(
      response$raw_json,
      con = file.path(output_directory, "ncbi_gene_esummary.json"),
      useBytes = TRUE
    )
    source_mode <- "live NCBI Gene ESummary records"
  }

  transcript_results <- load_project_01_results(project_directory)
  results <- combine_gene_and_transcript_data(annotations, transcript_results)

  csv_path <- write_results_csv(results, output_directory)
  plot_path <- write_structure_plot(results, output_directory)
  report_path <- write_markdown_report(results, output_directory, source_mode)

  print_console_summary(results)
  cat("\nFiles created:\n")
  cat(sprintf("- %s\n", csv_path))
  cat(sprintf("- %s\n", plot_path))
  cat(sprintf("- %s\n", report_path))

  invisible(results)
}

# When sourced by the real-data test suite, expose the functions without making
# an internet request or generating project outputs automatically.
if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(error) {
      message("Analysis stopped: ", conditionMessage(error))
      quit(status = 1L)
    }
  )
}
