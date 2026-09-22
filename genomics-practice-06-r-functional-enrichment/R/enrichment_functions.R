# Reusable functions for Project 06 functional enrichment analysis.

SOURCE_FILES <- list(
  go_annotations = list(
    filename = "goa_human.gaf.gz",
    url = "https://current.geneontology.org/annotations/goa_human.gaf.gz",
    bytes = 15041996,
    md5 = "3d16e91814f541afc3cdb03153ae0f55",
    sha256 = "db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30"
  ),
  go_ontology = list(
    filename = "go-basic.obo",
    url = "https://current.geneontology.org/ontology/go-basic.obo",
    bytes = 32227785,
    md5 = "26de643e50af9a1588f9a31062348cc4",
    sha256 = "b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc"
  ),
  reactome_mapping = list(
    filename = "NCBI2Reactome_All_Levels.txt",
    url = "https://reactome.org/download/97/NCBI2Reactome_All_Levels.txt",
    bytes = 98024930,
    md5 = "8ec3d8afe228764cf3111dc2eac6e838",
    sha256 = "53bbe42b357655e341630898547af8d23badab592e6050c34d4317fbffc04ce4"
  )
)

DERIVED_INPUTS <- list(
  universe = list(
    filename = "project05_analysis_universe.csv",
    rows = 20848L,
    md5 = "31e0a187d35d6841a81ee40bbfb030cf",
    sha256 = "75eb4546bb2d11a89dfbd2015b2764dea85109a14cc7f51a411c1a5e2176fa11"
  ),
  robust_genes = list(
    filename = "project05_robust_gene_set.csv",
    rows = 1328L,
    md5 = "8f4372dcc9e864066345a92445153b4c",
    sha256 = "2ff57aa27a8c74516b055cfd725ceb04bd1c38f2dddc52cca9a3ca2dfb439882"
  )
)

EXPERIMENTAL_GO_EVIDENCE <- c(
  "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
  "HTP", "HDA", "HMP", "HGI", "HEP"
)


assert_condition <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}


sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE)
  } else if (nzchar(sha256sum)) {
    output <- system2(sha256sum, shQuote(path), stdout = TRUE)
  } else {
    stop("Neither shasum nor sha256sum is available for source validation.", call. = FALSE)
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}


validate_file <- function(path, specification, validate_size = TRUE) {
  assert_condition(file.exists(path), paste("Required file not found:", path))
  if (validate_size) {
    observed_size <- unname(file.info(path)$size)
    assert_condition(
      identical(as.numeric(observed_size), as.numeric(specification$bytes)),
      sprintf(
        "Unexpected byte count for %s: %s observed; %s expected.",
        basename(path), observed_size, specification$bytes
      )
    )
  }
  observed_md5 <- unname(tools::md5sum(path))
  assert_condition(
    identical(observed_md5, specification$md5),
    sprintf("Unexpected MD5 for %s: %s", basename(path), observed_md5)
  )
  observed_sha256 <- sha256_file(path)
  assert_condition(
    identical(observed_sha256, specification$sha256),
    sprintf("Unexpected SHA-256 for %s: %s", basename(path), observed_sha256)
  )
  invisible(TRUE)
}


download_and_validate <- function(path, specification) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- paste0(path, ".download")
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  download.file(
    specification$url,
    temporary_path,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  validate_file(temporary_path, specification)
  installed <- file.rename(temporary_path, path)
  assert_condition(installed, paste("Could not install validated source file:", path))
  invisible(path)
}


ensure_sources <- function(data_dir, refresh = FALSE) {
  paths <- list()
  for (source_name in names(SOURCE_FILES)) {
    specification <- SOURCE_FILES[[source_name]]
    path <- file.path(data_dir, specification$filename)
    if (refresh || !file.exists(path)) {
      message(sprintf("Downloading %s (%s)...", source_name, specification$filename))
      download_and_validate(path, specification)
    }
    validate_file(path, specification)
    paths[[source_name]] <- path
  }
  paths
}


read_project05_inputs <- function(data_dir) {
  paths <- lapply(DERIVED_INPUTS, function(specification) {
    file.path(data_dir, specification$filename)
  })
  for (input_name in names(paths)) {
    specification <- DERIVED_INPUTS[[input_name]]
    validate_file(paths[[input_name]], specification, validate_size = FALSE)
  }

  universe <- read.csv(paths$universe, stringsAsFactors = FALSE, colClasses = "character")
  genes <- read.csv(paths$robust_genes, stringsAsFactors = FALSE, colClasses = "character")
  numeric_columns <- c(
    "mean_log2_fold_change", "fdr_bh",
    "sensitivity_log2_fold_change", "sensitivity_fdr_bh"
  )
  genes[numeric_columns] <- lapply(genes[numeric_columns], as.numeric)

  assert_condition(nrow(universe) == DERIVED_INPUTS$universe$rows, "Unexpected universe row count.")
  assert_condition(nrow(genes) == DERIVED_INPUTS$robust_genes$rows, "Unexpected robust-gene row count.")
  assert_condition(!anyDuplicated(universe$gene_id), "Universe gene IDs must be unique.")
  assert_condition(!anyDuplicated(universe$gene_symbol), "Universe gene symbols must be unique.")
  assert_condition(!anyDuplicated(genes$gene_id), "Robust gene IDs must be unique.")
  assert_condition(all(genes$gene_id %in% universe$gene_id), "A robust gene is absent from the universe.")
  assert_condition(
    identical(
      as.integer(table(factor(genes$direction, levels = c("Higher in tumor", "Lower in tumor")))),
      c(430L, 898L)
    ),
    "Unexpected Project 05 direction counts."
  )

  list(universe = universe, genes = genes, paths = paths)
}


read_go_terms <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  term_starts <- which(lines == "[Term]")
  rows <- vector("list", length(term_starts))

  for (index in seq_along(term_starts)) {
    start <- term_starts[[index]] + 1L
    next_stanza <- if (index < length(term_starts)) term_starts[[index + 1L]] - 1L else length(lines)
    block <- lines[start:next_stanza]
    id_line <- grep("^id: GO:", block, value = TRUE)
    name_line <- grep("^name: ", block, value = TRUE)
    namespace_line <- grep("^namespace: ", block, value = TRUE)
    is_a_lines <- grep("^is_a: GO:", block, value = TRUE)
    part_of_lines <- grep("^relationship: part_of GO:", block, value = TRUE)
    obsolete <- any(block == "is_obsolete: true")
    if (length(id_line) == 1L && length(name_line) == 1L && length(namespace_line) == 1L && !obsolete) {
      is_a_parents <- sub("^is_a: (GO:[0-9]+).*$", "\\1", is_a_lines)
      part_of_parents <- sub(
        "^relationship: part_of (GO:[0-9]+).*$",
        "\\1",
        part_of_lines
      )
      rows[[index]] <- data.frame(
        term_id = sub("^id: ", "", id_line),
        term_name = sub("^name: ", "", name_line),
        namespace = sub("^namespace: ", "", namespace_line),
        parent_ids = paste(unique(c(is_a_parents, part_of_parents)), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  assert_condition(!anyDuplicated(result$term_id), "GO term identifiers are not unique.")
  result
}


read_go_annotations <- function(path, go_terms) {
  column_classes <- c(
    "NULL", "NULL", "character", "character", "character", "NULL",
    "character", "NULL", "character", "NULL", "NULL", "NULL",
    "character", "NULL", "NULL", "NULL", "NULL"
  )
  annotations <- read.delim(
    gzfile(path),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "!",
    colClasses = column_classes,
    stringsAsFactors = FALSE
  )
  names(annotations) <- c(
    "gene_symbol", "qualifier", "term_id", "evidence_code", "aspect", "taxon"
  )
  annotations <- annotations[
    annotations$aspect == "P"
    & annotations$taxon == "taxon:9606"
    & !grepl("(^|\\|)NOT($|\\|)", annotations$qualifier),
    ,
    drop = FALSE
  ]
  annotations <- unique(annotations[c("gene_symbol", "term_id", "evidence_code")])
  valid_bp_terms <- go_terms$term_id[go_terms$namespace == "biological_process"]
  annotations <- annotations[annotations$term_id %in% valid_bp_terms, , drop = FALSE]
  assert_condition(nrow(annotations) > 0L, "No usable human GO biological-process annotations were found.")
  annotations
}


propagate_go_annotations <- function(annotations, go_terms) {
  parent_map <- setNames(strsplit(go_terms$parent_ids, ";", fixed = TRUE), go_terms$term_id)
  parent_map <- lapply(parent_map, function(values) values[nzchar(values)])
  cache <- new.env(parent = emptyenv())

  resolve_ancestors <- function(term_id, trail = character()) {
    if (exists(term_id, envir = cache, inherits = FALSE)) {
      return(get(term_id, envir = cache, inherits = FALSE))
    }
    assert_condition(!term_id %in% trail, paste("Cycle detected in GO hierarchy at", term_id))
    parents <- parent_map[[term_id]]
    parents <- parents[parents %in% go_terms$term_id]
    inherited <- if (length(parents) == 0L) {
      character()
    } else {
      unique(unlist(
        lapply(
          parents,
          resolve_ancestors,
          trail = c(trail, term_id)
        ),
        use.names = FALSE
      ))
    }
    result <- unique(c(term_id, parents, inherited))
    assign(term_id, result, envir = cache)
    result
  }

  direct_terms <- sort(unique(annotations$term_id))
  propagation_map <- do.call(rbind, lapply(direct_terms, function(term_id) {
    data.frame(
      direct_term_id = term_id,
      propagated_term_id = resolve_ancestors(term_id),
      stringsAsFactors = FALSE
    )
  }))
  expanded <- merge(
    annotations,
    propagation_map,
    by.x = "term_id",
    by.y = "direct_term_id",
    all = FALSE,
    sort = FALSE
  )
  expanded$term_id <- expanded$propagated_term_id
  unique(expanded[c("gene_symbol", "term_id", "evidence_code")])
}


read_reactome_annotations <- function(path) {
  mappings <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    colClasses = c("character", "character", "NULL", "character", "character", "character"),
    stringsAsFactors = FALSE
  )
  names(mappings) <- c("gene_id", "term_id", "term_name", "evidence_code", "species")
  mappings <- unique(mappings[mappings$species == "Homo sapiens", 1:4, drop = FALSE])
  assert_condition(nrow(mappings) > 0L, "No human Reactome mappings were found.")
  mappings
}


run_over_representation <- function(
    input_genes,
    universe_genes,
    mappings,
    term_metadata,
    gene_column,
    database,
    direction,
    label_lookup,
    minimum_term_size = 10L,
    maximum_term_size = 500L,
    minimum_overlap = 3L
) {
  mapping <- unique(mappings[c(gene_column, "term_id")])
  names(mapping)[[1]] <- "gene"
  mapping <- mapping[mapping$gene %in% universe_genes, , drop = FALSE]
  annotated_universe <- sort(unique(mapping$gene))
  tested_input <- sort(intersect(input_genes, annotated_universe))
  term_members <- split(mapping$gene, mapping$term_id)
  term_members <- lapply(term_members, unique)
  term_sizes <- vapply(term_members, length, integer(1))
  term_members <- term_members[
    term_sizes >= minimum_term_size & term_sizes <= maximum_term_size
  ]

  background_size <- length(annotated_universe)
  input_size <- length(tested_input)
  rows <- lapply(names(term_members), function(term_id) {
    members <- term_members[[term_id]]
    overlap <- intersect(tested_input, members)
    if (length(overlap) < minimum_overlap) return(NULL)
    term_size <- length(members)
    expected <- input_size * term_size / background_size
    p_value <- phyper(
      length(overlap) - 1L,
      term_size,
      background_size - term_size,
      input_size,
      lower.tail = FALSE
    )
    labels <- unname(label_lookup[overlap])
    labels[is.na(labels) | labels == ""] <- overlap[is.na(labels) | labels == ""]
    data.frame(
      database = database,
      direction = direction,
      term_id = term_id,
      background_size = background_size,
      input_size = input_size,
      term_size = term_size,
      overlap_count = length(overlap),
      expected_overlap = expected,
      enrichment_ratio = (length(overlap) / input_size) / (term_size / background_size),
      p_value = p_value,
      overlapping_gene_keys = paste(sort(overlap), collapse = ";"),
      overlapping_gene_symbols = paste(sort(unique(labels)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  assert_condition(!is.null(results) && nrow(results) > 0L, paste("No terms passed overlap filtering for", database, direction))
  results$fdr_bh <- p.adjust(results$p_value, method = "BH")
  results <- merge(results, term_metadata, by = "term_id", all.x = TRUE, sort = FALSE)
  results <- results[order(results$fdr_bh, -results$enrichment_ratio, results$term_id), ]
  rownames(results) <- NULL
  attr(results, "term_members") <- term_members
  attr(results, "annotated_universe_size") <- background_size
  attr(results, "tested_input_size") <- input_size
  results
}


attach_go_sensitivity <- function(primary_results, experimental_results) {
  sensitivity <- experimental_results[c(
    "term_id", "input_size", "term_size", "overlap_count",
    "enrichment_ratio", "p_value", "fdr_bh"
  )]
  names(sensitivity)[-1] <- paste0("experimental_", names(sensitivity)[-1])
  merged <- merge(primary_results, sensitivity, by = "term_id", all.x = TRUE, sort = FALSE)
  merged$experimental_supported <- !is.na(merged$experimental_fdr_bh) & merged$experimental_fdr_bh <= 0.05
  merged <- merged[order(merged$fdr_bh, -merged$enrichment_ratio, merged$term_id), ]
  rownames(merged) <- NULL
  merged
}


bind_rows_by_name <- function(frames) {
  all_columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  aligned <- lapply(frames, function(frame) {
    missing_columns <- setdiff(all_columns, names(frame))
    for (column in missing_columns) frame[[column]] <- NA
    frame[all_columns]
  })
  do.call(rbind, aligned)
}


jaccard_similarity <- function(first, second) {
  union_size <- length(union(first, second))
  if (union_size == 0L) return(0)
  length(intersect(first, second)) / union_size
}


select_representative_terms <- function(results, maximum_terms = 12L, threshold = 0.5) {
  candidates <- results[results$fdr_bh <= 0.05, , drop = FALSE]
  if (nrow(candidates) == 0L) return(candidates)
  selected <- integer(0)
  selected_genes <- list()
  for (index in seq_len(nrow(candidates))) {
    genes <- strsplit(candidates$overlapping_gene_symbols[[index]], ";", fixed = TRUE)[[1]]
    redundant <- any(vapply(
      selected_genes,
      function(existing) jaccard_similarity(genes, existing) >= threshold,
      logical(1)
    ))
    if (!redundant) {
      selected <- c(selected, index)
      selected_genes[[length(selected_genes) + 1L]] <- genes
    }
    if (length(selected) >= maximum_terms) break
  }
  representatives <- candidates[selected, , drop = FALSE]
  representatives$representative_rank <- seq_len(nrow(representatives))
  representatives$redundancy_threshold <- threshold
  representatives
}


write_source_manifest <- function(paths, derived_paths, output_path) {
  timestamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  source_rows <- lapply(names(paths), function(source_name) {
    specification <- SOURCE_FILES[[source_name]]
    data.frame(
      source_name = source_name,
      record_or_release = switch(
        source_name,
        go_annotations = "GOA human GAF generated 2026-05-20; declared GO version 2026-04-27",
        go_ontology = "GO release 2026-07-26",
        reactome_mapping = "Reactome release 97, released 2026-06-23"
      ),
      source_url = specification$url,
      file_name = basename(paths[[source_name]]),
      file_bytes = file.info(paths[[source_name]])$size,
      md5 = unname(tools::md5sum(paths[[source_name]])),
      sha256 = sha256_file(paths[[source_name]]),
      validated_utc = timestamp,
      stringsAsFactors = FALSE
    )
  })
  derived_rows <- lapply(names(derived_paths), function(input_name) {
    specification <- DERIVED_INPUTS[[input_name]]
    data.frame(
      source_name = paste0("project05_", input_name),
      record_or_release = "Derived from GSE19804/GPL570 Project 05 verified results",
      source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804",
      file_name = basename(derived_paths[[input_name]]),
      file_bytes = file.info(derived_paths[[input_name]])$size,
      md5 = unname(tools::md5sum(derived_paths[[input_name]])),
      sha256 = sha256_file(derived_paths[[input_name]]),
      validated_utc = timestamp,
      stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, c(source_rows, derived_rows)), output_path, row.names = FALSE)
}


summarize_results <- function(all_results, representatives, input_data, output_path) {
  combinations <- unique(all_results[c("database", "direction")])
  rows <- lapply(seq_len(nrow(combinations)), function(index) {
    database <- combinations$database[[index]]
    direction <- combinations$direction[[index]]
    subset <- all_results[
      all_results$database == database & all_results$direction == direction,
      ,
      drop = FALSE
    ]
    reps <- representatives[
      representatives$database == database & representatives$direction == direction,
      ,
      drop = FALSE
    ]
    data.frame(
      database = database,
      direction = direction,
      annotated_background_genes = unique(subset$background_size),
      annotated_input_genes = unique(subset$input_size),
      tested_terms = nrow(subset),
      significant_terms = sum(subset$fdr_bh <= 0.05),
      experimental_supported_terms = if (database == "GO Biological Process") {
        sum(subset$experimental_supported, na.rm = TRUE)
      } else {
        NA_integer_
      },
      representative_terms = nrow(reps),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary$project05_universe_genes <- nrow(input_data$universe)
  summary$project05_robust_genes <- nrow(input_data$genes)
  write.csv(summary, output_path, row.names = FALSE)
  summary
}


truncate_label <- function(value, maximum = 58L) {
  ifelse(nchar(value) > maximum, paste0(substr(value, 1L, maximum - 1L), "…"), value)
}


plot_enrichment <- function(representatives, output_path, device = c("png", "svg")) {
  device <- match.arg(device)
  if (device == "png") {
    png(output_path, width = 2400, height = 1800, res = 200)
  } else {
    svg(output_path, width = 14, height = 10)
  }
  on.exit(dev.off(), add = TRUE)

  combinations <- expand.grid(
    database = c("GO Biological Process", "Reactome"),
    direction = c("Higher in tumor", "Lower in tumor"),
    stringsAsFactors = FALSE
  )
  par(mfrow = c(2, 2), mar = c(4.5, 15, 3.5, 1), oma = c(1, 0, 3, 0))
  for (index in seq_len(nrow(combinations))) {
    database <- combinations$database[[index]]
    direction <- combinations$direction[[index]]
    panel <- representatives[
      representatives$database == database & representatives$direction == direction,
      ,
      drop = FALSE
    ]
    panel <- panel[order(panel$fdr_bh, decreasing = TRUE), ]
    if (nrow(panel) == 0L) {
      plot.new()
      title(main = paste(database, "—", direction))
      text(0.5, 0.5, "No terms reached FDR ≤ 0.05")
      next
    }
    x <- -log10(pmax(panel$fdr_bh, .Machine$double.xmin))
    point_size <- 0.8 + 0.18 * sqrt(panel$overlap_count)
    color <- if (direction == "Higher in tumor") "#D55E00" else "#0072B2"
    panel_title <- paste(
      ifelse(database == "GO Biological Process", "GO BP", database),
      "—",
      direction
    )
    plot(
      x,
      seq_len(nrow(panel)),
      pch = 19,
      cex = point_size,
      col = adjustcolor(color, alpha.f = 0.78),
      yaxt = "n",
      ylab = "",
      xlab = "−log10(BH FDR)",
      main = "",
      xlim = c(0, max(x) * 1.08)
    )
    title(main = panel_title, cex.main = 0.9)
    axis(
      2,
      at = seq_len(nrow(panel)),
      labels = truncate_label(panel$term_name, maximum = 46L),
      las = 2,
      cex.axis = 0.62
    )
    abline(v = -log10(0.05), lty = 2, col = "#555555")
  }
  mtext(
    "GSE19804 robust differential-expression functional enrichment",
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.35,
    font = 2
  )
}


write_report <- function(summary, representatives, output_path) {
  lines <- c(
    "# Project 06 Functional-Enrichment Report",
    "",
    "## Scope",
    "",
    "Over-representation analysis was performed separately for the 430 genes higher in tumor and 898 genes lower in tumor from Project 05. Backgrounds were restricted to genes measured on GPL570 and represented in each annotation source.",
    "",
    "## Result counts",
    "",
    "| Database | Direction | Annotated input | Tested terms | Significant terms | Experimentally supported GO terms |",
    "| --- | --- | ---: | ---: | ---: | ---: |"
  )
  for (index in seq_len(nrow(summary))) {
    experimental <- ifelse(
      is.na(summary$experimental_supported_terms[[index]]),
      "Not applicable",
      as.character(summary$experimental_supported_terms[[index]])
    )
    lines <- c(
      lines,
      sprintf(
        "| %s | %s | %s | %s | %s | %s |",
        summary$database[[index]],
        summary$direction[[index]],
        summary$annotated_input_genes[[index]],
        summary$tested_terms[[index]],
        summary$significant_terms[[index]],
        experimental
      )
    )
  }
  lines <- c(
    lines,
    "",
    "## Representative enriched terms",
    "",
    "| Database | Direction | Term | Overlap | Enrichment ratio | FDR |",
    "| --- | --- | --- | ---: | ---: | ---: |"
  )
  ordered <- representatives[order(representatives$database, representatives$direction, representatives$representative_rank), ]
  for (index in seq_len(nrow(ordered))) {
    lines <- c(
      lines,
      sprintf(
        "| %s | %s | %s | %s | %.2f | %.3g |",
        ordered$database[[index]],
        ordered$direction[[index]],
        ordered$term_name[[index]],
        ordered$overlap_count[[index]],
        ordered$enrichment_ratio[[index]],
        ordered$fdr_bh[[index]]
      )
    )
  }
  lines <- c(
    lines,
    "",
    "## Interpretation boundary",
    "",
    "Enrichment identifies annotation categories represented more often than expected under the declared microarray background. It does not show that every gene in a term is causal, that pathway activity was directly measured, or that the terms are independent. Bulk-tissue composition and annotation density can influence these results.",
    ""
  )
  writeLines(lines, output_path, useBytes = TRUE)
}
