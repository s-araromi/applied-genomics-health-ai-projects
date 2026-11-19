# Validation and statistical helpers. These functions never download data.
require_packages <- function() {
  needed <- c("DESeq2", "SummarizedExperiment", "digest")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
}

verify_manifest <- function(root = ".") {
  m <- read.csv(file.path(root, "data/source_manifest.csv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(m))) {
    path <- file.path(root, m$path[i])
    if (!file.exists(path) || digest::digest(file = path, algo = "sha256") != m$sha256[i])
      stop("Source identity mismatch: ", m$path[i])
  }
  invisible(m)
}

validate_inputs <- function(counts, samples) {
  if (!is.matrix(counts) || !is.numeric(counts) || !nrow(counts) ||
      any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)) ||
      any(counts > .Machine$integer.max)) stop("Counts must be finite non-negative integers")
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      anyNA(rownames(counts)) || any(!nzchar(rownames(counts)))) stop("Invalid gene IDs")
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts)) ||
      anyNA(colnames(counts))) stop("Invalid library IDs")
  needed <- c("run", "donor", "dex", "albuterol")
  if (all(needed %in% names(samples)))
    samples[needed] <- lapply(samples[needed], as.character)
  if (!all(needed %in% names(samples)) || anyNA(samples[, needed]) ||
      anyDuplicated(samples$run) || any(!nzchar(samples$donor)) ||
      !setequal(samples$run, colnames(counts))) stop("Sample metadata mismatch")
  samples <- samples[match(colnames(counts), samples$run), , drop = FALSE]
  if (!all(samples$albuterol == "untrt") ||
      !all(samples$dex %in% c("untrt", "trt"))) stop("Unexpected treatment")
  samples$donor <- factor(samples$donor, levels = sort(unique(samples$donor)))
  samples$dex <- factor(samples$dex, levels = c("untrt", "trt"))
  pairing <- table(samples$donor, samples$dex)
  if (nrow(pairing) < 3 || any(pairing != 1L)) stop("Require one complete pair per donor, at least three donors")
  rownames(samples) <- samples$run
  design <- model.matrix(~ donor + dex, samples)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design))
    stop("Design is not estimable")
  storage.mode(counts) <- "integer"
  list(counts = counts, samples = samples, design = design)
}

read_inputs <- function(root = ".") {
  verify_manifest(root)
  tab <- read.csv(file.path(root, "data/counts.csv"), check.names = FALSE)
  counts <- as.matrix(tab[, -1, drop = FALSE]); rownames(counts) <- tab$gene_id
  samples <- read.csv(file.path(root, "data/samples.csv"), stringsAsFactors = FALSE)
  x <- validate_inputs(counts, samples)
  genes <- read.csv(file.path(root, "data/gene_annotation.csv"), stringsAsFactors = FALSE)
  if (anyDuplicated(genes$gene_id) || !identical(genes$gene_id, rownames(x$counts)))
    stop("Gene annotation alignment mismatch")
  x$genes <- genes
  x
}

expression_filter <- function(counts, min_count = 10L, min_libraries = 4L) {
  if (min_count < 1 || min_libraries < 1 || min_libraries > ncol(counts))
    stop("Invalid expression-filter parameters")
  rowSums(counts >= min_count) >= min_libraries
}

# Count unavailable p values as non-rejections in the declared family. Keep
# their exported q values unavailable, rather than interpreting them as nulls.
bh_family <- function(p) {
  if (any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))) stop("Invalid p value")
  q <- p.adjust(replace(p, is.na(p), 1), method = "BH")
  q[is.na(p)] <- NA_real_
  q
}

fit_counts <- function(counts, samples, paired = TRUE) {
  samples$donor <- droplevels(samples$donor)
  form <- if (paired) ~ donor + dex else ~ dex
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, design = form)
  # Standard median-of-ratios normalization; no outlier count replacement.
  dds <- DESeq2::DESeq(dds, fitType = "parametric", betaPrior = FALSE,
                      minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
  if (!identical(attr(DESeq2::dispersionFunction(dds), "fitType"), "parametric"))
    stop("Dispersion fit changed from the declared parametric method")
  if (any(!is.finite(DESeq2::sizeFactors(dds)))) stop("Invalid size factors")
  dds
}

extract_results <- function(dds, threshold = log2(1.5)) {
  args <- list(object = dds, contrast = c("dex", "trt", "untrt"),
               independentFiltering = FALSE, alpha = .05)
  point <- as.data.frame(do.call(DESeq2::results, args))
  # Project 10 threshold-method compatibility correction
  # Preserve the original statistical procedure across the DESeq2 1.44 change.
  threshold_method <- if (utils::packageVersion("DESeq2") >= "1.44.0")
    "greaterAbs2014" else "greaterAbs"
  effect <- as.data.frame(do.call(DESeq2::results,
                         c(args, list(lfcThreshold = threshold, altHypothesis = threshold_method))))
  out <- data.frame(gene_id = rownames(point), base_mean = point$baseMean,
                    log2_fold_change = point$log2FoldChange, lfc_se = point$lfcSE,
                    wald_ci_lower = point$log2FoldChange - qnorm(.975) * point$lfcSE,
                    wald_ci_upper = point$log2FoldChange + qnorm(.975) * point$lfcSE,
                    p_point_null = point$pvalue, q_point_null = bh_family(point$pvalue),
                    p_above_1_5_fold = effect$pvalue,
                    q_above_1_5_fold = bh_family(effect$pvalue),
                    beta_converged = as.logical(S4Vectors::mcols(dds)$betaConv))
  # Non-converged coefficients must never be selected as findings.
  if (any(!out$beta_converged)) stop("Non-converged model coefficients; inspect before interpreting")
  out$primary_selected <- !is.na(out$q_above_1_5_fold) & out$q_above_1_5_fold < .05
  out
}

pca_counts <- function(vst, top_n = 500L) {
  vars <- apply(vst, 1, var)
  ids <- order(-vars, rownames(vst))
  ids <- ids[seq_len(min(top_n, sum(vars > 0)))]
  if (length(ids) < 2) stop("Insufficient variable genes for PCA")
  p <- prcomp(t(vst[ids, , drop = FALSE]), center = TRUE, scale. = FALSE)
  # Orient each PC so its largest-magnitude gene loading is positive.
  for (j in seq_len(ncol(p$rotation))) {
    anchor <- which.max(abs(p$rotation[, j]))
    if (p$rotation[anchor, j] < 0) {
      p$rotation[, j] <- -p$rotation[, j]; p$x[, j] <- -p$x[, j]
    }
  }
  list(scores = p$x, variance = p$sdev^2 / sum(p$sdev^2),
       genes = rownames(vst)[ids])
}

paired_log_ratios <- function(normalized, samples, pseudocount = .5) {
  donors <- levels(samples$donor)
  z <- vapply(donors, function(d) {
    trt <- which(samples$donor == d & samples$dex == "trt")
    ctl <- which(samples$donor == d & samples$dex == "untrt")
    log2((normalized[, trt] + pseudocount) / (normalized[, ctl] + pseudocount))
  }, numeric(nrow(normalized)))
  dimnames(z) <- list(rownames(normalized), donors)
  z
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
