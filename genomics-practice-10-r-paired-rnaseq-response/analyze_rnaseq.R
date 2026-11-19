#!/usr/bin/env Rscript
# Run from the project directory. An optional first argument selects an output
# directory, enabling a clean rerun without overwriting delivered results.
options(stringsAsFactors = FALSE, scipen = 0)
Sys.setenv(TZ = "UTC")
source("R/rnaseq_functions.R")
require_packages()
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) write_csv(x, file.path(outdir, name))
x <- read_inputs()
keep <- expression_filter(x$counts)
cts <- x$counts[keep, , drop = FALSE]; meta <- x$samples
message(sum(keep), " of ", nrow(x$counts), " genes retained")
write_out(data.frame(gene_id = rownames(x$counts), retained = keep,
                     libraries_ge_10 = rowSums(x$counts >= 10)), "gene_filter.csv")

paired <- fit_counts(cts, meta)
res <- extract_results(paired)
res$gene_name <- x$genes$gene_name[match(res$gene_id, x$genes$gene_id)]
res$gene_biotype <- x$genes$gene_biotype[match(res$gene_id, x$genes$gene_id)]
normalized <- counts(paired, normalized = TRUE)
qc <- data.frame(meta, assigned_fragments = colSums(x$counts),
                  expressed_genes = colSums(x$counts > 0),
                  retained_gene_fragments = colSums(cts),
                  size_factor = sizeFactors(paired))
write_out(qc, "sample_qc.csv")
write_out(data.frame(gene_id = rownames(cts), normalized, check.names = FALSE),
          "normalized_counts.csv")
mc <- as.data.frame(mcols(paired))
write_out(data.frame(gene_id = rownames(cts), base_mean = mc$baseMean,
                     dispersion_gene = mc$dispGeneEst, dispersion_trend = mc$dispFit,
                     dispersion_final = mc$dispersion), "dispersion_estimates.csv")
cooks <- assays(paired)[["cooks"]]
p <- ncol(model.matrix(~ donor + dex, meta)); n <- nrow(meta)
cook_cut <- qf(.99, p, n-p)
write_out(data.frame(gene_id = rownames(cts), max_cooks = apply(cooks, 1, max),
                     reference_cutoff = cook_cut,
                     point_p_available = !is.na(res$p_point_null)), "influence_diagnostics.csv")

# VST is used only for sample exploration, never as input to count models.
vst <- assay(varianceStabilizingTransformation(paired, blind = FALSE))
pc <- pca_counts(vst)
write_out(data.frame(run = rownames(pc$scores), pc$scores, check.names = FALSE), "pca_scores.csv")
write_out(data.frame(component = seq_along(pc$variance), variance_fraction = pc$variance),
          "pca_variance.csv")
write_out(data.frame(gene_id = pc$genes), "pca_features.csv")
write_out(data.frame(run = rownames(meta), as.matrix(dist(t(vst))), check.names = FALSE),
          "sample_distances.csv")

message("Fitting unpaired design sensitivity")
unpaired <- extract_results(fit_counts(cts, meta, paired = FALSE))
write_out(unpaired, "unpaired_design_results.csv")
comparison <- data.frame(gene_id = res$gene_id, paired_lfc = res$log2_fold_change,
                          unpaired_lfc = unpaired$log2_fold_change,
                          paired_q = res$q_above_1_5_fold,
                          unpaired_q = unpaired$q_above_1_5_fold,
                          paired_selected = res$primary_selected,
                          unpaired_selected = unpaired$primary_selected)
write_out(comparison, "design_comparison.csv")

# Refit all model parameters after removing BOTH libraries from one donor.
# The original retained-gene family is fixed across these descriptive refits.
lodo <- lapply(levels(meta$donor), function(d) {
  message("Omitting donor ", d)
  take <- meta$donor != d
  m <- droplevels(meta[take, , drop = FALSE])
  z <- extract_results(fit_counts(cts[, take, drop = FALSE], m))
  data.frame(omitted_donor = d, z)
})
lodo <- do.call(rbind, lodo)
write_out(lodo, "leave_one_donor_out.csv")
lodo_lfc <- vapply(levels(meta$donor), function(d)
  lodo$log2_fold_change[lodo$omitted_donor == d], numeric(nrow(cts)))
lodo_selected <- vapply(levels(meta$donor), function(d)
  lodo$primary_selected[lodo$omitted_donor == d], logical(nrow(cts)))
res$lodo_same_direction <- rowSums(sign(lodo_lfc) == sign(res$log2_fold_change))
res$lodo_selected_refits <- rowSums(lodo_selected)
res$lodo_min_abs_lfc <- apply(abs(lodo_lfc), 1, min)
res$lodo_max_abs_change <- apply(abs(lodo_lfc - res$log2_fold_change), 1, max)
res$direction_stable_all_refits <- res$lodo_same_direction == nlevels(meta$donor)
write_out(res, "paired_differential_expression.csv")
selected <- res[res$primary_selected, , drop = FALSE]
selected <- selected[order(selected$q_above_1_5_fold, selected$gene_id), , drop = FALSE]
write_out(selected, "selected_response_genes.csv")
ratios <- paired_log_ratios(normalized, meta)
write_out(data.frame(gene_id = rownames(ratios), ratios, check.names = FALSE),
          "donor_log2_ratios.csv")

summary <- data.frame(metric = c("source_gene_rows", "samples", "donors", "retained_genes",
 "point_null_fdr_below_0_05", "primary_1_5_fold_fdr_below_0_05", "primary_up", "primary_down",
 "primary_direction_stable_all_refits", "primary_selected_in_all_refits",
 "unpaired_primary_selected", "paired_unpaired_overlap", "primary_missing_p",
 "pca_pc1_variance_fraction", "pca_pc2_variance_fraction", "lfc_spearman_paired_unpaired"),
 value = c(nrow(x$counts), ncol(cts), nlevels(meta$donor), nrow(cts),
 sum(res$q_point_null < .05, na.rm = TRUE), nrow(selected),
 sum(selected$log2_fold_change > 0), sum(selected$log2_fold_change < 0),
 sum(selected$direction_stable_all_refits), sum(selected$lodo_selected_refits == 4),
 sum(unpaired$primary_selected), sum(res$primary_selected & unpaired$primary_selected),
 sum(is.na(res$p_above_1_5_fold)), pc$variance[1:2],
 cor(res$log2_fold_change, unpaired$log2_fold_change, method = "spearman")))
write_out(summary, "analysis_summary.csv")
write_out(data.frame(model = c("paired", "unpaired", paste0("omit_", levels(meta$donor))),
                    samples = c(8, 8, rep(6,4)),
                    design_columns = c(5,2,rep(4,4)), residual_df = c(3,6,rep(2,4)),
                    dispersion_fit = rep("parametric", 6)),
          "model_specification.csv")
# Actual fit types are checked in the fitting helper; do not silently switch.
source("R/plots.R")
plot_outputs(outdir, pc, meta, mc, res, ratios, selected)
writeLines(sub("[ \t]+$", "", capture.output(sessionInfo())),
           file.path(outdir, "session_info.txt"))
report <- c("# Paired RNA-seq response analysis", "",
 "Primary contrast: dexamethasone versus untreated; donor-adjusted negative-binomial model.",
 "Primary null: absolute log2 fold change <= log2(1.5); BH FDR < 0.05.", "",
 paste0("- ", summary$metric, ": ", format(summary$value, digits = 8)), "",
 "The leave-one-donor-out results measure sensitivity within this dataset, not independent validation.",
 "No clinical response, patient outcome, or predictive performance was estimated.")
writeLines(report, file.path(outdir, "analysis_report.md"))
print(summary, row.names = FALSE)
