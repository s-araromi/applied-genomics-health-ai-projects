#!/usr/bin/env Rscript
# Project 08: observed reference-panel frequencies and genotype-structure audit.
source("R/genotype_functions.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) stop("Run without arguments from the project directory")
dir.create("outputs", showWarnings = FALSE)
manifest <- read.csv("data/source_manifest.csv", stringsAsFactors = FALSE)
for (i in seq_len(nrow(manifest))) {
  path <- manifest$file[i]
  if (!file.exists(path) || unname(tools::md5sum(path)) != manifest$md5[i])
    stop(paste("Missing or changed pinned source:", path))
}
data <- read_inputs("data/rs4988235.phase3.vcf", "data/phase3_samples.panel")
summaries <- population_summaries(data)
pop <- summaries$populations
sp <- summaries$super
pooled <- summaries$pooled
mixture <- rbind(mixture_decomposition(pop, pooled), mixture_decomposition(sp, pooled))
stopifnot(sum(pop$n_called) == pooled$n_called, sum(pop$A_count) == pooled$A_count,
          sum(pop$n_GA) == pooled$n_GA)
write_table(pop, "outputs/population_genotype_summary.csv")
write_table(sp, "outputs/super_population_summary.csv")
write_table(pooled, "outputs/pooled_genotype_summary.csv")
write_table(mixture, "outputs/heterozygosity_decomposition.csv")
qc <- data.frame(metric = c("source_samples", "called_samples", "missing_genotypes",
  "populations", "super_populations", "ALT_alleles", "total_alleles",
  "informative_population_HWE_tests", "population_tests_BH_q_below_0.05"),
  value = c(nrow(data), pooled$n_called, pooled$n_missing, nrow(pop), nrow(sp),
    pooled$A_count, pooled$allele_number, sum(pop$hwe_informative),
    sum(pop$hwe_BH_q < 0.05, na.rm = TRUE)))
write_table(qc, "outputs/source_qc.csv")

colors <- c(AFR = "#278875", AMR = "#D58A22", EAS = "#688CB5", EUR = "#8B5DA9", SAS = "#CF5966")
draw <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.45, 1))
  par(mar = c(4.5, 4.5, 3.3, 1))
  ordered <- pop[order(pop$super_population, pop$A_frequency), ]
  ypos <- barplot(ordered$A_frequency, names.arg = ordered$group, horiz = TRUE,
    las = 1, col = colors[ordered$super_population], border = NA, xlim = c(0, 1),
    xlab = "GRCh37 forward-strand ALT A frequency",
    main = "rs4988235: 26 source populations", cex.names = 0.7)
  legend("bottomright", legend = names(colors), fill = colors, border = NA,
         bty = "n", cex = 0.85)
  par(mar = c(7.5, 4.5, 3.3, 1))
  values <- c(pooled$expected_heterozygosity,
              mixture$weighted_within_group_HWE_heterozygosity[1],
              pooled$observed_heterozygosity)
  at <- barplot(values, names.arg = c("Pooled HWE", "26-group HWE", "Observed"),
    col = c("#BDA9CE", "#688CB5", "#278875"), border = NA, las = 2,
    ylab = "Heterozygote proportion", ylim = c(0, max(values) * 1.2),
    main = "Population mixture and heterozygosity")
  text(at, values, sprintf("%.4f", values), pos = 3, cex = 0.9)
  mtext("Reference panel; frequencies describe sampled individuals", side = 1, line = 6, cex = 0.65)
}
png("outputs/population_structure_overview.png", width = 1900, height = 1300, res = 160)
draw(); dev.off()
svg("outputs/population_structure_overview.svg", width = 12, height = 8)
draw(); dev.off()
# Trim serializer whitespace so every generated text file passes Git checks.
svg_file <- "outputs/population_structure_overview.svg"
writeLines(sub("[[:blank:]]+$", "", readLines(svg_file)), svg_file)

fmt <- function(x) format(x, scientific = TRUE, digits = 5)
lines <- c("# Project 08: Population Genotype Structure Audit", "",
  "## Verified source and sample accounting", "",
  "1000 Genomes Phase 3 integrated v5a, release directory 20130502; GRCh37.",
  "Single-locus subset: rs4988235, chr2:136608646 G>A.",
  sprintf("All %d source samples aligned to 26 populations and five super-population labels.", nrow(data)),
  sprintf("Called genotypes: %d; missing: %d. ALT count: %d/%d (frequency %.6f).",
          pooled$n_called, pooled$n_missing, pooled$A_count, pooled$allele_number, pooled$A_frequency), "",
  "## Genotype and equilibrium results", "",
  sprintf("Pooled GG/GA/AA counts: %d/%d/%d.", pooled$n_GG, pooled$n_GA, pooled$n_AA),
  paste("Pooled probability-ordered exact HWE p-value:", fmt(pooled$hwe_exact_p)),
  sprintf("Observed heterozygosity %.6f; pooled HWE expectation %.6f.",
          pooled$observed_heterozygosity, pooled$expected_heterozygosity),
  sprintf("Of %d polymorphic population tests, %d had BH q<0.05.",
    sum(pop$hwe_informative), sum(pop$hwe_BH_q < 0.05, na.rm = TRUE)), "",
  "## Mixture decomposition", "",
  sprintf("26-population between-group component: %.6f; within-group residual: %.6f.",
          mixture$between_group_component[1], mixture$within_group_residual[1]),
  sprintf("Between-group component / total heterozygote deficit: %.2f%%.",
          100 * mixture$between_component_fraction_of_deficit[1]), "",
  "## Interpretation and limits", "",
  "Pooling populations with different allele frequencies changes the equilibrium expectation.",
  "A small pooled HWE p-value does not establish a genotype-calling error or selection.",
  "Group labels come from the source panel and are not genetically inferred ancestry or clinical labels.",
  "Carrier intervals assume independent individuals and hypothetical binomial sampling; the reference panel is not a random prevalence survey.",
  "The allele is not a universal predictor of lactase persistence; other regulatory variants occur in Africa.",
  "No phenotype, lactose-intolerance diagnosis, dietary advice, association analysis, or selection inference was produced.")
writeLines(lines, "outputs/population_structure_report.md")
session_lines <- c(capture.output(sessionInfo()), "", "Input manifest:",
                   capture.output(print(manifest)))
writeLines(sub("[[:blank:]]+$", "", session_lines), "outputs/session_info.txt")
print(qc, row.names = FALSE)
print(pooled)
print(mixture)
