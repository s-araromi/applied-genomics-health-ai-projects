#!/usr/bin/env Rscript
# Tests use the authentic committed data. Altered copies are rejection checks.
source("R/genotype_functions.R")
data <- read_inputs("data/rs4988235.phase3.vcf", "data/phase3_samples.panel")
s <- population_summaries(data)
passed <- 0L
test <- function(name, code) {
  force(code)
  passed <<- passed + 1L
  cat(sprintf("ok %02d - %s\n", passed, name))
}
near <- function(x, y, tolerance = 1e-12) {
  stopifnot(length(x) == length(y), all(abs(x - y) < tolerance))
}
rejects <- function(code) {
  rejected <- inherits(tryCatch({force(code); NULL}, error = identity), "error")
  stopifnot(rejected)
}

test("authentic pinned file identities", {
  manifest <- read.csv("data/source_manifest.csv")
  stopifnot(all(unname(tools::md5sum(manifest$file)) == manifest$md5),
            all(file.info(manifest$file)$size == manifest$bytes))
})
test("complete sample and source-label coverage", {
  stopifnot(nrow(data) == 2504L, !anyDuplicated(data$sample),
            length(unique(data$pop)) == 26L, length(unique(data$super_pop)) == 5L,
            !anyNA(data$dosage))
})
test("observed genotype and source allele counts", {
  stopifnot(identical(tabulate(data$dosage + 1L, 3L), c(1893L, 414L, 197L)),
            sum(data$dosage) == 808L, s$pooled$allele_number == 5008L)
  near(s$pooled$A_frequency, 808 / 5008)
})
test("source super-population frequencies", {
  stopifnot(identical(s$super$group, c("AFR", "AMR", "EAS", "EUR", "SAS")))
  near(s$super$A_frequency, c(36/1322, 150/694, 0, 511/1006, 111/978))
})
test("population partitions reconcile all genotype classes", {
  for (column in c("n_called", "n_GG", "n_GA", "n_AA", "A_count", "allele_number"))
    stopifnot(sum(s$populations[[column]]) == s$pooled[[column]])
})
test("exact test agrees with exhaustive allele-position enumeration", {
  # Independently count every assignment of the observed A copies among the
  # chromosomes of six real participants. These are mathematical null states,
  # not fabricated study records or analysis input.
  d <- data$dosage[1:6]
  n <- length(d); m <- sum(d)
  assignments <- combn(2 * n, m)
  heterozygotes <- apply(assignments, 2, function(a) {
    copies <- tabulate((a + 1L) %/% 2L, nbins = n)
    sum(copies == 1L)
  })
  probabilities <- table(heterozygotes) / ncol(assignments)
  observed <- probabilities[as.character(sum(d == 1L))]
  exhaustive_p <- sum(probabilities[probabilities <= observed * (1 + 1e-10)])
  counts <- tabulate(d + 1L, nbins = 3L)
  near(hwe_exact(counts[1], counts[2], counts[3]), unname(exhaustive_p))
})
test("exact test is invariant to allele relabeling", {
  for (i in seq_len(nrow(s$populations))) {
    g <- s$populations[i, ]
    near(hwe_exact(g$n_GG, g$n_GA, g$n_AA), hwe_exact(g$n_AA, g$n_GA, g$n_GG))
  }
})
test("monomorphic populations are flagged rather than informative tests", {
  z <- s$populations[!s$populations$hwe_informative, ]
  stopifnot(nrow(z) == 9L, all(z$hwe_exact_p == 1),
            all(is.na(z$hwe_BH_q)), all(is.na(z$heterozygote_deficit_F)))
})
test("BH correction uses only the 17 informative population tests", {
  informative <- s$populations$hwe_informative
  stopifnot(sum(informative) == 17L)
  near(s$populations$hwe_BH_q[informative],
       p.adjust(s$populations$hwe_exact_p[informative], "BH"))
  stopifnot(!any(s$populations$hwe_BH_q < 0.05, na.rm = TRUE))
})
test("heterozygosity mixture identity and recorded effect size", {
  x <- mixture_decomposition(s$populations, s$pooled)
  near(x$between_group_component, 2 * x$weighted_frequency_variance)
  near(x$total_heterozygote_deficit, x$between_group_component + x$within_group_residual)
  near(x$between_component_fraction_of_deficit, 0.9364266364103802)
})
test("coarser source groups retain additional frequency heterogeneity", {
  fine <- mixture_decomposition(s$populations, s$pooled)
  coarse <- mixture_decomposition(s$super, s$pooled)
  stopifnot(fine$between_group_component > coarse$between_group_component,
            fine$within_group_residual < coarse$within_group_residual)
})
test("carrier interval boundaries remain meaningful at zero observations", {
  d <- data$dosage[data$super_pop == "EAS"]
  ci <- wilson_interval(sum(d > 0), length(d))
  stopifnot(ci[1] == 0, ci[2] > 0, ci[2] < 0.01)
  rejects(wilson_interval(length(d) + 1, length(d)))
})
test("invalid genotype formats and exact-test counts are rejected", {
  rejects(parse_dosage("0|2")); rejects(parse_dosage("1"))
  rejects(hwe_exact(-1, 1, 1)); rejects(hwe_exact(1.5, 1, 1))
  stopifnot(is.na(hwe_exact(0, 0, 0)))
})
test("missing calls do not contribute chromosomes or carriers", {
  stopifnot(all(is.na(parse_dosage(c(".", "./.", "0|.", "./1")))))
})
test("a changed genomic build or inconsistent source allele count fails", {
  authentic <- readLines("data/rs4988235.phase3.vcf")
  temp <- tempfile(fileext = ".vcf")
  changed <- sub("136608646", "135851076", authentic, fixed = TRUE)
  writeLines(changed, temp)
  rejects(read_inputs(temp, "data/phase3_samples.panel"))
  changed <- sub("AC=808;", "AC=809;", authentic, fixed = TRUE)
  writeLines(changed, temp)
  rejects(read_inputs(temp, "data/phase3_samples.panel"))
  changed <- sub("hs37d5.fa.gz", "GRCh38.fa.gz", authentic, fixed = TRUE)
  writeLines(changed, temp)
  rejects(read_inputs(temp, "data/phase3_samples.panel"))
  unlink(temp)
})
test("sample metadata order does not change genotype assignment", {
  lines <- readLines("data/phase3_samples.panel")
  temp <- tempfile(fileext = ".panel")
  writeLines(c(lines[1], rev(lines[-1])), temp)
  reordered <- read_inputs("data/rs4988235.phase3.vcf", temp)
  stopifnot(identical(reordered$sample, data$sample),
            identical(reordered$dosage, data$dosage), identical(reordered$pop, data$pop))
  unlink(temp)
})
cat(sprintf("\n%d tests passed; authentic source fixture verified.\n", passed))
