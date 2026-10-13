# Project 08: validated diploid counts and population-structure diagnostics.

parse_dosage <- function(gt) {
  # Phase is preserved in the VCF but irrelevant to single-locus genotype counts.
  field <- sub(":.*$", "", gt)
  missing <- field %in% c(".", "./.", ".|.", "0/.", "./0", "1/.", "./1",
                           "0|.", ".|0", "1|.", ".|1")
  valid <- grepl("^[01][/|][01]$", field)
  if (any(!valid & !missing)) stop("Non-biallelic or non-diploid GT encountered")
  dosage <- rep(NA_integer_, length(field))
  dosage[valid] <- as.integer(substr(field[valid], 1, 1)) +
                   as.integer(substr(field[valid], 3, 3))
  dosage
}

read_inputs <- function(vcf_path, panel_path) {
  lines <- readLines(vcf_path, warn = FALSE)
  reference <- lines[startsWith(lines, "##reference=")]
  if (length(reference) != 1L || !grepl("hs37d5.fa.gz", reference, fixed = TRUE) ||
      !"##fileDate=20150218" %in% lines ||
      !"##source=1000GenomesPhase3Pipeline" %in% lines)
    stop("Unexpected reference assembly or historical call-set metadata")
  hdr <- lines[startsWith(lines, "#CHROM\t")]
  rows <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (length(hdr) != 1L || length(rows) != 1L) stop("Expected one header and one variant")
  fields <- strsplit(rows, "\t", fixed = TRUE)[[1]]
  header <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
  if (length(fields) != length(header)) stop("VCF field count differs from header")
  if (!identical(fields[1:5], c("2", "136608646", "rs4988235", "G", "A")))
    stop("Wrong variant, build, position, or allele orientation")
  if (fields[7] != "PASS" || fields[9] != "GT") stop("Expected PASS with GT-only calls")
  samples <- header[-(1:9)]
  if (anyDuplicated(samples)) stop("Duplicate VCF sample IDs")
  # The original panel header has two empty trailing fields; preserve the file
  # unchanged and read only its four populated fields with fill = TRUE.
  panel <- read.table(panel_path, header = TRUE, sep = "\t", fill = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
  panel <- panel[, c("sample", "pop", "super_pop", "gender")]
  if (anyDuplicated(panel$sample)) stop("Duplicate panel sample IDs")
  if (!setequal(samples, panel$sample)) stop("VCF and panel sample sets differ")
  panel <- panel[match(samples, panel$sample), ]
  if (nrow(panel) != 2504L || length(unique(panel$pop)) != 26L ||
      length(unique(panel$super_pop)) != 5L) stop("Unexpected Phase 3 panel dimensions")
  if (anyNA(panel) || any(panel$pop == "") || any(panel$super_pop == ""))
    stop("Missing population assignment")
  panel$dosage <- parse_dosage(fields[-(1:9)])
  info <- strsplit(fields[8], ";", fixed = TRUE)[[1]]
  value <- function(key) as.numeric(sub(paste0("^", key, "="), "",
                          info[startsWith(info, paste0(key, "="))]))
  if (length(value("AC")) != 1L || length(value("AN")) != 1L ||
      sum(panel$dosage, na.rm = TRUE) != value("AC") ||
      2 * sum(!is.na(panel$dosage)) != value("AN")) stop("GT totals disagree with VCF AC/AN")
  if (abs(mean(panel$dosage, na.rm = TRUE) / 2 - value("AF")) > 5e-6)
    stop("GT frequency disagrees with rounded source AF")
  panel
}

hwe_exact <- function(n0, n1, n2) {
  # Conditional enumeration of all genotype configurations with fixed allele count.
  # Equivalent probability-ordered two-sided test to Wigginton et al. (2005).
  counts <- c(n0, n1, n2)
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)))
    stop("Genotype counts must be non-negative integers")
  n <- sum(counts)
  if (n == 0) return(NA_real_)
  m <- n1 + 2 * n2
  h <- seq(m %% 2, min(m, 2 * n - m), by = 2)
  a <- (m - h) / 2
  g <- n - h - a
  logw <- h * log(2) - lgamma(g + 1) - lgamma(h + 1) - lgamma(a + 1)
  prob <- exp(logw - max(logw))
  prob <- prob / sum(prob)
  observed <- match(n1, h)
  min(1, sum(prob[prob <= prob[observed] * (1 + 1e-10)]))
}

wilson_interval <- function(successes, n, confidence = 0.95) {
  # Applied to independent individuals' carrier status, not pooled chromosomes.
  if (n <= 0 || successes < 0 || successes > n) stop("Invalid binomial counts")
  z <- qnorm(1 - (1 - confidence) / 2)
  phat <- successes / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  c(lower = if (successes == 0) 0 else max(0, center - half),
    upper = if (successes == n) 1 else min(1, center + half))
}

summarize_group <- function(dosage, group, super_population, level) {
  called <- dosage[!is.na(dosage)]
  counts <- tabulate(called + 1L, nbins = 3L)
  n <- length(called)
  if (!n) stop("No called genotypes in group")
  p <- mean(called) / 2
  he <- 2 * p * (1 - p)
  ho <- counts[2] / n
  ci <- wilson_interval(sum(called > 0), n)
  data.frame(level = level, group = group, super_population = super_population,
    n_total = length(dosage), n_called = n, n_missing = sum(is.na(dosage)),
    n_GG = counts[1], n_GA = counts[2], n_AA = counts[3],
    A_count = sum(called), allele_number = 2 * n, A_frequency = p,
    carriers = sum(called > 0), carrier_fraction = mean(called > 0),
    carrier_wilson_lower = ci[1], carrier_wilson_upper = ci[2],
    observed_heterozygosity = ho, expected_heterozygosity = he,
    heterozygote_deficit_F = if (he > 0) 1 - ho / he else NA_real_,
    hwe_exact_p = hwe_exact(counts[1], counts[2], counts[3]),
    hwe_informative = he > 0, stringsAsFactors = FALSE, row.names = NULL)
}

population_summaries <- function(data) {
  by_pop <- lapply(sort(unique(data$pop)), function(pop) {
    x <- data[data$pop == pop, ]
    if (length(unique(x$super_pop)) != 1L) stop("Population maps to multiple super-populations")
    summarize_group(x$dosage, pop, x$super_pop[1], "population")
  })
  populations <- do.call(rbind, by_pop)
  populations$hwe_BH_q <- NA_real_
  eligible <- populations$hwe_informative
  populations$hwe_BH_q[eligible] <- p.adjust(populations$hwe_exact_p[eligible], "BH")
  super <- do.call(rbind, lapply(sort(unique(data$super_pop)), function(sp) {
    summarize_group(data$dosage[data$super_pop == sp], sp, sp, "super_population")
  }))
  pooled <- summarize_group(data$dosage, "ALL", "ALL", "pooled")
  list(populations = populations, super = super, pooled = pooled)
}

mixture_decomposition <- function(groups, pooled) {
  w <- groups$n_called / sum(groups$n_called)
  pbar <- sum(w * groups$A_frequency)
  if (abs(pbar - pooled$A_frequency) > 1e-12) stop("Partition frequency reconciliation failed")
  h_pooled <- 2 * pbar * (1 - pbar)
  h_within <- sum(w * groups$expected_heterozygosity)
  h_observed <- pooled$observed_heterozygosity
  variance <- sum(w * (groups$A_frequency - pbar)^2)
  deficit <- h_pooled - h_observed
  between <- h_pooled - h_within
  residual <- h_within - h_observed
  if (abs(between - 2 * variance) > 1e-12 ||
      abs(deficit - between - residual) > 1e-12) stop("Mixture identity failed")
  data.frame(partition = groups$level[1], n_groups = nrow(groups),
    pooled_HWE_heterozygosity = h_pooled,
    weighted_within_group_HWE_heterozygosity = h_within,
    observed_heterozygosity = h_observed,
    total_heterozygote_deficit = deficit,
    between_group_component = between,
    within_group_residual = residual,
    weighted_frequency_variance = variance,
    between_component_fraction_of_deficit = if (deficit > 0) between / deficit else NA_real_)
}

write_table <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", eol = "\n")
}
