#!/usr/bin/env Rscript
# Tests use authentic committed data. Deliberate corruptions below test rejection
# paths; they are never analyzed or distributed as replacement datasets.
source("R/rnaseq_functions.R")
require_packages()
n <- 0L
test <- function(label, expr) {
  force(expr); n <<- n + 1L; cat(sprintf("ok %02d - %s\n", n, label))
}
rejects <- function(expr) {
  got <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  stopifnot(got)
}
x <- read_inputs(); cts <- x$counts; m <- x$samples

test("all committed source identities match SHA-256 pins", {
  manifest <- verify_manifest(); stopifnot(nrow(manifest) == 5L)
})
test("complete authentic gene and sample coverage", {
  stopifnot(identical(dim(cts), c(63677L,8L)), nrow(x$genes) == 63677,
            nlevels(m$donor) == 4L, all(table(m$donor,m$dex) == 1))
})
test("source library totals agree with the pinned airway object", {
  stopifnot(identical(as.numeric(colSums(cts)),
    c(20637971,18809481,25348649,15163415,24448408,30818215,19126151,21164133)))
})
test("GEO status snapshot is public and includes selected samples", {
  text <- readLines("data/GSE52778_status_2026-10-19.soft", warn = FALSE)
  stopifnot(any(text == "^SERIES = GSE52778"),
            any(grepl("^!Series_status = Public",text)),
            all(paste0("!Series_sample_id = ",m$geo_sample) %in% text),
            sum(grepl("^!Series_sample_id = ",text)) == 16L)
})
test("metadata reordering preserves count-to-sample assignment", {
  y <- validate_inputs(cts, m[8:1,])
  stopifnot(identical(y$samples$run,m$run), identical(y$counts,cts))
})
test("fractional, missing and negative counts fail before fitting", {
  for (value in c(.5,NA,-1)) {bad <- cts; bad[1,1] <- value; rejects(validate_inputs(bad,m))}
})
test("duplicate gene and sample identities are rejected", {
  bad <- cts; rownames(bad)[2] <- rownames(bad)[1]; rejects(validate_inputs(bad,m))
  bad <- m; bad$run[2] <- bad$run[1]; rejects(validate_inputs(cts,bad))
})
test("missing pair and albuterol contamination are rejected", {
  rejects(validate_inputs(cts[,-1],m[-1,]))
  bad <- m; bad$albuterol[1] <- "trt"; rejects(validate_inputs(cts,bad))
})
test("unrecognized treatment and sample mismatch are rejected", {
  bad <- m; bad$dex <- as.character(bad$dex); bad$dex[1] <- "other"
  rejects(validate_inputs(cts,bad))
  bad <- m; bad$run[1] <- "unexpected"; rejects(validate_inputs(cts,bad))
})
test("paired design is estimable with three residual degrees of freedom", {
  stopifnot(qr(x$design)$rank == 5L, nrow(x$design)-ncol(x$design) == 3L)
})
keep <- expression_filter(cts)
test("count filter retains the verified authentic family", {
  stopifnot(sum(keep) == 16139L, all(rowSums(cts[keep,] >= 10) >= 4),
            all(rowSums(cts[!keep,] >= 10) < 4))
})
test("gene annotation preserves stable IDs and historical source labels", {
  stopifnot(identical(x$genes$gene_id,rownames(cts)),
            x$genes$gene_name[1] == "TSPAN6",
            x$genes$gene_seq_start[1] == 99883667,
            x$genes$seq_name[1] == "X")
})
res <- read.csv("outputs/paired_differential_expression.csv")
primary <- read.csv("outputs/selected_response_genes.csv")
test("result rows reconcile exactly to the retained gene family", {
  stopifnot(identical(res$gene_id,rownames(cts)[keep]), all(res$beta_converged),
            identical(sort(primary$gene_id),sort(res$gene_id[res$primary_selected])))
})
test("threshold hypothesis is stricter than the point-null hypothesis", {
  valid <- complete.cases(res[,c("p_point_null","p_above_1_5_fold")])
  stopifnot(all(res$p_above_1_5_fold[valid] + 1e-14 >= res$p_point_null[valid]),
            all(abs(primary$log2_fold_change) > log2(1.5)),
            all(primary$q_above_1_5_fold < .05))
})
test("full-family BH adjustment includes unavailable tests as non-rejections", {
  p <- res$p_above_1_5_fold
  reference <- p.adjust(p[!is.na(p)],method="BH",n=nrow(res))
  stopifnot(isTRUE(all.equal(res$q_above_1_5_fold[!is.na(p)],reference,tolerance=1e-12)))
  if (anyNA(p)) stopifnot(all(is.na(bh_family(p)[is.na(p)])))
  rejects(bh_family(c(p,2)))
})
test("Wald confidence intervals bracket each reported effect", {
  stopifnot(all(res$wald_ci_lower <= res$log2_fold_change),
            all(res$wald_ci_upper >= res$log2_fold_change), all(res$lfc_se > 0))
})
qc <- read.csv("outputs/sample_qc.csv")
test("normalization factors and fragment totals are valid", {
  stopifnot(identical(qc$run,m$run), all(qc$size_factor > 0),
            abs(exp(mean(log(qc$size_factor)))-1) < .01,
            all(qc$retained_gene_fragments <= qc$assigned_fragments))
})
pcs <- read.csv("outputs/pca_scores.csv"); pv <- read.csv("outputs/pca_variance.csv")
test("PCA uses 500 genuine genes and preserves centered sample geometry", {
  features <- read.csv("outputs/pca_features.csv")
  stopifnot(nrow(features) == 500L, all(features$gene_id %in% res$gene_id),
            identical(pcs$run,m$run), abs(sum(pv$variance_fraction)-1)<1e-12,
            max(abs(colMeans(as.matrix(pcs[,-1])))) < 1e-8)
})
test("sample distances are symmetric with a zero diagonal", {
  ds <- read.csv("outputs/sample_distances.csv",check.names=FALSE)
  z <- as.matrix(ds[,-1]);stopifnot(all(diag(z)==0),max(abs(z-t(z)))<1e-12)
})
lodo <- read.csv("outputs/leave_one_donor_out.csv")
test("each donor omission covers the same complete testing family", {
  stopifnot(nrow(lodo) == 4L*nrow(res))
  for (d in levels(m$donor)) {
    z <- lodo[lodo$omitted_donor == d,]
    stopifnot(identical(z$gene_id,res$gene_id))
    take <- m$donor != d
    y <- validate_inputs(cts[,take],droplevels(m[take,]))
    stopifnot(nlevels(y$samples$donor)==3L,nrow(y$design)-ncol(y$design)==2L)
  }
})
test("direction and significance sensitivity agree with exported refits", {
  z <- vapply(levels(m$donor),function(d) lodo$log2_fold_change[lodo$omitted_donor==d],numeric(nrow(res)))
  s <- vapply(levels(m$donor),function(d) lodo$primary_selected[lodo$omitted_donor==d],logical(nrow(res)))
  stopifnot(all(res$lodo_same_direction == rowSums(sign(z)==sign(res$log2_fold_change))),
            all(res$lodo_selected_refits == rowSums(s)),
            all(res$direction_stable_all_refits == (res$lodo_same_direction==4)))
})
test("descriptive donor ratios reproduce normalized count contrasts", {
  tab <- read.csv("outputs/normalized_counts.csv",check.names=FALSE)
  z <- as.matrix(tab[,-1]);rownames(z)<-tab$gene_id
  actual <- read.csv("outputs/donor_log2_ratios.csv",check.names=FALSE)
  stopifnot(isTRUE(all.equal(unname(as.matrix(actual[,-1])),
                            unname(paired_log_ratios(z,m)),tolerance=1e-12)))
})
test("source corruption is detected before analysis", {
  temp <- tempfile();dir.create(file.path(temp,"data"),recursive=TRUE)
  file.copy("data/source_manifest.csv",file.path(temp,"data/source_manifest.csv"))
  rejects(verify_manifest(temp));unlink(temp,recursive=TRUE)
})
# Project 10 threshold-method compatibility correction
test("threshold p-values preserve the original method in all six fits", {
  for (tab in list(res, read.csv("outputs/unpaired_design_results.csv"), lodo)) {
    valid <- !is.na(tab$p_above_1_5_fold)
    expected <- pmin(1, 2 * pnorm(
      (abs(tab$log2_fold_change[valid]) - log2(1.5)) / tab$lfc_se[valid],
      lower.tail = FALSE))
    stopifnot(length(expected) > 0L, all(is.finite(expected)),
              max(abs(tab$p_above_1_5_fold[valid] - expected)) < 1e-12)
  }
})
cat(sprintf("\n%d tests passed; authentic source and analysis outputs verified.\n",n))
