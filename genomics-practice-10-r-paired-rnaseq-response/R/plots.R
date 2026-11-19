# Base-R graphics keep the plotting dependency surface small.
plot_outputs <- function(outdir, pc, meta, mc, res, ratios, selected) {
  navy <- "#172F4D"; teal <- "#087E8B"; orange <- "#CC5B31"; gray <- "#A9B3BF"
  sample_colors <- ifelse(meta$dex == "trt", orange, teal)
  donor_shapes <- seq_len(nlevels(meta$donor)) + 14L
  draw_overview <- function() {
    par(mfrow = c(2,2), mar = c(4.5,4.8,3.5,1.2), oma = c(1,0,2,0),
        col.axis = navy, col.lab = navy, col.main = navy, family = "sans")
    a <- pc$scores
    plot(a[,1], a[,2], type = "n", xlim = range(a[,1])*1.15,
         ylim = range(a[,2])*1.2,
         xlab = sprintf("PC1 (%.1f%%)", 100*pc$variance[1]),
         ylab = sprintf("PC2 (%.1f%%)", 100*pc$variance[2]),
         main = "A  |  VST sample structure: 500 variable genes")
    for (d in levels(meta$donor)) {
      take <- which(meta$donor == d)
      segments(a[take[1],1], a[take[1],2], a[take[2],1], a[take[2],2], col = gray)
    }
    points(a[,1], a[,2], pch = donor_shapes[as.integer(meta$donor)],
           col = sample_colors, cex = 1.5)
    legend("topright", c("Untreated", "Dexamethasone"), col = c(teal,orange), pch = 16,
           bty = "n", cex = .8)
    legend("bottomleft", levels(meta$donor), pch = donor_shapes, bty = "n", cex = .7)
    good <- is.finite(mc$dispGeneEst) & mc$dispGeneEst > 0
    plot(mc$baseMean[good], mc$dispGeneEst[good], log = "xy", pch = 16, cex = .25,
         col = adjustcolor(gray, .4), xlab = "Mean normalized fragments",
         ylab = "Negative-binomial dispersion", main = "B  |  Dispersion estimation")
    ord <- order(mc$baseMean)
    lines(mc$baseMean[ord], mc$dispFit[ord], col = orange, lwd = 2)
    legend("topright", c("Gene estimate", "Fitted trend"), col = c(gray,orange),
           pch = c(16,NA), lty = c(NA,1), bty = "n", cex = .8)
    colors <- rep(adjustcolor(gray,.45), nrow(res))
    colors[res$primary_selected & res$log2_fold_change > 0] <- orange
    colors[res$primary_selected & res$log2_fold_change < 0] <- teal
    plot(res$base_mean, res$log2_fold_change, log = "x", pch = 16, cex = .35,
         col = colors, xlab = "Mean normalized fragments", ylab = "Log2 fold change: treated / untreated",
         main = "C  |  Donor-adjusted treatment response")
    abline(h = c(-1,1)*log2(1.5), lty = 2, col = navy)
    legend("bottomright", c("Selected increase", "Selected decrease"),
           col = c(orange,teal), pch = 16, bty = "n", cex = .8)
    freq <- table(factor(selected$lodo_selected_refits, levels = 0:4))
    barplot(freq, names.arg = 0:4, col = teal, border = NA,
            xlab = "Omission refits retaining primary significance (of 4)",
            ylab = "Genes selected in the full paired model",
            main = "D  |  Leave-one-donor-out sensitivity")
    mtext("GSE52778  |  Eight libraries, four paired donors", outer = TRUE,
          side = 3, font = 2, cex = 1.25, col = navy)
    mtext("Primary test: |log2 fold change| > log2(1.5), BH FDR < 0.05. Omission refits are not external validation.",
          side = 1, outer = TRUE, cex = .8, col = navy)
  }
  png(file.path(outdir,"rnaseq_response_overview.png"), width = 2000, height = 1540, res = 170)
  draw_overview(); dev.off()
  svg(file.path(outdir,"rnaseq_response_overview.svg"), width = 12, height = 9.2)
  draw_overview(); dev.off()
  top <- head(selected, 12)
  z <- ratios[match(top$gene_id, rownames(ratios)), , drop = FALSE]
  labels <- ifelse(nzchar(top$gene_name), top$gene_name, top$gene_id)
  draw_genes <- function() {
    par(mar = c(5,8,4,2), col.axis = navy, col.lab = navy, col.main = navy)
    lim <- max(abs(z))
    image(seq_len(ncol(z)), seq_len(nrow(z)), t(z[nrow(z):1,,drop=FALSE]),
          zlim = c(-lim,lim), axes = FALSE,
          col = colorRampPalette(c(teal,"#F7F9FC",orange))(101),
          xlab = "Donor (source cell-line label)", ylab = "",
          main = "Top 12 response genes | descriptive within-donor log2 ratios")
    axis(1, at = seq_len(ncol(z)), labels = colnames(z))
    axis(2, at = seq_len(nrow(z)), labels = rev(labels), las = 2, cex.axis = .9)
    for (i in seq_len(ncol(z))) for (j in seq_len(nrow(z)))
      text(i, j, sprintf("%.2f", z[nrow(z)+1-j,i]), cex = .8,
           col = ifelse(abs(z[nrow(z)+1-j,i]) > lim*.65,"white",navy))
    mtext("log2[(normalized treated + 0.5) / (normalized untreated + 0.5)]; genes ranked by primary q value",
          side = 1, line = 3.5, cex = .72)
    box()
  }
  png(file.path(outdir,"donor_response_top_genes.png"), width = 1600, height = 1100, res = 160)
  draw_genes(); dev.off()
  svg(file.path(outdir,"donor_response_top_genes.svg"), width = 10, height = 7)
  draw_genes(); dev.off()
}
