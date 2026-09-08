#!/usr/bin/env Rscript
# propeller half of the method comparison (Supplementary Table S7).
#
#   Rscript scripts/propeller.R
#
# Reads the donor x cell-type counts written by 05_statistics.py --part methods
# and applies propeller (speckle) with the arcsine square-root transform and a
# robust moderated t-test, i.e. the tool's own defaults for a two-group design.
# Each cohort is tested separately because the two differ in tissue composition
# and in sequencing chemistry.

suppressPackageStartupMessages({
  library(speckle)
  library(limma)
})

root <- "."
counts_file <- file.path(root, "data", "interim", "donor_celltype_counts.csv")
out_file <- file.path(root, "results", "tables", "method_comparison_propeller.csv")

d <- read.csv(counts_file, check.names = FALSE)
states <- c("B", "CD4T", "CD8T", "Macrophage", "Mast", "Mono", "NK", "Plasma", "Treg", "pDC")

out <- list()
for (ds in unique(d$dataset)) {
  s <- d[d$dataset == ds, ]
  if (length(unique(s$condition)) < 2) next
  if (min(table(s$condition)) < 3) next
  # propeller wants one row per cell; expand the counts back out
  long <- do.call(rbind, lapply(states, function(g)
    data.frame(clusters = g, sample = s$gsm, group = s$condition, n = s[[g]])))
  long <- long[rep(seq_len(nrow(long)), long$n), c("clusters", "sample", "group")]
  r <- propeller(clusters = long$clusters, sample = long$sample, group = long$group,
                 transform = "asin", trend = FALSE, robust = TRUE)
  r$dataset <- ds
  r$CellType <- rownames(r)
  out[[ds]] <- r
}

res <- do.call(rbind, out)
keep <- intersect(c("dataset", "CellType", "PropMean.Control", "PropMean.Endometriosis",
                    "Tstatistic", "P.Value", "FDR"), colnames(res))
res <- res[order(res$P.Value), keep]
print(res, row.names = FALSE, digits = 3)
cat("\nsignificant at FDR < 0.05:", sum(res$FDR < 0.05), "of", nrow(res), "\n")
write.csv(res, out_file, row.names = FALSE)
cat("written:", out_file, "\n")
