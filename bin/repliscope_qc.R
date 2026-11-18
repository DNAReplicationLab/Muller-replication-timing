#!/usr/bin/env Rscript
# Per-sample quality control: per-chromosome count distribution (aneuploidy,
# segmental duplication) and coverage evenness along each chromosome.
suppressPackageStartupMessages(library(Repliscope))
args <- commandArgs(trailingOnly = TRUE)
bedFile <- args[1]; sampleName <- args[2]

bed <- loadBed(bedFile, genome = sacCer3[["genome"]], name = sampleName)
med <- median(bed$score)
perChrom <- aggregate(score ~ chrom, data = bed, FUN = median)
perChrom$relative <- round(perChrom$score / med, 3)

writeLines(c(
  sprintf("sample\t%s", sampleName),
  sprintf("bins\t%d", nrow(bed)),
  sprintf("median_count\t%.1f", med),
  sprintf("mean_count\t%.1f", mean(bed$score)),
  sprintf("cv_percent\t%.2f", 100 * sd(bed$score) / mean(bed$score)),
  paste0("chrom_median_relative\t",
         paste(perChrom$chrom, perChrom$relative, sep = "=", collapse = ",")),
  sprintf("max_chrom_deviation\t%.3f", max(abs(perChrom$relative - 1)))
), paste0(sampleName, ".qc.txt"))

pdf(paste0(sampleName, ".qc.pdf"), width = 11, height = 5)
print(plotBed(bed, plotting = FALSE))
print(plotCoverage(bed, plotting = FALSE))
invisible(dev.off())
