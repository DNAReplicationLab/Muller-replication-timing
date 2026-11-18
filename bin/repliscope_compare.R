#!/usr/bin/env Rscript
# Compare two replication profiles bin by bin with Repliscope's z-score test.
#
# compareRatios() asks, for each bin, whether the two profiles differ by more than their
# sampling noise allows, and returns a p value per bin. Use it to test a mutant against
# its wild type, or - as here - to test whether two independent measurements of the same
# sample agree.
suppressPackageStartupMessages(library(Repliscope))

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }

f1 <- getArg("--ratio1"); f2 <- getArg("--ratio2")
l1 <- getArg("--label1", "profile1"); l2 <- getArg("--label2", "profile2")
prefix <- getArg("--prefix", "comparison")
alpha  <- as.numeric(getArg("--alpha", "0.01"))

read1 <- function(p) {
  d <- read.table(p, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d$chrom <- factor(d$chrom, levels = unique(d$chrom))
  d
}
r1 <- read1(f1); r2 <- read1(f2)

cmp <- compareRatios(r1, r2)
write.table(cmp, paste0(prefix, ".compared.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- cmp$p.value[!is.na(cmp$p.value)]
# one p value per bin, so correct for testing every bin in the genome
adj <- p.adjust(p, method = "BH")
out <- c(
  sprintf("profile1\t%s", l1),
  sprintf("profile2\t%s", l2),
  sprintf("bins_tested\t%d", length(p)),
  sprintf("alpha\t%.3f", alpha),
  sprintf("bins_differing_raw\t%d", sum(p < alpha)),
  sprintf("bins_differing_BH\t%d", sum(adj < alpha)),
  sprintf("percent_differing_BH\t%.2f", 100 * sum(adj < alpha) / length(p))
)
writeLines(out, con = stdout())
writeLines(out, paste0(prefix, ".compared.stats.txt"))

pdf(paste0(prefix, ".compared.pdf"), width = 12, height = 6)
print(plotGenome(cmp, plotting = FALSE, region = "chrX", lines = sacCer3[["cen"]]))
invisible(dev.off())
