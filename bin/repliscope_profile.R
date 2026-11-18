#!/usr/bin/env Rscript
# Build a normalised, smoothed relative-copy-number profile from two binned count
# files, following Nature Protocols 2020 (doi:10.1038/s41596-019-0287-7) step 61A.
suppressPackageStartupMessages(library(Repliscope))

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

repBed     <- getArg("--rep")
nonRepBed  <- getArg("--nonRep")
repName    <- getArg("--repName")
nonRepName <- getArg("--nonRepName")
prefix     <- getArg("--prefix")
normMode   <- getArg("--norm", "auto")          # auto | rFactor | upperLimit
normValue  <- as.numeric(getArg("--value", NA))
loLim      <- as.numeric(getArg("--loLim", "0.25"))

stopifnot(!is.null(repBed), !is.null(nonRepBed), !is.null(prefix))

# ---- load and clean -------------------------------------------------------
# sacCer3[["genome"]] masks to the 16 nuclear chromosomes: the mitochondrial genome
# is multi-copy and not meaningful on a relative copy number scale.
loadClean <- function(path, name) {
  bed <- loadBed(path, genome = sacCer3[["genome"]], name = name)
  n0  <- nrow(bed)
  # Remove bins holding less than a quarter of the genome-wide median count.
  # These are repeats and low-mappability regions and they give noisy ratios.
  bed <- rmOutliers(bed, "median", loLim = loLim)
  message(sprintf("  %s: %d bins loaded, %d kept after outlier removal", name, n0, nrow(bed)))
  bed
}

message("Loading counts")
bedRep    <- loadClean(repBed,    repName)
bedNonRep <- loadClean(nonRepBed, nonRepName)

# ---- ratio ----------------------------------------------------------------
# makeRatio divides replicating by non-replicating counts bin-by-bin and corrects
# for the difference in total read number, giving a genome-wide mean close to 1.
ratioDF <- makeRatio(bedRep, bedNonRep)
message(sprintf("Ratio: %d shared bins, median %.4f", nrow(ratioDF), median(ratioDF$ratio)))

# ---- normalisation --------------------------------------------------------
# The scale is relative copy number: 1 = unreplicated, 2 = fully replicated.
ratioDF <- switch(normMode,
  # Sort-seq / HU: minimise the sum of values falling outside the 1-2 scale.
  "auto"       = normaliseRatio(ratioDF),
  # Sync-seq: rFactor = 1 + (bulk genome replication as a decimal), from flow cytometry.
  "rFactor"    = { stopifnot(!is.na(normValue)); normaliseRatio(ratioDF, rFactor = normValue) },
  # MFA-seq: maximum relative copy number = 1 + (proportion of cells in S phase).
  "upperLimit" = { stopifnot(!is.na(normValue)); normaliseRatio(ratioDF, upperLimit = normValue) },
  stop("unknown --norm mode: ", normMode)
)
rFactorUsed <- as.numeric(as.character(ratioDF$ratioFactor[1]))
message(sprintf("Normalisation: mode=%s factor=%.4f range=%.3f-%.3f",
                normMode, rFactorUsed, min(ratioDF$ratio), max(ratioDF$ratio)))

# ---- smoothing ------------------------------------------------------------
ratioDF <- smoothRatio(ratioDF)

# ---- outputs --------------------------------------------------------------
write.table(ratioDF, paste0(prefix, ".ratio.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# BED5 of the smoothed profile, for bigWig conversion and browser display
sm <- ratioDF[!is.na(ratioDF$splineSmooth), ]
write.table(data.frame(sm$chrom, sm$chromStart, sm$chromEnd,
                       paste0(repName, "_vs_", nonRepName), round(sm$splineSmooth, 5)),
            paste0(prefix, ".smooth.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
rw <- ratioDF[!is.na(ratioDF$ratio), ]
write.table(data.frame(rw$chrom, rw$chromStart, rw$chromEnd,
                       paste0(repName, "_vs_", nonRepName), round(rw$ratio, 5)),
            paste0(prefix, ".ratio.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

writeLines(c(
  sprintf("comparison\t%s_vs_%s", repName, nonRepName),
  sprintf("norm_mode\t%s", normMode),
  sprintf("ratio_factor\t%.4f", rFactorUsed),
  sprintf("bins\t%d", nrow(ratioDF)),
  sprintf("ratio_min\t%.4f", min(ratioDF$ratio)),
  sprintf("ratio_max\t%.4f", max(ratioDF$ratio)),
  sprintf("ratio_median\t%.4f", median(ratioDF$ratio))
), paste0(prefix, ".normstats.txt"))

pdf(paste0(prefix, ".profile.pdf"), width = 12, height = 7)
print(plotGenome(ratioDF, plotting = FALSE, region = "chrX",
                 circles = sacCer3[["ori"]], lines = sacCer3[["cen"]]))
print(plotGenome(ratioDF, plotting = FALSE))
invisible(dev.off())
message("Done: ", prefix)
