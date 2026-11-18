#!/usr/bin/env Rscript
# Recover the normalisation factor that was applied to a previously published profile.
#
# Repliscope's normaliseRatio multiplies the depth-normalised ratio by a scalar to put
# it on the relative copy number scale. For sort-seq that scalar is found
# algorithmically, but for sync-seq it comes from a flow-cytometry measurement of bulk
# replication, and for MFA-seq from the proportion of cells in S phase. When those
# measurements were not published, the factor can still be recovered from a published
# normalised profile of the same sample: build the unnormalised ratio here, regress the
# published values on it over shared bins, and the slope is the factor that was used.
#
# The regression is reported with its correlation so a poor match is visible rather
# than silently producing a confident-looking number.
suppressPackageStartupMessages(library(Repliscope))

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }

repBed    <- getArg("--rep")
nonRepBed <- getArg("--nonRep")
reference <- getArg("--reference")      # published profile, BED, same assembly and bins
label     <- getArg("--label", "sample")
loLim     <- as.numeric(getArg("--loLim", "0.25"))
outFile   <- getArg("--out")

bedRep    <- rmOutliers(loadBed(repBed,    genome = sacCer3[["genome"]], name = "rep"),
                        "median", loLim = loLim)
bedNonRep <- rmOutliers(loadBed(nonRepBed, genome = sacCer3[["genome"]], name = "nonRep"),
                        "median", loLim = loLim)
ratioDF <- makeRatio(bedRep, bedNonRep)      # depth-corrected, genome mean near 1

ref <- read.table(reference, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(ref)[1:3] <- c("chrom", "chromStart", "chromEnd")
ref$refValue <- as.numeric(ref[[ncol(ref)]])
ref <- ref[, c("chrom", "chromStart", "chromEnd", "refValue")]

m <- merge(ratioDF, ref, by = c("chrom", "chromStart", "chromEnd"))
m <- m[is.finite(m$ratio) & is.finite(m$refValue), ]
if (nrow(m) < 1000) stop("only ", nrow(m), " shared bins - check assembly and bin size")

# Through-the-origin fit: both profiles are copy number, so a scale factor and no offset
slope    <- sum(m$ratio * m$refValue) / sum(m$ratio^2)
fit      <- lm(refValue ~ ratio, data = m)       # with intercept, as a sanity check
r        <- cor(m$ratio, m$refValue)
meanRef  <- mean(m$refValue)
meanOurs <- mean(m$ratio)

out <- c(
  sprintf("label\t%s", label),
  sprintf("shared_bins\t%d", nrow(m)),
  sprintf("pearson_r\t%.5f", r),
  sprintf("mean_unnormalised_ratio\t%.5f", meanOurs),
  sprintf("mean_reference\t%.5f", meanRef),
  sprintf("rFactor_through_origin\t%.4f", slope),
  sprintf("rFactor_mean_ratio\t%.4f", meanRef / meanOurs),
  sprintf("lm_slope\t%.4f", coef(fit)[2]),
  sprintf("lm_intercept\t%.4f", coef(fit)[1])
)
writeLines(out, con = stdout())
if (!is.null(outFile)) writeLines(out, outFile)
