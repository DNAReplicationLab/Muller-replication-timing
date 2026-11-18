#!/usr/bin/env Rscript
# Median replication time (Trep) from a sync-seq time course.
# Repliscope's calcTrep fits y ~ 2 + (1-2)/(1+exp((x-Trep)/slope)) per bin, with an
# anchor point (t=0, copy number 1). Its reference implementation scans the whole
# dataframe for every bin, which is O(n^2); the fit below is identical but indexes
# the bins first. --validate re-runs the package function on one chromosome and
# asserts the two agree.
suppressPackageStartupMessages(library(Repliscope))

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }

ratioFiles <- strsplit(getArg("--ratios"), ",")[[1]]
times      <- as.numeric(strsplit(getArg("--times"), ",")[[1]])
prefix     <- getArg("--prefix")
validChrom <- getArg("--validate", "chrI")
stopifnot(length(ratioFiles) == length(times))

message("Loading ", length(ratioFiles), " timepoints")
ratioDFs <- do.call(rbind, lapply(ratioFiles, function(f)
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
ratioDFs$chrom <- factor(ratioDFs$chrom, levels = unique(ratioDFs$chrom))

fitTrep <- function(df, times) {
  key <- paste(df$chrom, df$chromStart, df$chromEnd, sep = ":")
  bins <- unique(data.frame(chrom = df$chrom, chromStart = df$chromStart,
                            chromEnd = df$chromEnd, key = key, stringsAsFactors = FALSE))
  samples <- unique(df[, c("name.rep", "name.nonRep")])
  # wide matrix: one row per bin, one column per timepoint, in the supplied time order
  mat <- matrix(NA_real_, nrow = nrow(bins), ncol = nrow(samples))
  idx <- match(key, bins$key)
  scol <- match(paste(df$name.rep, df$name.nonRep),
                paste(samples$name.rep, samples$name.nonRep))
  mat[cbind(idx, scol)] <- df$ratio

  x <- c(0, times); a <- 2; b <- 1
  out <- vapply(seq_len(nrow(mat)), function(i) {
    y <- c(1, mat[i, ])
    fit <- NULL
    try(fit <- nls(y ~ a + (b - a) / (1 + exp((x - c) / d)),
                   start = list(c = 25, d = 1)), silent = TRUE)
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    p <- summary(fit)$parameters
    c(round(p[1, 1], 1), round(p[1, 2], 2))
  }, numeric(2))

  res <- data.frame(chrom = bins$chrom, chromStart = bins$chromStart,
                    chromEnd = bins$chromEnd, Trep = out[1, ], TrepErr = out[2, ])
  # Repliscope's sanity check: a Trep outside the sampled window is not supported
  bogus <- which(res$Trep < min(times) | res$Trep > max(times))
  res$Trep[bogus] <- NA; res$TrepErr[bogus] <- NA
  res
}

message("Fitting Trep for ", length(unique(paste(ratioDFs$chrom, ratioDFs$chromStart))), " bins")
TrepDF <- fitTrep(ratioDFs, times)

if (!is.na(validChrom) && validChrom != "none") {
  message("Validating against Repliscope::calcTrep on ", validChrom)
  ref  <- calcTrep(subset(ratioDFs, chrom == validChrom), times = times)
  mine <- subset(TrepDF, chrom == validChrom)
  ref  <- ref[order(ref$chromStart), ];  mine <- mine[order(mine$chromStart), ]
  stopifnot(nrow(ref) == nrow(mine), all(ref$chromStart == mine$chromStart))
  agree <- isTRUE(all.equal(ref$Trep, mine$Trep, tolerance = 1e-8))
  message(sprintf("  %d bins, identical: %s (NA counts %d vs %d)",
                  nrow(ref), agree, sum(is.na(ref$Trep)), sum(is.na(mine$Trep))))
  if (!agree) stop("fast Trep fit disagrees with Repliscope::calcTrep")
}

write.table(TrepDF, paste0(prefix, ".trep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ok <- TrepDF[!is.na(TrepDF$Trep), ]
write.table(data.frame(ok$chrom, ok$chromStart, ok$chromEnd, "Trep", ok$Trep),
            paste0(prefix, ".trep.bed"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
message(sprintf("Trep: %d of %d bins fitted (%.1f%%), range %.1f-%.1f min",
                nrow(ok), nrow(TrepDF), 100 * nrow(ok) / nrow(TrepDF),
                min(ok$Trep), max(ok$Trep)))

pdf(paste0(prefix, ".trep.pdf"), width = 12, height = 5)
print(plotTrep(TrepDF, plotting = FALSE))
invisible(dev.off())
