#!/usr/bin/env Rscript
# Repliscope is published on CRAN but has no conda-forge/bioconda package,
# so it is installed into the environment's R library after env creation.
if (!requireNamespace("Repliscope", quietly = TRUE)) {
  install.packages("Repliscope", repos = "https://cloud.r-project.org")
}
cat("Repliscope", as.character(packageVersion("Repliscope")), "\n")
