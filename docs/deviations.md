# Deviations from the original pipeline

Every point where this workflow differs from `localMapper.sh` + the 2020 Nature Protocols
text, or from the 2013 analysis of the Müller *et al.* data, with the reason.

## Mapping

**`bowtie2 -a` dropped for single-end reads.** `localMapper.sh` passes `-a` (report every
alignment found) when mapping single-end data. Downstream, the counting step keeps only
reads with MAPQ ≥ 30 and no `XS:i:` tag, which excludes every multi-mapping read, so `-a`
cannot change the counts — it only makes the search exhaustive, costing orders of magnitude
in runtime and producing very large intermediate SAM files. Bowtie 2's default
(best-alignment) mode is used instead. Paired-end mapping is unchanged.

**`samtools markdup` instead of Picard `MarkDuplicates`.** Same operation — mark
coordinate-duplicate read pairs so the `-F 3840` filter drops them at the counting step —
but without a JVM dependency, and it runs natively on Apple Silicon. Requires the
`samtools fixmate -m` pass that is in the alignment module.

**5′-end intervals are written as proper 1 bp BED features.** `localMapper.sh` converts
`bedtools genomecov -5 -d` output with `print $1,$2,$2,...`, which makes a zero-length
interval at a 1-based coordinate. Here it is `print $1,$2-1,$2,...` — a valid 1 bp,
0-based half-open BED feature at the same position. This moves reads sitting exactly on a
bin boundary into the bin they belong to; it affects of order 0.1% of reads and leaves the
genome-wide profile unchanged.

## Colourspace data

Fourteen of the eighteen Müller *et al.* runs were sequenced on an AB SOLiD 5500 and are
stored in the SRA natively in colourspace (`cs_native="true"` in the run XML; the base
histogram is `0/1/2/3/.`, not ACGT). Bowtie 2 cannot read colourspace.

These runs are decoded to base space by `fastq-dump -B` and then mapped with Bowtie 2 like
everything else, so the whole dataset passes through one pipeline. The alternative —
mapping colourspace natively with Bowtie 1.2.3 `-C` — was rejected because that release is
unmaintained and x86-only, which would make the pipeline unreproducible on current hardware.

The cost is a lower alignment rate: a colour error propagates through the decoded read,
so some reads that a colourspace-aware aligner would place are lost. Because this analysis
counts reads in 1 kb bins and then takes a *ratio* between two samples prepared the same
way, a uniform loss of reads cancels. The measured rates are in `docs/muller2014.md`.

## Analysis

**Trep fitting is indexed, not scanned.** `Repliscope::calcTrep` searches the whole
dataframe for every bin and every timepoint, which is O(n²) and takes hours genome-wide.
`bin/repliscope_trep.R` reshapes to a bin × timepoint matrix first and fits the same
`nls` model with the same starting values, the same anchor point at (t=0, copy number 1)
and the same out-of-range sanity check. The script re-runs the package function on one
chromosome at every execution and fails if the two disagree.

**Trep model differs from the 2014 paper.** Müller *et al.* fitted
`F[t] = c + (d − c)/(1 + exp(b(t − e)))^f`, a five-parameter Richards curve.
Repliscope, and therefore this pipeline, fits the two-parameter Boltzmann sigmoid
`y = 2 + (1 − 2)/(1 + exp((t − Trep)/slope))` that the 2020 protocol describes. The
protocol is the specification being implemented here; the difference is noted because
Trep values will not be numerically identical to the 2014 figures.

**bigWig written with pyBigWig, not `bedGraphToBigWig`.** The bioconda osx-arm64 build of
the UCSC tool links against a Homebrew libpng and does not run. pyBigWig produces the same
format.

## Reference

UCSC sacCer3 (`chrI`–`chrXVI` plus `chrM`, 12,157,105 bp). The 2013 analysis used sacCer1
(`chr1`–`chr16`) — the published profiles in GEO are labelled `db=sacCer1`. Roman-numeral
naming was chosen to match Repliscope's bundled `sacCer3` object and the coordinates used
by OriDB.

Ratios are computed on the 16 nuclear chromosomes only. `Repliscope::sacCer3[["genome"]]`
masks out `chrM`, which is multi-copy and has no meaningful position on a 1–2 relative
copy number scale.

## Data acquisition

**Reads come from ENA, not sra-tools.** The protocol assumes FASTQ files already in hand.
Fetching them from the archive is new here, and it does not go through `prefetch` /
`fastq-dump`: sra-tools resolves a run to an S3 object, and an environment that cannot
reach that endpoint gets only `Failed to call external services` with no fallback. ENA
mirrors every SRA run, publishes an MD5 for each file and serves plain FASTQ over HTTPS,
so `bin/ena_fetch.py` resolves through the ENA portal API and transfers with curl.

Transfers of this size stall silently often enough to matter - a run that has moved no
bytes for ten minutes looks exactly like a slow one - so the fetcher uses curl's stall
detection (`--speed-limit`/`--speed-time`), resumes with `--continue-at`, verifies the
published MD5 before the file is used, and keeps raw downloads in a persistent cache
directory so an interrupted run resumes rather than restarts.

## Normalisation of the bundled dataset

**MFA-seq is normalised by a recovered factor, not by `upperLimit`.** The protocol sets
the MFA-seq ceiling to `1 + proportion of cells in S phase`, but that proportion is not
reported for the T7107 exponential culture. Setting `upperLimit` to 1.43, the maximum of
the published profile, makes Repliscope's minimiser return a factor of 1.216, whereas the
published profile actually carries 1.137 - the 2014 analysis predates Repliscope and used
"empirical normalisation to give a minimum value of 1". The recovered factor is used so
the MFA-seq profile sits on the same scale as the published data and as every other track
in the release. `norm_mode: upperLimit` remains available for datasets where the S-phase
proportion is known.

**Sync-seq factors are recovered, not measured.** Same reasoning: the flow-cytometry bulk
replication values the 2014 analysis used are not in the paper, but they are preserved in
the normalised profiles deposited in GEO. See `bin/recover_rfactor.R`.
