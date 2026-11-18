# Muller-replication-timing

Genome-wide replication timing profiles from DNA copy number, measured by high-throughput
sequencing: the **sort-seq**, **sync-seq** and **MFA-seq** family of methods.

This is a Nextflow DSL2 reimplementation of the pipeline published as

> Batrakou DG, Müller CA, Wilson RHC, Nieduszynski CA (2020).
> *DNA copy-number measurement of genome replication dynamics by high-throughput sequencing:
> the sort-seq, sync-seq and MFA-seq family.* Nature Protocols 15:1255–1284.
> [doi:10.1038/s41596-019-0287-7](https://doi.org/10.1038/s41596-019-0287-7)

which combined two tools:

| Original | Role | Here |
|---|---|---|
| [localMapper](https://github.com/DNAReplicationLab/localMapper) | map reads, count 5′ ends in fixed bins | `modules/bowtie2_align.nf`, `modules/count_bins.nf` |
| [Repliscope](https://github.com/DNAReplicationLab/Repliscope) | ratios, normalisation, smoothing, Trep | `bin/repliscope_*.R` (calls the CRAN package) |

The analysis logic is unchanged — Repliscope is called as a library, not reimplemented.
What is new is that acquisition, mapping, counting and analysis run as one resumable,
samplesheet-driven workflow, so a new dataset is added by adding rows to two CSV files.

## What it produces

For every sample: a BED file of uniquely-mapped read 5′ ends counted in 1 kb bins, plus
per-sample QC (per-chromosome count distribution, coverage evenness).

For every comparison: a sequencing-depth-normalised ratio of replicating to non-replicating
counts, put onto the relative copy number scale (1 = unreplicated, 2 = fully replicated),
spline-smoothed, written as TSV and BED and exported as bigWig.

For a sync-seq time course: a per-bin median replication time (Trep) track.

The bigWig tracks are the hand-off point to
[oridb-data-hub](https://github.com/DNAReplicationLab/oridb-data-hub), which owns their
display configuration. This repository stops at producing them.

## Quick start

```bash
conda env create -f env/environment.yml
conda activate Muller-replication-timing
Rscript bin/install_repliscope.R

# fetch and index the reference (UCSC sacCer3)
bash bin/fetch_reference.sh

nextflow run main.nf -profile conda \
    --samplesheet assets/samplesheet_muller2014.csv \
    --comparisons assets/comparisons_muller2014.csv \
    --outdir results
```

## Input format

`samplesheet.csv` — one row per sequencing run:

| column | meaning |
|---|---|
| `sample` | name used throughout the outputs |
| `run` | SRA run accession, or blank when using `--fastq_dir` |
| `platform` | `ILLUMINA` or `ABI_SOLID` (the latter triggers colourspace decoding) |
| `layout` | `single` or `paired` |
| `variant` | `sortseq`, `syncseq`, `mfaseq` or `hu` |
| `role` | `rep` (replicating) or `nonRep` (non-replicating control) |
| `strain`, `timepoint`, `series` | metadata, carried through to the outputs |

`comparisons.csv` — one row per profile to build:

| column | meaning |
|---|---|
| `comparison` | output name |
| `rep`, `nonRep` | sample names from the samplesheet |
| `norm_mode` | `auto` (sort-seq, HU), `rFactor` (sync-seq), `upperLimit` (MFA-seq) |
| `norm_value` | required for `rFactor` and `upperLimit`; see the protocol, step 61A(viii) |

## Normalisation

The ratio is first corrected for total read number, giving a genome-wide mean near 1.
It is then scaled to relative copy number according to the experimental variant:

- **sort-seq** and **HU origin mapping** — factor found algorithmically, by minimising the
  sum of values falling outside the 1–2 scale (`norm_mode: auto`).
- **sync-seq** — factor = 1 + (bulk genome replication as a decimal), measured by flow
  cytometry for that timepoint (`norm_mode: rFactor`, `norm_value: e.g. 1.35`).
- **MFA-seq** — upper limit = 1 + (proportion of cells in S phase)
  (`norm_mode: upperLimit`).

## Bundled dataset

`assets/samplesheet_muller2014.csv` and `assets/comparisons_muller2014.csv` describe the
18 sequencing runs of

> Müller CA, Hawkins M, Retkute R, Malla S, Wilson R, Blythe MJ, Nakato R, Komata M,
> Shirahige K, de Moura APS, Nieduszynski CA (2014).
> *The dynamics of genome replication using deep sequencing.* Nucleic Acids Research 42:e3.
> [doi:10.1093/nar/gkt878](https://doi.org/10.1093/nar/gkt878)

(GEO GSE42243 and GSE48212), remapped to sacCer3. See `docs/muller2014.md`.

## Documentation

- [`docs/method.md`](docs/method.md) — what each step does and why
- [`docs/deviations.md`](docs/deviations.md) — every point where this differs from the
  original pipeline or the 2013 analysis, and the reason
- [`docs/muller2014.md`](docs/muller2014.md) — the bundled dataset

## Licence

GPL-3, matching Repliscope.
