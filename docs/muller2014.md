# Muller et al. 2014 - the bundled dataset

> Muller CA, Hawkins M, Retkute R, Malla S, Wilson R, Blythe MJ, Nakato R, Komata M,
> Shirahige K, de Moura APS, Nieduszynski CA (2014). *The dynamics of genome replication
> using deep sequencing.* Nucleic Acids Research 42:e3. doi:10.1093/nar/gkt878

Eighteen sequencing runs across two GEO series, reprocessed onto sacCer3.

## Samples

| Sample | Run | Platform | Layout | Variant | Role | Strain |
|---|---|---|---|---|---|---|
| T7107_G2_SOLiD | SRR616663 | AB SOLiD 5500 | single | sort-seq | non-replicating | T7107 |
| T7107_S_SOLiD | SRR616664 | AB SOLiD 5500 | single | sort-seq | replicating | T7107 |
| T9475_G2_SOLiD | SRR616665 | AB SOLiD 5500 | single | sort-seq | non-replicating | T9475 |
| T9475_S_SOLiD | SRR616666 | AB SOLiD 5500 | single | sort-seq | replicating | T9475 |
| AUY077_alpha | SRR922015 | AB SOLiD 5500 | single | sync-seq | non-replicating | AUY077 |
| AUY077_25min to 90min | SRR922016-22 | AB SOLiD 5500 | single | sync-seq | replicating | AUY077 |
| CAY106_alpha | SRR922023 | AB SOLiD 5500 | single | HU | non-replicating | CAY106 |
| CAY106_HU | SRR922024 | AB SOLiD 5500 | single | HU | replicating | CAY106 |
| T7107_exp | SRR926343 | HiSeq 2000 | paired 2x51 | MFA-seq | replicating | T7107 |
| T7107_stat | SRR926344 | HiSeq 2000 | paired 2x51 | MFA-seq | non-replicating | T7107 |
| T9475_G2_Illumina | SRR926345 | HiSeq 2000 | paired 2x51 | sort-seq | non-replicating | T9475 |
| T9475_S_Illumina | SRR926346 | HiSeq 2000 | paired 2x51 | sort-seq | replicating | T9475 |

GSE42243 holds the four SOLiD sort-seq runs; GSE48212 holds the other fourteen.

T9475 appears twice, sequenced on SOLiD and again on Illumina. The paper uses that pair
to show the two platforms give the same replication profile, which makes it a useful
internal control on the reprocessing here.

## Colourspace

Fourteen of the eighteen runs are SOLiD and are archived in colourspace. They are decoded
to base space before mapping; see [deviations](deviations.md) for why, and
`tests/test_decode_colourspace.sh` for the check that the decoder matches NCBI's own
conversion byte for byte.

## Normalisation factors

The 2014 analysis normalised each sync-seq timepoint "to the bulk DNA content as measured
by flow cytometry". Those per-timepoint flow measurements are not reported in the paper,
but they are recoverable: the normalised profiles deposited in GEO carry them, so the
genome-wide mean of each published profile is the normalisation factor that was applied.
The values in `assets/comparisons_muller2014.csv` were obtained by regressing the
depth-normalised ratios produced here against the published profile for the same sample,
after lifting the published sacCer1 coordinates to sacCer3.

## Published profiles

GEO also holds the original normalised profiles as sacCer1 wig files at 1 kb resolution.
They are used here only for validation, correlating the new sacCer3 profiles against
them, and are not part of the released data.
