#!/usr/bin/env bash
# The colourspace decoder must reproduce NCBI's own base-space conversion exactly.
#
# NCBI serves a decoded FASTQ for any run through the Run Browser backend; ENA serves
# the same run in its archived colourspace form. Decoding the ENA file with --pad-n
# (emit N from the first missing colour call, as fastq-dump -B does, rather than
# truncating) must give byte-identical sequence and quality strings.
#
# Verified for SRR616663 (T7107 G2, AB SOLiD 5500): 20000/20000 reads identical.
set -euo pipefail

RUN="${1:-SRR616663}"
N_READS="${2:-20000}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ENA always returns run_accession as the first column, so fastq_ftp is column 2.
ENA_URL=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=fastq_ftp&format=tsv" \
          | tail -1 | cut -f2 | cut -d';' -f1)
if [ -z "$ENA_URL" ]; then echo "FAIL: ENA returned no FASTQ path for $RUN"; exit 1; fi

curl -s "https://${ENA_URL}" | gzip -dc 2>/dev/null | head -n $((N_READS * 4)) > "$WORK/cs.fastq" || true
curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${RUN}" \
     | gzip -dc 2>/dev/null | head -n $((N_READS * 4)) > "$WORK/ncbi.fastq" || true

for f in cs ncbi; do
    if [ ! -s "$WORK/$f.fastq" ]; then echo "FAIL: no reads retrieved for $f"; exit 1; fi
done

python3 "$HERE/../bin/decode_colourspace.py" --in "$WORK/cs.fastq" --out "$WORK/mine.fastq" --pad-n --min-length 0

if cmp -s <(awk 'NR%4==2' "$WORK/mine.fastq") <(awk 'NR%4==2' "$WORK/ncbi.fastq") &&
   cmp -s <(awk 'NR%4==0' "$WORK/mine.fastq") <(awk 'NR%4==0' "$WORK/ncbi.fastq"); then
    echo "PASS: decoder matches NCBI base-space conversion for $RUN ($N_READS reads)"
else
    echo "FAIL: decoder output differs from NCBI base-space conversion for $RUN"
    exit 1
fi
