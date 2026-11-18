#!/usr/bin/env bash
# Fetch UCSC sacCer3, assemble a single FASTA with Roman-numeral chromosome names
# (chrI-chrXVI plus chrM), write chrom.sizes and build the Bowtie 2 index.
set -euo pipefail

REF_DIR="${1:-reference}"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
ORDER="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI"

mkdir -p "$REF_DIR/chroms"
cd "$REF_DIR"

if [ ! -f chromFa.tar.gz ]; then
    echo "Downloading $URL"
    curl -fsSL -o chromFa.tar.gz "$URL"
fi
# Members sit at the top level of this archive; do not use --strip-components.
tar -xzf chromFa.tar.gz -C chroms

: > sacCer3.fa
for r in $ORDER; do cat "chroms/chr${r}.fa" >> sacCer3.fa; done
cat chroms/chrM.fa >> sacCer3.fa

awk '/^>/ {if (n) print n"\t"l; n=substr($1,2); l=0; next} {l+=length($0)} END {print n"\t"l}' \
    sacCer3.fa > sacCer3.chrom.sizes

echo "Assembled $(grep -c '^>' sacCer3.fa) sequences, $(awk '{s+=$2} END {print s}' sacCer3.chrom.sizes) bp"
md5sum sacCer3.fa 2>/dev/null || md5 sacCer3.fa

bowtie2-build --threads "${NCPU:-4}" -f sacCer3.fa sacCer3
samtools faidx sacCer3.fa
echo "Reference ready in $(pwd)"
