#!/usr/bin/env bash
set -euo pipefail

out_dir="${1:-raw}"
mkdir -p "$out_dir/GSE50452"

curl -L --fail --retry 3 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50452/suppl/GSE50452_RAW.tar \
  -o "$out_dir/GSE50452_RAW.tar"
tar -xf "$out_dir/GSE50452_RAW.tar" -C "$out_dir/GSE50452"

curl -L --fail --retry 3 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73057/suppl/GSE73057_Huh7_miRNA_chimera_interactions.xlsx \
  -o "$out_dir/GSE73057_Huh7_miRNA_chimera_interactions.xlsx"

curl -L --fail --retry 3 \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz \
  -o "$out_dir/hg19_ensGene.txt.gz"

curl -L --fail --retry 3 \
  https://www.mirbase.org/download/mature.fa \
  -o "$out_dir/mirbase_mature.fa"

sha256sum "$out_dir/GSE50452_RAW.tar" \
  "$out_dir/GSE73057_Huh7_miRNA_chimera_interactions.xlsx" \
  "$out_dir/hg19_ensGene.txt.gz" \
  "$out_dir/mirbase_mature.fa" > "$out_dir/SHA256SUMS"
