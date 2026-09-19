#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
raw_dir="${script_dir}/raw"
mkdir -p "${raw_dir}"

curl --fail --location --retry 3 \
  --output "${raw_dir}/hg19_refGene.txt.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
curl --fail --location --retry 3 \
  --output "${raw_dir}/hg19.2bit" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"
curl --fail --location --retry 3 \
  --output "${raw_dir}/twoBitToFa" \
  "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa"
chmod +x "${raw_dir}/twoBitToFa"

(
  cd "${raw_dir}"
  sha256sum hg19_refGene.txt.gz hg19.2bit twoBitToFa > HG19_REFERENCE_SHA256SUMS
)

printf 'Downloaded hg19 RefSeq coordinates, genome, and twoBitToFa to %s\n' "${raw_dir}"
