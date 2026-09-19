#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
raw_dir="${script_dir}/raw"
mkdir -p "${raw_dir}"

download() {
  local url="$1"
  local name="$2"
  curl --fail --location --retry 3 --output "${raw_dir}/${name}" "${url}"
}

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1268nnn/GSM1268942/suppl/GSM1268942_hela.15.sumCM.bed.gz" \
  "GSM1268942_hela.15.sumCM.bed.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1268nnn/GSM1268943/suppl/GSM1268943_hek293.15.sumCM.bed.gz" \
  "GSM1268943_hek293.15.sumCM.bed.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1268nnn/GSM1268944/suppl/GSM1268944_huh7.15.sumCM.bed.gz" \
  "GSM1268944_huh7.15.sumCM.bed.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1268nnn/GSM1268945/suppl/GSM1268945_imr90.15.sumCM.bed.gz" \
  "GSM1268945_imr90.15.sumCM.bed.gz"

for cell in HEK293 HUH7 HeLa IMR90; do
  download \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52530/suppl/GSE52530_${cell}.expData.txt.gz" \
    "GSE52530_${cell}.expData.txt.gz"
  download \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52530/suppl/GSE52530_${cell}.expData.qn.txt.gz" \
    "GSE52530_${cell}.expData.qn.txt.gz"
done

(
  cd "${raw_dir}"
  sha256sum ./*.gz > SHA256SUMS
)

printf 'Downloaded GSE52527/GSE52530 processed files to %s\n' "${raw_dir}"
