#!/usr/bin/env bash
set -euo pipefail

out_dir="${1:-raw}"
mkdir -p "$out_dir"
output="$out_dir/miRTarBase_human_reporter_2026-09-08.csv"

curl -L --fail --retry 5 --retry-all-errors --retry-delay 20 \
  -A 'Mozilla/5.0' \
  'https://awi.cuhk.edu.cn/miRTarBase/search/results/download/?mode=method&species=hsa&methods=reporter_assay' \
  -o "$output"

sha256sum "$output" > "$out_dir/SHA256SUMS"
