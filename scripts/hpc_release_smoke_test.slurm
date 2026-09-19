#!/usr/bin/env bash
#SBATCH --job-name=apamirank_smoke
#SBATCH --partition=mid
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=apamirank_smoke_%j.out
#SBATCH --error=apamirank_smoke_%j.err

set -euo pipefail

repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
validation_root="${APAMIRANK_VALIDATION_ROOT:-/scratch/users/hpc-ddemirel/APAmiRank_release_validation}"

case "$validation_root" in
  /scratch/users/hpc-ddemirel/*) ;;
  *)
    echo "Refusing to write outside /scratch/users/hpc-ddemirel/: $validation_root" >&2
    exit 2
    ;;
esac

cd "$repo_root"
for required in pyproject.toml environment.yml workflow/Snakefile config/config.thermodynamic.example.yaml; do
  if [[ ! -f "$required" ]]; then
    echo "Submit this job from the APAmiRank repository root; missing: $required" >&2
    exit 2
  fi
done

run_id="$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID:-manual}"
run_dir="$validation_root/runs/$run_id"
results_dir="$run_dir/results"
env_prefix="${APAMIRANK_ENV_PREFIX:-$validation_root/envs/apamirank_v0.2.18_${SLURM_JOB_ID:-manual}}"

if [[ -e "$run_dir" ]]; then
  echo "Refusing to reuse an existing run directory: $run_dir" >&2
  exit 2
fi
mkdir -p "$run_dir" "$validation_root/envs"

if [[ -n "${APAMIRANK_ENV_PREFIX:-}" ]]; then
  if [[ ! -x "$env_prefix/bin/python" ]]; then
    echo "APAMIRANK_ENV_PREFIX does not contain an executable Python: $env_prefix" >&2
    exit 2
  fi
else
  if [[ -e "$env_prefix" ]]; then
    echo "Refusing to overwrite an existing environment: $env_prefix" >&2
    exit 2
  fi
  if command -v mamba >/dev/null 2>&1; then
    env_manager=mamba
  elif command -v conda >/dev/null 2>&1; then
    env_manager=conda
  else
    echo "Neither mamba nor conda is available. Load the cluster Conda module and resubmit." >&2
    exit 2
  fi
  "$env_manager" env create --prefix "$env_prefix" --file environment.yml
fi

# Snakemake rules execute in child shells. Export the selected environment so
# the APAmiRank CLI and pinned RNAhybrid/RNAup executables remain discoverable.
export PATH="$env_prefix/bin:$PATH"

cp config/config.thermodynamic.example.yaml "$run_dir/config.thermodynamic.yaml"
"$env_prefix/bin/python" - "$run_dir/config.thermodynamic.yaml" "$results_dir" "$repo_root/config" <<'PY'
from pathlib import Path
import sys
import yaml

config_path = Path(sys.argv[1])
results_dir = Path(sys.argv[2]).resolve()
template_base = Path(sys.argv[3]).resolve()
config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
for key, value in config.get("input", {}).items():
    if value and not Path(value).is_absolute():
        config["input"][key] = str((template_base / value).resolve())
config["output_dir"] = str(results_dir)
config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
PY

{
  echo "run_id=$run_id"
  echo "utc_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "hostname=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-not_available}"
  echo "slurm_partition=${SLURM_JOB_PARTITION:-not_available}"
  echo "slurm_cpus=${SLURM_CPUS_PER_TASK:-4}"
  echo "repository=$repo_root"
  echo "environment=$env_prefix"
  if command -v git >/dev/null 2>&1 && git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_commit=$(git -C "$repo_root" rev-parse HEAD)"
  elif [[ -n "${APAMIRANK_GIT_COMMIT:-}" ]]; then
    echo "git_commit=$APAMIRANK_GIT_COMMIT"
  else
    echo "git_commit=not_available"
  fi
} > "$run_dir/runtime_provenance.txt"

"$env_prefix/bin/python" --version > "$run_dir/python_version.txt" 2>&1
"$env_prefix/bin/pip" freeze > "$run_dir/pip_freeze.txt"
"$env_prefix/bin/apamirank" --version > "$run_dir/apamirank_version.txt"
{
  "$env_prefix/bin/RNAhybrid" --help 2>&1 || true
} > "$run_dir/rnahybrid_help.txt"
{
  "$env_prefix/bin/RNAup" --version 2>&1 || "$env_prefix/bin/RNAup" --help 2>&1 || true
} > "$run_dir/rnaup_version.txt"
{
  "$env_prefix/bin/snakemake" --version
} > "$run_dir/snakemake_version.txt"

"$env_prefix/bin/pytest" -q 2>&1 | tee "$run_dir/pytest.log"

"$env_prefix/bin/snakemake" \
  --snakefile workflow/Snakefile \
  --configfile "$run_dir/config.thermodynamic.yaml" \
  --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --dry-run 2>&1 | tee "$run_dir/snakemake_dry_run.log"

"$env_prefix/bin/snakemake" \
  --snakefile workflow/Snakefile \
  --configfile "$run_dir/config.thermodynamic.yaml" \
  --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --printshellcmds \
  --rerun-incomplete 2>&1 | tee "$run_dir/snakemake_run.log"

"$env_prefix/bin/python" scripts/validate_release_smoke_outputs.py \
  --results "$results_dir" \
  --report "$run_dir/validation_report.json" 2>&1 | tee "$run_dir/validation.log"

echo "utc_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$run_dir/runtime_provenance.txt"
echo "APAmiRank release smoke test PASS: $run_dir"
