# KUACC release smoke test

This test executes the pinned RNAhybrid 2.1.2 and ViennaRNA/RNAup 2.7.2 tools
through the Snakemake entry point on the synthetic example. It is a software
integration test and has no biological interpretation.

## Submit from the repository root

Upload or clone the complete APAmiRank repository under
`/scratch/users/hpc-ddemirel/`, enter its root directory, and submit:

```bash
sbatch scripts/hpc_release_smoke_test.slurm
```

If the uploaded copy does not contain its `.git` directory, pass the exact
GitHub commit SHA so it is retained in the provenance record:

```bash
APAMIRANK_GIT_COMMIT=PASTE_FULL_COMMIT_SHA \
  sbatch scripts/hpc_release_smoke_test.slurm
```

The script requests the `mid` partition for four hours, four CPU cores, and
16 GB RAM. It performs no compute on the login node. It creates a fresh Conda
environment and a time-stamped result directory below:

```text
/scratch/users/hpc-ddemirel/APAmiRank_release_validation/
```

It never deletes or overwrites an existing run or environment. To use an
already-created full environment, submit with an explicit prefix:

```bash
APAMIRANK_ENV_PREFIX=/scratch/users/hpc-ddemirel/path/to/environment \
  sbatch scripts/hpc_release_smoke_test.slurm
```

## Monitor and inspect

```bash
squeue -u "$USER"
tail -f apamirank_smoke_JOBID.out
```

After completion, locate the newest run and inspect its validation report:

```bash
ls -1dt /scratch/users/hpc-ddemirel/APAmiRank_release_validation/runs/* | head -1
cat /scratch/users/hpc-ddemirel/APAmiRank_release_validation/runs/RUN_ID/validation_report.json
```

A release-ready run must report `"status": "PASS"`, zero RNAhybrid/RNAup
parse failures, matching site/evidence counts, and successful verification of
all hashes recorded in `run_manifest.json`. Keep the complete time-stamped run
directory as release provenance.

If the job fails, do not weaken strict parsing. Preserve the `.out`, `.err`,
and time-stamped run directory so the executable output or parser incompatibility
can be diagnosed.
