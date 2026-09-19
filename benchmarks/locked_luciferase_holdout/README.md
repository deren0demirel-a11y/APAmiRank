# Rank-independent locked luciferase holdout

This directory defines the next validation stage after the v0.2.5 manual pilot.
It creates a deterministic pair-level holdout without accessing APAmiRank ranks,
then separates literature screening from score evaluation.

Generate the locked candidate list from the repository root:

```bash
python benchmarks/locked_luciferase_holdout/lock_candidate_manifest.py \
  --evidence benchmarks/reporter_validation/derived/reporter_evidence.tsv \
  --prior-audit benchmarks/strict_luciferase/curated/strict_luciferase_audit.tsv \
  --output benchmarks/locked_luciferase_holdout/locked/candidate_manifest.tsv
```

The authoritative screening rules are in `SCREENING_PROTOCOL.md`. Do not join
the manifest to rank outputs until the screening form is complete and locked.

Generate rank-free PubMed candidates and an auditable query log:

```bash
python benchmarks/locked_luciferase_holdout/discover_pubmed_candidates.py \
  --manifest benchmarks/locked_luciferase_holdout/locked/candidate_manifest.tsv \
  --output benchmarks/locked_luciferase_holdout/screening/pubmed_candidate_articles.tsv \
  --query-log benchmarks/locked_luciferase_holdout/screening/pubmed_query_log.tsv
```

When PMC full-text access is available, extract passages that mention both a
reporter assay and construct-design language:

```bash
python benchmarks/locked_luciferase_holdout/extract_pmc_evidence.py \
  --candidates benchmarks/locked_luciferase_holdout/screening/pubmed_candidate_articles.tsv \
  --output benchmarks/locked_luciferase_holdout/screening/pmc_evidence_passages.tsv
```

The extractor is only a navigation aid and never assigns eligibility. Every
decision still requires manual verification under `SCREENING_PROTOCOL.md`.
Current progress and access limitations are recorded in `SCREENING_PROGRESS.md`.

The completed screening form was frozen before rank access with SHA-256
`8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5`.
Reproduce the fail-closed evaluation from the repository root with:

```bash
python benchmarks/locked_luciferase_holdout/evaluate_locked_holdout.py \
  --manifest benchmarks/locked_luciferase_holdout/locked/candidate_manifest.tsv \
  --screening benchmarks/locked_luciferase_holdout/screening/screening_form.tsv \
  --screening-sha256 8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5 \
  --ranks benchmarks/reporter_validation/evaluation/primary_reporter_positive_consensus_ranks.tsv \
  --output-dir benchmarks/locked_luciferase_holdout/evaluation
```

The result and its interpretation are documented in `RESULTS.md`. The generated
pair-level table retains unobservable candidates as missing; it does not impute
or manufacture ranks.

The post-unblinding secondary uncertainty analysis is separately labelled in
`SECONDARY_ANALYSIS_PROTOCOL.md`. Run it with:

```bash
python benchmarks/locked_luciferase_holdout/evaluate_holdout_uncertainty.py \
  --pair-ranks benchmarks/locked_luciferase_holdout/evaluation/locked_holdout_pair_ranks.tsv \
  --universe-summary benchmarks/reporter_validation/evaluation/reporter_consensus_summary.tsv \
  --output-dir benchmarks/locked_luciferase_holdout/evaluation \
  --iterations 50000 \
  --seed 20260910

python benchmarks/locked_luciferase_holdout/plot_locked_holdout.py \
  --pair-ranks benchmarks/locked_luciferase_holdout/evaluation/locked_holdout_pair_ranks.tsv \
  --screening benchmarks/locked_luciferase_holdout/screening/screening_form.tsv \
  --null-draws benchmarks/locked_luciferase_holdout/evaluation/holdout_rank_uniform_null.tsv.gz \
  --uncertainty-summary benchmarks/locked_luciferase_holdout/evaluation/holdout_uncertainty_summary.tsv \
  --output-prefix benchmarks/locked_luciferase_holdout/evaluation/locked_holdout_figure
```
