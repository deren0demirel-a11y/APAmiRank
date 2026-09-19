# Strict WT/mutant 3′UTR luciferase audit

This benchmark adds construct-level manual review to the broad miRTarBase
Reporter Assay benchmark. A pair enters the primary set only when accessible
primary full text verifies all of the following:

1. the tested sequence is a 3′UTR;
2. the named miRNA represses the wild-type reporter;
3. mutation or deletion of the cognate site removes or materially weakens that
   reporter effect; and
4. the comparison is specific to the named mature miRNA rather than inferred
   only from a seed family.

`curated/strict_luciferase_audit.tsv` retains included, excluded, and unresolved
records so selection decisions are inspectable. `construct_corroborated` records
are never mixed into the primary full-text scope; they appear only in a named
sensitivity scope.

## Reproduce

From the repository root:

```bash
python benchmarks/strict_luciferase/evaluate_strict_luciferase.py \
  --audit benchmarks/strict_luciferase/curated/strict_luciferase_audit.tsv \
  --ranks benchmarks/reporter_validation/evaluation/primary_reporter_positive_consensus_ranks.tsv \
  --output-dir benchmarks/strict_luciferase/evaluation
```

The rank table is frozen output from the v0.2.4 seven-miRNA reporter panel. The
evaluator performs a deterministic join and descriptive summaries; it does not
recompute APAmiRank scores.

## Interpretation boundary

This is a construct-validation pilot, not an unbiased performance test. Candidate
pairs were manually inspected after the broader reporter results were available,
creating selection bias. The very small number of eligible pairs also makes
inferential p-values misleading. Therefore the evaluator intentionally emits no
hypothesis test or confidence interval. A publication-grade next step must freeze
the literature-search protocol and eligibility decisions before ranks are opened,
then include every eligible pair found by that protocol.
