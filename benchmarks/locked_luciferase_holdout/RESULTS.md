# Locked luciferase holdout results

## Design and screening

The manifest contains 33 previously unaudited miRNA–target pairs selected
deterministically without APAmiRank rank or score columns. Original-paper
screening was completed under the frozen construct-level protocol before the
rank table was opened. Twenty-one pairs met every inclusion rule and 12 were
excluded. Exclusions are protocol or provenance failures and must not be read
as biological negatives.

The final screening-form SHA-256 is
`8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5`.
The evaluator rejects any different checksum, incomplete decision, altered
manifest membership, missing required inclusion field, or nonconforming
biological field.

## Prespecified evaluation

| Metric | Result |
|---|---:|
| Locked pairs | 33 |
| Strict inclusions | 21 |
| Observable strict inclusions | 12 |
| Candidate coverage | 0.5714 |
| Mean APA-exposure percentile | 0.7847 |
| Median APA-exposure percentile | 0.8185 |
| Top-decile fraction | 0.3333 |
| Top-quintile fraction | 0.5833 |
| Median APA minus unweighted percentile | 0.0651 |

The results are descriptive. The primary signal is that independently screened
direct targets tend to rank above the middle of their prespecified candidate
universes, while the comparison is limited by 57.1% coverage. Nine strict
inclusions lacked a row in the frozen rank table and remained missing. No score
was imputed, no candidate universe was changed, and no post-screening threshold
was introduced.

The positive median paired difference of 0.0651 describes the observed pairs;
it is not an inferential claim that APA exposure universally improves ranking.
The small, selected and incompletely observable holdout does not estimate a
calibrated target probability, sensitivity, specificity, or clinical utility.

## Exploratory secondary uncertainty analysis

After the descriptive holdout result was available, a secondary rank-uniform
null and miRNA-cluster bootstrap were added. This analysis is explicitly
post-unblinding and is not part of the prespecified primary evaluation. Its
complete status and assumptions are documented in
`SECONDARY_ANALYSIS_PROTOCOL.md`.

| Metric | Observed | Null mean | Raw P | Holm P | Cluster-equal estimate (95% bootstrap interval) |
|---|---:|---:|---:|---:|---:|
| Mean APA percentile | 0.7847 | 0.4992 | 0.00020 | 0.00080 | 0.7915 (0.7030–0.8847) |
| Median APA percentile | 0.8185 | 0.4993 | 0.00388 | 0.01164 | 0.7901 (0.7016–0.8847) |
| Top-decile fraction | 0.3333 | 0.0996 | 0.02534 | 0.02534 | 0.3750 (0.1250–0.6667) |
| Top-quintile fraction | 0.5833 | 0.1991 | 0.00438 | 0.01164 | 0.5833 (0.3333–0.8333) |
| Median APA minus unweighted | 0.0651 | — | — | — | 0.0423 (−0.0136–0.0899) |

Conditional on observability, the validated targets rank higher than expected
under the miRNA-matched rank-uniform null. The APA-minus-unweighted cluster
bootstrap interval crosses zero. Thus this holdout supports target-ranking
ability but does not independently establish an incremental advantage from APA
weighting. The P values remain exploratory because the inferential procedure
was designed after unblinding, the null is theoretical rather than an empirical
negative set, coverage is 57.1%, and only six miRNA clusters contribute.

## Reproducible outputs

| File | SHA-256 |
|---|---|
| `evaluation/locked_holdout_pair_ranks.tsv` | `7a7459ae95fe210488a6d63e80cdcee625ba6eed9b80ab6d083f8860d22477bf` |
| `evaluation/locked_holdout_summary.tsv` | `a16806c7712ea8192010acfec32376d760fc83c5b0de19bce1c7e0c90db46a1a` |

The secondary analysis additionally writes `holdout_uncertainty_summary.tsv`,
`holdout_mirna_cluster_summary.tsv`, compressed null/bootstrap draws, and the
publication figure in PNG, SVG and PDF formats. Their checksums are recorded in
the release verification document.

The pair-level table is the auditable source for all reported metrics. The
screening evidence and reasons for every inclusion or exclusion are retained in
`screening/screening_form.tsv` and summarized in `SCREENING_PROGRESS.md`.
