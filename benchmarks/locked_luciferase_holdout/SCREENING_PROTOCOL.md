# Locked luciferase holdout screening protocol

Protocol date: 2026-09-08  
Protocol identifier: `APAmiRank-locked-luciferase-holdout-v1`

## Status and claim boundary

This is a rank-independent, pair-level holdout within the already studied
seven-miRNA panel. It is not described as a fully prospective or investigator-
blinded validation because the broader APAmiRank reporter results existed before
this protocol. The candidate-selection script does not read any APAmiRank rank,
score, percentile, or outcome table.

## Source population

The source is the dated 2026-09-08 human miRTarBase Reporter Assay snapshot
already frozen in `benchmarks/reporter_validation`. Eligible records must:

1. match one of the seven prespecified mature human miRNAs;
2. be human-miRNA to human-target records;
3. have Reporter Assay evidence;
4. have Western blot or qPCR evidence;
5. be miRTarBase Tier A or Tier B;
6. report at least one paper in the snapshot; and
7. not occur anywhere in the v0.2.5 manual pilot audit, irrespective of that
   pilot record's inclusion status.

Within each miRNA, candidates are ordered by SHA-256 of
`seed + tab + miRNA + tab + gene`. The first five are selected. If fewer than
five candidates are eligible, every eligible candidate is retained. The seed
and full hash are stored in the manifest.

## Paper-level screening

Screeners must work only from `screening/screening_form.tsv`; APAmiRank rank
files must not be joined until every row has a final decision. For every pair,
retrieve the original article(s) underlying the listed MTI identifier and record
the DOI, PMID/PMCID, species, reporter vector, inserted region, cell context,
wild-type result, mutant design, mutant result, and orthogonal endogenous assay.

### Primary inclusion criteria

All criteria are required:

- the named mature miRNA, or an unambiguously identical mature sequence, is
  experimentally manipulated;
- the reporter insert derives from the target gene's 3′UTR;
- a wild-type construct is compared with a cognate site-mutant or site-deletion
  construct;
- the named miRNA represses the wild-type reporter relative to its control;
- mutation materially weakens or removes that reporter response; and
- the result can be verified in the original paper, supplement, or an official
  full construct record linked to that paper.

### Exclusion codes

`NOT_3UTR`, `NO_MUTANT`, `MUTANT_RESULT_UNCLEAR`, `FAMILY_ONLY`,
`WRONG_MIRNA`, `WRONG_TARGET`, `NON_ORIGINAL_SOURCE`, `RETRACTED`,
`FULL_TEXT_UNAVAILABLE`, and `OTHER_WITH_NOTE`.

Corrections, expressions of concern, and retractions must be checked and
recorded. A corrected article is not automatically excluded, but the correction
must be inspected for impact on the reporter evidence.

## Lock and analysis gate

The completed screening table must have no `pending` decisions and must receive
a SHA-256 checksum before evaluation. The evaluation script must fail closed if
any row remains pending, a required provenance field is missing for an included
record, or a selected pair is absent. Only then may the frozen APAmiRank rank
table be joined.

Primary reporting will include candidate accounting, observable-pair coverage,
APA percentile distributions, top-decile/top-quintile recall, and paired
APA-minus-unweighted percentile differences. Because multiple pairs share a
miRNA and literature selection remains non-random, pooled gene-pair p-values
will not be treated as independent-observation inference. Per-miRNA estimates
and cluster-aware uncertainty will be used only if the final eligible count is
adequate.
