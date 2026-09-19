# Benchmark protocol (draft, freeze before outcome analysis)

## Claim under test

APA-aware isoform exposure improves the ranking of functional miRNA targets
beyond sequence architecture alone. The benchmark is not intended to claim
that an APA-exposed seed match is sufficient proof of direct targeting.

## Dataset roles

| Role | Dataset | Use |
|---|---|---|
| Primary functional benchmark | GSE52527 + GSE52530 | Cell-line 3P-seq context paired with miR-124/miR-155 perturbation RNA-seq |
| Direct-interaction benchmark | GSE73057 CLEAR-CLIP (primary) and GSE50452 CLASH (secondary) | Positive--unlabeled orthogonal physical-interaction support; see `benchmarks/direct_interactions/README.md` |
| Reporter benchmark | Dated official miRTarBase Reporter Assay export | Positive--unlabeled low-throughput validation; primary positives additionally require Western blot or qPCR; see `benchmarks/reporter_validation/README.md` |
| Strict reporter audit | Manually checked original articles | WT repression plus loss/weakening after cognate 3′UTR-site mutation; full-text and construct-corroborated scopes kept separate; see `benchmarks/strict_luciferase/README.md` |
| Locked reporter holdout | Deterministic 33-pair sample from the dated miRTarBase snapshot | Rank-independent pair selection, original-paper screening, and a fail-closed evaluation gate; see `benchmarks/locked_luciferase_holdout/SCREENING_PROTOCOL.md` |
| Curated strong evidence | miRTarBase, release to be frozen | Reporter-assay-supported miRNA–target pairs; evidence types retained separately |
| External baseline | TargetScan, release to be frozen | Established sequence-based ranking comparator |

Dataset releases, genome builds, mature miRNA sequences, download dates, URLs,
and checksums must be recorded before analysis.

## Frozen model variants

1. Longest annotated 3′UTR, sequence score only.
2. Dominant 3P-seq isoform, sequence score only.
3. APAmiRank architecture score without usage weighting.
4. APAmiRank exposure-weighted score.
5. TargetScanHuman Release 8.0 total context++ and cumulative weighted
   context++ scores as external comparators, reported separately. TargetScan
   7/8 already uses 3′UTR isoform profiles and must not be described as an
   APA-naive method.

No score weights may be optimized against the held-out perturbation, chimera,
or reporter labels. Post hoc annotations must remain outside primary ranking.

## Primary functional analysis

- Analyze each cell line and transfected miRNA separately.
- For Huh7, compare each miRNA with its own matched pUC19 control pair; do not
  pool the two control pairs.
- Preserve replicate-level normalized expression. Specify the low-expression
  filter and any pseudocount before calculating fold changes.
- Primary endpoint: improvement in held-out prediction of repression
  (`-log2FC`) from adding isoform exposure to a prespecified baseline model.
- The prespecified covariate block is log10 mean control expression, log1p
  maximum benchmark 3′UTR length, log1p canonical site count, and mean local
  AU fraction. Compare covariates only, covariates plus unweighted score,
  covariates plus exposure-weighted score, and covariates plus both scores on
  identical gene sets.
- Use 100 deterministic repetitions of 10-fold cross-validation for the
  GSE52530 checkpoint. Percentiles across repeated splits are stability
  intervals, not confidence intervals, because genes recur across repeats.
- Secondary endpoints: Spearman correlation, top-k repression, and enrichment
  of repressed genes across rank quantiles.
- Report effect sizes and bootstrap confidence intervals over genes. Use
  expression- and 3′UTR-length-matched permutations for enrichment tests.
- Use leave-one-cell-line-out evaluation where the model contains fitted
  parameters. Keep all transcripts from the same gene in one fold.

## Direct and curated evidence analyses

- Count a chimera-supported pair as positive only when the mature miRNA and
  mapped target locus agree after versioned identifier harmonization.
- Stratify curated pairs by assay. The strict primary subset requires a direct
  reporter assay. The construct-level pilot additionally requires an exact named
  miRNA, a 3′UTR reporter, WT repression, and loss or weakening after cognate
  site mutation. Full-text-verified and construct-corroborated records are
  reported separately.
- Treat the current strict luciferase audit as post-selection descriptive work:
  its literature candidates were inspected after the broad reporter ranks were
  available. A generalization analysis requires a prospectively frozen search
  and screening protocol before APAmiRank ranks are opened.
- Database or assay absence is **unlabeled**, not a confirmed negative.
  Therefore report recall/enrichment against a defined candidate universe and
  avoid interpreting ordinary AUROC against unlabeled pairs as specificity.
- Deduplicate evidence at the gene–miRNA level for gene-level analyses while
  retaining site-level records for localization analyses.

## Leakage and sensitivity controls

- Freeze candidate generation and ranking before joining outcome labels.
- Repeat analyses with common-only and extension-only sites.
- Repeat with/without thermodynamic features, local AU, and opportunity
  adjustment.
- Report identifier mapping loss, genes without measurable expression, and
  loci without usable APA estimates in a CONSORT-like accounting table.
- Treat the single 3P-seq library per cell line in GSE52527 as a limitation;
  it cannot estimate biological uncertainty in isoform usage.

## Minimum success criterion

The paper's central claim requires a reproducible improvement of the
exposure-weighted model over the unweighted architecture model in held-out
functional data, with directionally consistent support in at least one direct
interaction or reporter-evidence analysis. Failure to meet this criterion is
reported as a negative benchmark rather than repaired by outcome-guided score
tuning.
