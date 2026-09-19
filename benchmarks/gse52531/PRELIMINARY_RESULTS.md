# Preliminary GSE52531 benchmark — 2026-09-07

These results are a pipeline-development checkpoint, not a frozen final
analysis. They use APAmiRank architecture mode (canonical seed class and local
AU), UCSC hg19 RefSeq coordinates, and GSE52530 quantile-normalized RPKM.

## Data preparation

GSE52527 clusters were interpreted as mapped 3P-tag intervals. The PAS boundary
was BED `end` for plus-strand records and BED `start` for minus-strand records.
Clusters with fewer than two tags and assignments overlapping more than one
RefSeq terminal-UTR locus were excluded. A locus was retained as APA when it had
at least 10 mapped tags and at least two PASs each reaching 1% usage. Minor PASs
within qualifying loci remained in the denominator.

| Cell line | Conservative APA loci | Retained PAS isoforms |
|---|---:|---:|
| HeLa | 5,105 | 19,480 |
| HEK293 | 5,348 | 21,013 |
| Huh7 | 3,354 | 11,014 |
| IMR90 | 3,655 | 11,690 |

The common atlas required at least 10 tags in every cell line. PAS boundaries
within a maximum 30-nt span were collapsed, yielding 3,053 loci and 13,947 PAS
isoforms.

## Per-cell-line preliminary comparison

For each gene, the most highly expressed RefSeq record in the matched control
was selected independently of miRNA outcome. Genes required mean control RPKM
≥1 and positive mean expression in treatment and control. `Repression` was
defined as `-log2(treatment/control)` without a pseudocount.

| Cell line | miRNA | Genes | Unweighted R² | Exposure-weighted R² | ΔR² (95% paired bootstrap CI) |
|---|---|---:|---:|---:|---:|
| HeLa | miR-124 | 2,734 | 0.0371 | 0.0798 | +0.0428 (0.0251, 0.0607) |
| HeLa | miR-155 | 1,573 | 0.0806 | 0.1612 | +0.0807 (0.0564, 0.1062) |
| HEK293 | miR-124 | 2,968 | 0.0638 | 0.1054 | +0.0415 (0.0289, 0.0552) |
| HEK293 | miR-155 | 1,773 | 0.1266 | 0.2137 | +0.0871 (0.0632, 0.1115) |
| Huh7 | miR-124 | 1,613 | 0.0449 | 0.0703 | +0.0254 (0.0082, 0.0439) |
| Huh7 | miR-155 | 776 | 0.0742 | 0.0715 | −0.0027 (−0.0303, 0.0249) |
| IMR90 | miR-124 | 1,922 | 0.0822 | 0.0976 | +0.0154 (0.0005, 0.0329) |

IMR90–miR-124 is supporting rather than part of the six-combination primary
grid. IMR90–miR-155 is excluded from the confirmatory summary because the main
paper describes the human primary design as both miRNAs in HeLa, HEK293 and
Huh7, with IMR90–miR-124 as an additional analysis.

## Common-atlas cognate-profile control

On the common PAS atlas, every one of the six primary combinations had higher
Pearson R² with the cognate exposure-weighted score than with the unweighted
score. However, the cognate cell-line profile was the best of the four exposure
profiles in only 2/6 combinations (both HeLa outcomes). Therefore this
checkpoint supports a general benefit of isoform-abundance weighting but does
**not** yet support a strong claim that the cognate cell-line AIR profile is
consistently superior to non-cognate profiles.

## Limitations to resolve before manuscript-level inference

- The original study's exact historical RefSeq snapshot and full wContext+
  implementation have not yet been reproduced.
- The architecture-only score is deliberately simpler than TargetScan
  context+/context++ and is not a like-for-like replication.
- The 30-nt cross-sample PAS harmonization window requires sensitivity analysis.
- Transcript-to-gene outcome collapse currently selects the highest-control
  RefSeq record; alternative collapse rules must be tested.
- Covariate-adjusted and held-out models controlling for site number, 3′UTR
  length, baseline expression, and AU composition remain to be run.
- Huh7–miR-155 does not show a robust improvement in the per-cell-line analysis.

The correct current interpretation is that this benchmark is promising and
technically coherent, but not yet sufficient to establish the final biological
performance claim.
