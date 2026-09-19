# Covariate-adjusted GSE52531 checkpoint — 2026-09-08

This is a prerelease analysis checkpoint, not a frozen manuscript result.

## Prespecified models

All comparisons use identical gene sets within each cell-line/miRNA analysis.
The outcome is `-log2(treatment/control)` without a pseudocount. The primary
expression threshold is mean control RPKM >=1.

The covariate block contains:

- log10 mean control RPKM;
- log1p maximum benchmark 3′UTR length;
- log1p canonical site count;
- mean raw 30-nt-flank AU fraction across sites.

Four ordinary least-squares models were compared: covariates only (C), C plus
unweighted architecture score (C+U), C plus exposure-weighted score (C+W), and
C plus both scores (C+U+W). The primary predictive metric is the median pooled
out-of-fold R² from 100 deterministic repetitions of 10-fold cross-validation.
The 2.5th--97.5th percentile range across repeats is a stability interval and
must not be described as a confidence interval because repeats reuse genes.

## Primary-threshold results

| Cell line | miRNA | Genes | C CV R² | C+U CV R² | C+W CV R² | C+U+W CV R² | C+W minus C+U | W beyond U |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| HeLa | miR-124 | 2,734 | 0.0634 | 0.1056 | 0.1421 | 0.1501 | +0.0365 | +0.0446 |
| HeLa | miR-155 | 1,573 | 0.1232 | 0.1528 | 0.2246 | 0.2250 | +0.0720 | +0.0725 |
| HEK293 | miR-124 | 2,968 | 0.0620 | 0.1191 | 0.1519 | 0.1635 | +0.0328 | +0.0444 |
| HEK293 | miR-155 | 1,773 | 0.1230 | 0.1618 | 0.2438 | 0.2431 | +0.0821 | +0.0813 |
| Huh7 | miR-124 | 1,613 | 0.0498 | 0.0948 | 0.1032 | 0.1164 | +0.0084 | +0.0214 |
| Huh7 | miR-155 | 776 | 0.0960 | 0.1169 | 0.1161 | 0.1224 | -0.0010 | +0.0055 |
| IMR90 | miR-124 | 1,922 | 0.0447 | 0.1243 | 0.1118 | 0.1363 | -0.0124 | +0.0121 |

`W beyond U` is the cross-validated R² difference between C+U+W and C+U. It
tests predictive information rather than mechanistic independence. The nested
in-sample partial F-test for adding W to C+U was nominally significant in all
six primary comparisons (largest P=0.0105 for Huh7--miR-155), but this is not a
substitute for an independent dataset and no multiplicity-adjusted claim is
made at this checkpoint.

The C+W versus C+U comparison favored exposure weighting in five of six primary
combinations. Huh7--miR-155 remained the exception at RPKM >=1 and was also
inconsistent across expression thresholds. IMR90--miR-124 is supporting rather
than part of the primary grid; C+W alone was inferior to C+U there, although W
added information when both scores were retained.

The full standardized design condition number was 4.07--5.40 at the primary
threshold, which does not indicate severe numerical instability. This does not
eliminate biological redundancy between U and W.

## TargetScan status

The comparison is designed for TargetScanHuman Release 8.0 context++ output.
During this checkpoint, the official bulk-download endpoint returned HTTP 502
and the documented CyTargetLinker/Figshare Release 8.0 mirror returned HTTP 403
from both command-line and browser environments. No unofficial or fabricated
scores were substituted. TargetScan therefore remains pending and must be run
when the versioned source file can be obtained and checksummed.

This matters conceptually because TargetScan 7/8 already incorporates 3′UTR
isoform profiles. It is a state-of-the-art APA-aware comparator, not a purely
APA-naive baseline. A separate longest-UTR total-context++ comparison is needed
to isolate the incremental value of APA weighting.

## Interpretation boundary

These results support the narrower statement that APA exposure contains
predictive information beyond simple gene-level expression, UTR length, site
count, local AU, and the unweighted APAmiRank score in this dataset. They do not
establish causal cell-type specificity, superiority to TargetScan, or
generalization beyond GSE52527/GSE52530.
