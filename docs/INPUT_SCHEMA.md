# Input schemas

## miRNA

The configuration requires an identifier and a mature miRNA sequence written
5′→3′. DNA `T` is converted to RNA `U`. Nucleotides 1–8 must be unambiguous.

## Isoform metadata TSV

Required columns:

| Column | Meaning |
|---|---|
| `isoform_id` | Unique identifier matching a FASTA record |
| `gene` | Gene symbol or stable gene identifier |
| `locus_id` | Alternative terminal-exon/3′UTR locus identifier |
| `pas_rank` | 1 for the shortest isoform, increasing toward distal PASs |
| `chrom` | Reference contig |
| `start0` | BED-like 0-based genomic start |
| `end` | BED-like exclusive genomic end |
| `strand` | `+` or `-` |

Ranks must be unique and contiguous within each `gene`/`locus_id`. Sequences
must be strand-correct transcript-oriented 5′→3′ sequences and their lengths
must equal `end - start0`.

## Isoform FASTA

One record per `isoform_id`. `apamirank prepare-qapa` can generate both the
metadata and FASTA from a QAPA BED7 file and an indexed genome FASTA.

## Condition-specific isoform usage TSV (optional)

Long-form table with one row per isoform and condition:

| Column | Meaning |
|---|---|
| `isoform_id` | Identifier present in the isoform metadata |
| `condition` | Sample, cell type, state, or treatment label |
| `usage_fraction` | Fraction from 0 to 1 assigned to this isoform within its locus and condition |

For every locus represented in the site ranking, each condition must contain a
row for every annotated isoform. Fractions must sum to 1 within the tolerance
set by `usage_sum_tolerance` (default 0.02). APAmiRank does not impute absent
rows or silently renormalize fractions.

The table can be constructed from QAPA PAU values or another validated 3′-end
quantification method, provided identifiers and reference builds are made
consistent before analysis. Biological replicates should be represented as
distinct condition labels in v0.2.2; replicate-aware uncertainty modelling is
planned for a later release.

## Eligible genes (optional)

Plain text with one gene identifier per line. Filtering occurs before site
enumeration. Define expression eligibility independently of target prediction.

## External site annotations (optional)

TSV with a unique `site_id` column. Other columns are added with an
`external_` prefix after the primary site ranking is fixed. This supports
external predictors such as TargetNet without mixing them into the primary
score.

## Post hoc gene annotations (optional)

TSV with a unique `gene` column. Differential-expression, rescue, pathway, or
experimental annotations may be supplied. Columns are added with a `posthoc_`
prefix after primary site and gene ranks are fixed.
