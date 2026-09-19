# Direct-interaction benchmark results

## Main result

The primary analysis used Huh-7.5 CLEAR-CLIP interactions annotated to 3'UTRs
and Huh7 APA profiles. Among direct-positive genes that also contained at least
one APAmiRank canonical 3'UTR site, the APA-exposure rank showed positive--
unlabeled AUCs of 0.648--0.671 for miR-16, miR-17, miR-25, and miR-186. All four
one-sided rank tests remained significant after Benjamini--Hochberg correction
across the five-miRNA family. miR-92a was weaker (AUC 0.551; q=0.159), although
its top-decile enrichment was 2.40-fold (hypergeometric P=0.0133).

| miRNA | Source 3'UTR genes | Observable genes | Coverage | PU-AUC | BH q | Top-10% enrichment | Top-10% P |
|---|---:|---:|---:|---:|---:|---:|---:|
| hsa-miR-16-5p | 165 | 40 | 24.2% | 0.671 | 2.76e-4 | 2.99 | 3.13e-4 |
| hsa-miR-17-5p | 360 | 109 | 30.3% | 0.652 | 2.73e-7 | 2.10 | 2.75e-4 |
| hsa-miR-25-3p | 36 | 10 | 27.8% | 0.667 | 0.0431 | 2.97 | 0.0708 |
| hsa-miR-92a-3p | 212 | 33 | 15.6% | 0.551 | 0.159 | 2.40 | 0.0133 |
| hsa-miR-186-5p | 96 | 29 | 30.2% | 0.648 | 0.00525 | 2.75 | 0.00603 |

Coverage is intentionally reported: APAmiRank ranks genes with a canonical seed
site in the modeled 3'UTR isoforms, whereas CLEAR-CLIP also captures noncanonical
interactions and genes outside that sequence universe.

## Independent CLASH check

After mapping GSE50452 transcript coordinates to hg19 `ensGene`, the HEK293
3'UTR subset showed weaker and heterogeneous effects. APA-exposure PU-AUC ranged
from 0.504 to 0.601; none of the five one-sided tests passed BH q<0.05. This is a
non-replication of the strong CLEAR-CLIP rank-distribution effect, not a technical
failure hidden from the report. Protocol differences, transcript-model mapping,
noncanonical interactions, and differences between the APA reference and the
CLASH experiment remain plausible contributors.

## APA-specific interpretation

The exposure-weighted additive rank was compared descriptively with an
unweighted additive rank. Median percentile gain was positive for all five
CLEAR-CLIP miRNAs, but the mean gain was negative for miR-17, miR-92a, and
miR-186. Therefore these data support direct-target prioritization most clearly;
they do **not** yet establish a uniform incremental benefit attributable solely
to APA weighting. A perturbation or matched condition with site-resolved 3'UTR
usage is still needed for that stronger claim.

## Legacy miR-124 / miR-155 audit

- GSE50452 contains two miR-124 chimeras, both mapped to CDS by the transcript
  coordinate audit; no 3'UTR direct benchmark is available.
- GSE73057 contains 11 gene-annotated miR-155 chimeras, but only one annotated
  3'UTR gene. That gene does not enter the canonical-site candidate universe.

The original miR-124/miR-155 perturbation benchmark and this five-miRNA direct
interaction benchmark therefore answer complementary questions and should be
reported as separate validation layers.
