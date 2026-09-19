# Reporter-assay validation results

## High-confidence primary set

The primary set required human--human miRTarBase Reporter Assay evidence plus
Western blot or qPCR support. In the at-least-two-context consensus universe,
APA-exposure PU-AUC ranged from 0.722 to 0.899 among miRNAs with at least five
observable positive genes. All six evaluable miRNAs passed BH q<0.001; miR-186
had only two observable positives and is not interpreted inferentially.

| miRNA | Source positives | Observable | Coverage | APA PU-AUC | BH q | Top-10% enrichment |
|---|---:|---:|---:|---:|---:|---:|
| hsa-miR-16-5p | 49 | 25 | 51.0% | 0.830 | 2.35e-8 | 3.99 |
| hsa-miR-17-5p | 45 | 24 | 53.3% | 0.782 | 2.23e-6 | 4.17 |
| hsa-miR-25-3p | 25 | 9 | 36.0% | 0.824 | 4.48e-4 | 6.64 |
| hsa-miR-92a-3p | 29 | 8 | 27.6% | 0.899 | 8.76e-6 | 6.23 |
| hsa-miR-124-3p | 74 | 34 | 45.9% | 0.802 | 4.82e-9 | 2.94 |
| hsa-miR-155-5p | 72 | 23 | 31.9% | 0.722 | 1.74e-4 | 2.17 |
| hsa-miR-186-5p | 10 | 2 | 20.0% | 0.898 | not tested | 5.00 |

## Increment attributable to APA exposure

APA-exposure PU-AUC exceeded unweighted-additive PU-AUC for every panel miRNA,
but paired bootstrap evidence was not uniform. The gain remained positive after
BH correction for miR-16, miR-25, and miR-124. Confidence intervals crossed
zero for miR-17, miR-92a, and miR-155; miR-186 was suppressed because n=2.

| miRNA | APA minus unweighted AUC | Bootstrap 95% CI | BH q |
|---|---:|---:|---:|
| hsa-miR-16-5p | +0.048 | 0.005 to 0.087 | 0.0340 |
| hsa-miR-17-5p | +0.021 | -0.063 to 0.092 | 0.276 |
| hsa-miR-25-3p | +0.070 | 0.017 to 0.128 | 0.0150 |
| hsa-miR-92a-3p | +0.040 | -0.028 to 0.104 | 0.182 |
| hsa-miR-124-3p | +0.060 | 0.026 to 0.093 | 0.00120 |
| hsa-miR-155-5p | +0.036 | -0.037 to 0.100 | 0.185 |
| hsa-miR-186-5p | not tested | fewer than five positives | not tested |

The defensible conclusion is that APAmiRank strongly prioritizes previously
reporter-supported targets within its canonical-site universe. These data also
suggest an incremental benefit from APA exposure for a subset of miRNAs, but
do not support a universal APA-improvement claim.

## Sensitivities

Using all reporter records rather than the primary orthogonally supported set
gave APA PU-AUCs of 0.744--0.838 for miR-16, miR-17, miR-25, miR-92a, miR-124,
and miR-155. miR-186 again had only three observable positives. The strict
all-four-context universe preserved strong rank effects where enough positives
remained, but eliminated all observable primary positives for miR-25 and
miR-186 and left only one for miR-92a; it is therefore informative as a
coverage sensitivity rather than the main benchmark.
