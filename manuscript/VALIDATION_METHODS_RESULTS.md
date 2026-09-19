# APAmiRank: validation methods, results, figure legend, and claim boundaries

**Manuscript-ready working draft — 10 September 2026**  
**Software version:** APAmiRank v0.2.18  
**Scope:** validation analyses only; title, abstract, introduction, software implementation, discussion, and author metadata remain to be integrated into the complete manuscript.

## Purpose and evidential scope

APAmiRank is an alternative-polyadenylation-aware framework for ranking candidate microRNA (miRNA) targets. The validation programme was designed to distinguish two questions that are related but not interchangeable: whether APAmiRank prioritizes experimentally supported targets, and whether APA exposure weighting adds predictive information beyond the corresponding unweighted site-architecture score. Five complementary evidence layers were used: miRNA-perturbation RNA sequencing, ligation-based direct-interaction assays, curated reporter evidence, a manually verified construct-level pilot, and a rank-independent locked luciferase holdout. Because most benchmarks are positive–unlabeled, the analyses evaluate prioritization among eligible canonical-site genes and do not estimate specificity, a false-positive rate, or a calibrated probability of targeting.

The combined evidence supports APAmiRank as a target-prioritization tool. Evidence for an incremental contribution from APA weighting is positive in the perturbation benchmark and for a subset of reporter-supported miRNAs, but is not uniform across miRNAs or independently established by the locked holdout. The current analyses do not support claims of superiority to TargetScan, universal improvement from APA weighting, causal cell-type specificity, or clinical utility.

## Methods

### APAmiRank candidate universe and ranking

Strand-correct 3′UTR isoforms were ordered from proximal to distal polyadenylation-site usage within each locus. Canonical 8mer, 7mer-m8, 7mer-A1, and 6mer motifs were generated from the supplied mature miRNA sequence and enumerated across eligible isoforms. Genomically identical sites present in multiple isoforms were collapsed while retaining their first and last available polyadenylation-site ranks and their classification as common or extension-specific sites. For condition-specific analyses, the exposure of a site was the summed usage fraction of all isoforms carrying that site. The fixed site-evidence score was multiplied by exposure, and exposure-weighted evidence was summed across sites to obtain a gene-level score. The unweighted comparator used the same eligible sites and evidence scores without the exposure multiplier. Scores were interpreted as relative prioritization statistics, not calibrated probabilities.

Percentile ranks were calculated within each miRNA-specific canonical-site candidate universe. Consequently, genes without an enumerated canonical 3′UTR site were outside the modeled universe rather than classified as negative. Unless stated otherwise, higher percentiles indicate higher APAmiRank priority.

### Functional perturbation benchmark (GSE52527/GSE52530)

The functional benchmark combined processed human 3P-seq poly(A)-site profiles from GSE52527 with quantile-normalized RNA-sequencing data after miR-124, miR-155, or mock transfection from GSE52530, the two relevant subseries of GSE52531. HeLa, HEK293, Huh7, and IMR90 contexts were processed on hg19. There was one processed 3P-seq profile per cell type; miRNA-transfection RNA sequencing comprised two biological replicates per treatment, and the separate Huh7 mock pairs corresponding to miR-124 and miR-155 were retained rather than pooled. The primary expression threshold was mean control RPKM ≥1. The response was −log2(treatment/control), calculated without a pseudocount.

Four ordinary least-squares models were compared on identical gene sets: covariates only (C), C plus the unweighted architecture score (C+U), C plus the exposure-weighted score (C+W), and C plus both scores (C+U+W). Covariates were log10 mean control RPKM, log1p maximum benchmark 3′UTR length, log1p canonical site count, and mean raw AU fraction in 30-nt site flanks. Predictive performance was the median pooled out-of-fold R² over 100 deterministic repetitions of 10-fold cross-validation. The 2.5th–97.5th percentile range across repetitions was treated as a stability interval, not a confidence interval, because genes were reused across repetitions. The difference between C+U+W and C+U quantified predictive information contributed by W conditional on U; it was not interpreted as mechanistic independence. In-sample partial F tests were secondary and were not used as substitutes for external validation.

### Direct-interaction benchmarks (GSE73057 and GSE50452)

Ligation-based interaction data were evaluated as positive–unlabeled benchmarks. The primary analysis used CLEAR-CLIP interactions from Huh-7.5 cells (GSE73057) restricted to annotated 3′UTRs and compared them with Huh7 APA profiles; Huh-7.5 and Huh7 were regarded as near-matched, not identical. The independent check used CLASH interactions from HEK293 cells (GSE50452). CLASH transcript coordinates were mapped to hg19 UCSC ensGene models and assigned to 5′UTR, coding sequence, or 3′UTR; the 3′UTR subset was the relevant comparison. The prespecified panel comprised miR-16-5p, miR-17-5p, miR-25-3p, miR-92a-3p, and miR-186-5p.

For each miRNA, direct-interaction genes were compared with all other genes in the same canonical-site universe. End points were positive–unlabeled area under the receiver-operating-characteristic curve (PU-AUC), a one-sided Mann–Whitney rank test, and top-decile and top-quintile hypergeometric enrichment. Mann–Whitney P values were adjusted by the Benjamini–Hochberg procedure across the five-miRNA family within each dataset, region, and rank-metric family. Absence from CLEAR-CLIP or CLASH was not treated as a confirmed negative.

### Curated reporter-assay benchmark

A human Reporter Assay export was obtained from the official miRTarBase portal on 8 September 2026. Because the download endpoint was not version tagged, the file was treated as a dated portal snapshot rather than asserted to be an immutable miRTarBase v10 file. Of 7,790 rows, 7,669 were exact human-miRNA–human-target records; the seven-miRNA panel contributed 513 unique harmonized pairs. The panel consisted of miR-16-5p, miR-17-5p, miR-25-3p, miR-92a-3p, miR-124-3p, miR-155-5p, and miR-186-5p.

The primary positive set required Reporter Assay evidence plus Western blot or qPCR support. Because the literature spanned heterogeneous cell types and the export did not provide a usable cell line for every record, interactions were not designated as cell matched. The primary consensus universe required each gene to be rankable in at least two of four APA contexts (HeLa, HEK293, Huh7, and IMR90); available percentile ranks were averaged and the consensus score was re-ranked. Requiring all four contexts and retaining all reporter records were sensitivity analyses.

Primary end points were candidate coverage, PU-AUC, one-sided Mann–Whitney tests with Benjamini–Hochberg adjustment across the seven-miRNA family, and top-decile/top-quintile enrichment. The incremental APA analysis compared APA-exposure and unweighted-additive PU-AUCs on identical genes using 5,000 paired stratified bootstrap samples of positive and unlabeled genes. Inferential results were suppressed when fewer than five positive genes were observable.

### Construct-level pilot

Candidate pairs were manually audited for a human target 3′UTR reporter, repression of a wild-type construct, and loss or material weakening of the response after mutation of the cognate miRNA site. Four pairs met all criteria: miR-16-5p–BCL2, miR-17-5p–BMPR2, miR-25-3p–PTEN, and miR-124-3p–RAB27A. A separately labeled sensitivity scope added miR-155-5p–ETS1 when construct-level corroboration was accepted. The audit was performed after the broad reporter benchmark had been examined and was therefore reported descriptively without a P value or confidence interval.

### Rank-independent locked luciferase holdout

To reduce outcome-aware selection, a deterministic manifest of 33 previously unaudited miRNA–target pairs was generated without rank or score columns. Original papers were screened manually under a frozen construct-level protocol before the rank table was opened. Eligibility required an exact mature miRNA, a human target, a reporter containing the cognate human target 3′UTR sequence, wild-type repression, and loss or material weakening after mutation of the cognate site. Fail-closed exclusions were retained as protocol or provenance failures and were not reclassified as biological negatives. The completed screening form was frozen with SHA-256 `8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5` before joining the eligible pairs to the prespecified consensus rank table.

The prespecified evaluation was descriptive: coverage, mean and median APA-exposure percentile, fractions in the top decile and top quintile, and the paired median difference between APA-exposure and unweighted percentiles. Missing pairs were retained as missing; no scores were imputed and the candidate universe was not altered.

After the descriptive results had been opened, a separately labeled exploratory uncertainty analysis was designed. A theoretical rank-uniform null preserved, for each of the six represented miRNAs, the number of observable holdout pairs and the size of that miRNA’s full consensus ranking universe. In each of 50,000 iterations, ordinal rank positions were sampled without replacement and converted to percentiles. One-sided Monte Carlo P values used a plus-one correction, and the four correlated tests were adjusted by Holm’s method. This null did not preserve empirical score-tie multiplicities and was conditional on observability.

For cluster-aware sensitivity analysis, miRNA-specific statistics were calculated first and the six miRNA clusters were sampled with replacement for 50,000 iterations with equal miRNA weight. Percentile 2.5% and 97.5% quantiles were reported. For APA-minus-unweighted ranking, the cluster statistic was the within-miRNA median paired difference. With six clusters, these intervals were treated as exploratory sensitivity estimates.

### Reproducibility and statistical interpretation

All input transformations, eligibility rules, model specifications, random seeds, multiple-testing families, rank tables, and pair-level joins are stored with the benchmark code. Candidate coverage was reported throughout. Tests used the prespecified candidate universe and did not manufacture negative labels. P values are nominal unless an adjustment method is stated explicitly. The locked holdout’s secondary P values are post-unblinding and exploratory even after Holm adjustment. All tests in the v0.2.18 software package passed at release.

## Results

### APA exposure contributed predictive information in the functional perturbation benchmark

At mean control RPKM ≥1, six primary cell-line/miRNA combinations contained 776–2,968 genes. Adding the exposure-weighted score to the unweighted model increased cross-validated R² in all six primary comparisons: the C+U+W minus C+U difference ranged from +0.0055 for Huh7–miR-155 to +0.0813 for HEK293–miR-155. Direct comparison of C+W with C+U favored exposure weighting in five of six primary combinations; Huh7–miR-155 was the exception (difference −0.0010). The largest improvements were observed for miR-155 in HeLa (+0.0720 for C+W versus C+U) and HEK293 (+0.0821). In the supporting IMR90–miR-124 comparison, C+W alone was inferior to C+U (−0.0124), whereas W retained a positive increment when added to U (+0.0121). Thus APA exposure contained predictive information beyond the prespecified covariates and unweighted score in this dataset, but the magnitude and form of benefit varied by context.

The full standardized model condition number was 4.07–5.40, providing no indication of severe numerical instability. In-sample partial F tests for adding W to C+U were nominally significant in the six primary comparisons (largest P=0.0105), but were not used to claim independent generalization or mechanism. TargetScanHuman Release 8.0 context++ comparison remained pending because a versioned bulk file could not be retrieved and no unofficial substitute was introduced.

### Direct-interaction prioritization replicated in CLEAR-CLIP but not CLASH

In the primary Huh-7.5 CLEAR-CLIP 3′UTR analysis, APA-exposure PU-AUC was 0.671 for miR-16, 0.652 for miR-17, 0.667 for miR-25, 0.551 for miR-92a, and 0.648 for miR-186. The one-sided rank tests remained significant after Benjamini–Hochberg correction for miR-16 (q=2.76×10⁻⁴), miR-17 (q=2.73×10⁻⁷), miR-25 (q=0.0431), and miR-186 (q=0.00525), but not miR-92a (q=0.159). Candidate coverage ranged from 15.6% to 30.3%, reflecting both the canonical-site restriction and interactions outside the modeled 3′UTR universe.

The independent HEK293 CLASH 3′UTR check was weaker and heterogeneous: APA-exposure PU-AUC ranged from 0.504 to 0.601, and none of the five one-sided rank tests passed q<0.05. This was treated as a non-replication of the strong CLEAR-CLIP rank-distribution effect. Protocol, annotation, cell-state, and noncanonical-interaction differences remain plausible explanations but were not tested as post hoc exclusions. Exposure-weighted and unweighted rankings did not show a uniform direction of difference in the direct-interaction data, precluding a general APA-increment claim from this benchmark.

### Reporter-supported targets were strongly prioritized, with miRNA-specific APA increments

In the orthogonally supported miRTarBase set and the at-least-two-context consensus universe, APA-exposure PU-AUC ranged from 0.722 to 0.899 among the six miRNAs with at least five observable positives. All six passed Benjamini–Hochberg q<0.001. The results were: miR-16, PU-AUC 0.830 (25/49 observable, q=2.35×10⁻⁸); miR-17, 0.782 (24/45, q=2.23×10⁻⁶); miR-25, 0.824 (9/25, q=4.48×10⁻⁴); miR-92a, 0.899 (8/29, q=8.76×10⁻⁶); miR-124, 0.802 (34/74, q=4.82×10⁻⁹); and miR-155, 0.722 (23/72, q=1.74×10⁻⁴). miR-186 had only two observable positives and was not interpreted inferentially.

APA-exposure PU-AUC exceeded unweighted-additive PU-AUC for every panel miRNA, but bootstrap evidence for the increment was nonuniform. Differences remained positive after Benjamini–Hochberg correction for miR-16 (+0.048, 95% bootstrap interval 0.005–0.087, q=0.0340), miR-25 (+0.070, 0.017–0.128, q=0.0150), and miR-124 (+0.060, 0.026–0.093, q=0.00120). Intervals crossed zero for miR-17, miR-92a, and miR-155; miR-186 was not tested. Sensitivity analyses using all reporter records preserved positive ranking performance, whereas the all-four-context requirement substantially reduced coverage for several miRNAs.

### Construct-level examples occupied the extreme upper tail

All four full-text-verified wild-type/cognate-mutant pairs were observable and fell in the top APA-exposure decile (mean percentile 0.9673; median 0.9733). Including construct-corroborated miR-155–ETS1 yielded five of five pairs in the top decile (mean 0.9668; median 0.9649). These results provide concrete biological examples of high-ranked validated targets but, because the audit followed inspection of the broad reporter benchmark and involved only four or five pairs, they were not used as an unbiased performance estimate.

### The locked holdout supported ranking ability conditional on 57.1% coverage

Of 33 locked pairs, 21 met every construct-level inclusion rule and 12 were excluded for protocol or provenance reasons. Twelve of the 21 strict inclusions were observable in the frozen rank table (coverage 57.1%); nine remained missing. Among observable pairs, the mean APA-exposure percentile was 0.7847 and the median was 0.8185. Four of 12 pairs were in the top decile and seven of 12 were in the top quintile. The median paired APA-minus-unweighted percentile was +0.0651.

In the explicitly post-unblinding matched-null analysis, the observed pooled mean percentile exceeded the theoretical null mean (0.7847 versus 0.4992; Holm P=0.00080). The corresponding Holm-adjusted P values were 0.01164 for the median percentile, 0.02534 for the top-decile fraction, and 0.01164 for the top-quintile fraction. Equal-miRNA cluster analysis estimated a mean percentile of 0.7915 with a 95% bootstrap interval of 0.7030–0.8847. In contrast, the cluster-equal APA-minus-unweighted estimate was +0.0423 and its interval crossed zero (−0.0136 to 0.0899). The holdout therefore supports ranking of validated targets conditional on observability, but does not independently establish added value from APA weighting.

### Cross-benchmark interpretation

Three results converge on the narrow primary conclusion: reporter-supported targets were strongly prioritized, the rank-independent holdout placed strict construct-validated targets above the middle of their candidate universes, and exposure-aware scores added predictive information in the matched perturbation dataset. The CLEAR-CLIP data provided additional support for direct-target prioritization. However, CLASH did not reproduce the strong rank-distribution result, candidate coverage was incomplete, and evidence for APA-over-unweighted gain varied by miRNA and context. The appropriate conclusion is therefore that APAmiRank is a reproducible candidate-prioritization framework with condition-dependent evidence for added APA information—not a universally superior target predictor.

## Figure legend

**Figure X. Rank-independent locked luciferase holdout and exploratory uncertainty analysis.** (A) APA-exposure percentile ranks of the 12 observable strict inclusions. The reference lines indicate the 50th, 80th, and 90th percentiles. (B) Null distribution of the pooled mean percentile from 50,000 theoretical miRNA-matched rank-uniform simulations. The vertical line marks the observed mean of 0.7847 (Holm-adjusted exploratory P=0.00080), and the dashed line marks the null mean. (C) Paired comparison of APA-exposure and unweighted-additive percentiles for the observable pairs. The diagonal denotes equal ranks; the observed median paired difference was +0.0651. (D) Screening and rank coverage. Of 33 rank-free locked candidates, 21 satisfied the frozen human 3′UTR wild-type/cognate-mutant reporter criteria and 12 were observable in the consensus rank table. Protocol or provenance exclusions are not biological negatives, and the nine unobservable strict inclusions were not imputed. Panel B is post-unblinding and exploratory; all rank summaries are conditional on 57.1% observability. The separate equal-miRNA cluster bootstrap estimated a mean APA percentile of 0.7915 (95% interval 0.7030–0.8847) and an APA-minus-unweighted difference of +0.0423 (−0.0136 to 0.0899).

## Limitations and claim boundary

1. **Positive–unlabeled outcomes.** Non-annotated genes were not proven negatives. PU-AUC and enrichment quantify prioritization relative to an unlabeled canonical-site universe and must not be described as specificity or diagnostic discrimination.
2. **Restricted eligibility and missingness.** APAmiRank evaluates genes with enumerated canonical 3′UTR sites. Noncanonical sites, non-3′UTR interactions, genes absent from the historical annotation, and unobservable holdout pairs are outside the evaluated universe. Locked-holdout coverage was 57.1%, and inference is conditional on the 12 observable pairs.
3. **Context mismatch.** The reporter literature spans heterogeneous biological contexts. Huh-7.5 CLEAR-CLIP was only near-matched to Huh7 APA profiles. A consensus across four APA contexts reduces context specificity rather than proving it.
4. **APA uncertainty.** GSE52527 provides a single processed 3P-seq library per cell type, so biological-replicate uncertainty in APA usage is not estimated. Historical hg19/RefSeq mapping and internal-priming filters may also affect isoform assignment.
5. **Selection bias.** Reporter and luciferase experiments are often selected after canonical seed prediction, which can favor a canonical-site ranking method. The four-pair pilot was reviewed after the broad benchmark and is descriptive only.
6. **Holdout design.** Candidate generation and construct screening were rank independent, but the study was not a fully prospective blinded trial because the broader ranking resource existed. The matched-null and cluster-bootstrap procedures were designed after unblinding and remain exploratory.
7. **Small cluster count and theoretical null.** Only six miRNA clusters contributed to holdout uncertainty estimates. The theoretical rank-uniform null does not reproduce empirical score ties or correct missingness mechanisms.
8. **Non-replication.** The strong CLEAR-CLIP effect did not replicate in the independent CLASH 3′UTR analysis. This result constrains generalization across direct-interaction technologies and contexts.
9. **Comparator gap.** TargetScanHuman Release 8.0 data were not obtained during the checkpoint. TargetScan already incorporates 3′UTR isoform profiles; without a version-locked comparison, neither superiority to TargetScan nor isolation of APA-specific benefit relative to that method can be claimed.
10. **No causal or clinical interpretation.** The scores do not establish direct endogenous binding, functional mediation, causal cell-type specificity, therapeutic response, or clinical utility.

### Defensible manuscript claim

> APAmiRank reproducibly prioritizes experimentally supported miRNA targets within its canonical 3′UTR site universe. Functional perturbation and reporter analyses indicate that APA exposure can add predictive information beyond an otherwise matched unweighted architecture score, but this increment is miRNA- and context-dependent and was not independently resolved in the locked holdout.

### Claims not supported by the current evidence

- APA weighting improves prediction for every miRNA or cell type.
- APAmiRank is superior to TargetScan or other state-of-the-art predictors.
- An unannotated or low-ranked gene is a true non-target.
- The reported PU-AUCs estimate specificity or clinical discrimination.
- The locked holdout was preregistered, fully prospective, or fully blinded.
- APAmiRank scores prove endogenous interaction, causal regulation, or therapeutic relevance.

## References

1. Nam J-W, Rissland OS, Koppstein D, et al. Global analyses of the effect of different cellular contexts on microRNA targeting. *Molecular Cell*. 2014;53:1031–1043. https://doi.org/10.1016/j.molcel.2014.02.013
2. Helwak A, Kudla G, Dudnakova T, Tollervey D. Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding. *Cell*. 2013;153:654–665. https://doi.org/10.1016/j.cell.2013.03.043
3. Moore MJ, Scheel TKH, Luna JM, et al. miRNA–target chimeras reveal miRNA 3′-end pairing as a major determinant of Argonaute target specificity. *Nature Communications*. 2015;6:8864. https://doi.org/10.1038/ncomms9864
4. Huang H-Y, Lin Y-C-D, Cui S, et al. miRTarBase update 2025: an informative resource for experimentally validated miRNA–target interactions. *Nucleic Acids Research*. 2025;53:D147–D156. https://doi.org/10.1093/nar/gkae1072
5. Ha KCH, Blencowe BJ, Morris Q. QAPA: a new method for the systematic analysis of alternative polyadenylation from RNA-seq data. *Genome Biology*. 2018;19:45. https://doi.org/10.1186/s13059-018-1414-4

## Data and code availability statement (working draft)

All analysis code, frozen protocols, source-accession manifests, checksums, random seeds, pair-level outputs, and publication figures described here are included in the APAmiRank v0.2.18 research package. Public-source datasets are not redistributed; scripts retrieve or transform the corresponding GEO and miRTarBase records. The release currently remains research software pending selection of a final open-source license.
