# Changelog

## 0.2.18 — 2026-09-10

- Added a manuscript-ready English validation narrative integrating the
  perturbation, direct-interaction, reporter, construct-level, and locked-
  holdout evidence layers.
- Separated evidence for target prioritization from evidence for incremental
  benefit attributable to APA exposure and retained the CLASH non-replication.
- Added a corrected four-panel figure legend, explicit limitations, a
  defensible claim statement, unsupported-claim boundaries, and references.
- Added a styled DOCX edition with the locked-holdout publication figure and
  completed an eight-page render-based visual inspection.
- Preserved the unresolved TargetScan Release 8.0 comparison as an explicit
  comparator gap rather than substituting an unofficial score source.
- Re-ran the complete regression suite; 33 tests pass.

## 0.2.17 — 2026-09-10

- Added a separately labelled post-unblinding secondary uncertainty analysis
  without altering the locked manifest, screening decisions or primary join.
- Ran 50,000 theoretical rank-uniform null draws matched by miRNA, ranking-
  universe size and observable-pair count; adjusted four tests by Holm's method.
- Found conditional enrichment of observable validated targets: mean percentile
  0.7847 versus null mean 0.4992 (Holm P=0.00080).
- Added a 50,000-iteration equal-miRNA cluster bootstrap. The incremental
  APA-minus-unweighted interval crossed zero (−0.0136 to 0.0899), so no
  independent APA-increment claim is made from this holdout.
- Documented that the secondary null is theoretical, does not preserve tie
  multiplicities, is conditional on 57.1% coverage, and uses only six clusters.
- Added deterministic secondary-analysis tests and a four-panel publication
  figure in PNG, SVG and PDF; 33 tests pass.

## 0.2.16 — 2026-09-10

- Completed all 33 rank-independent locked luciferase reviews: 21 included and
  12 excluded under the prespecified construct-level criteria.
- Included miR-92a-3p–CPEB2 and –CDH1 after original-paper verification of human
  3′UTR WT repression and loss of response in cognate mutants.
- Excluded NRF1 fail-closed because the mouse-model paper did not establish a
  human reporter-insert provenance; excluded STAT3 under `WRONG_TARGET` because
  its linked reporter tested RECK; excluded MTO1 under `FULL_TEXT_UNAVAILABLE`.
- Froze the final screening form by SHA-256, opened the evaluation gate, and
  joined only the 21 strict inclusions to the prespecified rank table.
- Recorded 12 observable pairs (57.1% coverage), mean/median APA percentiles of
  0.7847/0.8185, top-decile fraction 0.3333, top-quintile fraction 0.5833, and
  median APA-minus-unweighted percentile 0.0651.
- Added auditable pair-level and summary outputs and documented that missing
  candidates were not assigned or imputed scores.
- Added regression tests for a conforming completed gate and rejection of an
  included record with a non-human biological field.

## 0.2.15 — 2026-09-09

- Completed original-paper screening for five locked miR-25-3p pairs.
- Included FBXW7, MDM2 and TP53 after verifying WT 3′UTR repression and loss
  or material weakening of that response in cognate-site mutants.
- Excluded CCL26 and TCEAL1 under `WRONG_TARGET`: the linked source papers test
  TSC1 and BIM, respectively, rather than the genes named by those records.
- Kept the provenance exclusions distinct from biological negative results.
- Kept the evaluation gate closed with five of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.14 — 2026-09-09

- Completed original-paper screening for all three locked miR-186-5p pairs.
- Excluded AKAP12 and CSNK2A1 fail-closed under `FULL_TEXT_UNAVAILABLE` because
  the required cognate WT/mutant evidence could not be verified from the
  accessible original-source material.
- Excluded MOB1A under `NOT_3UTR`: its open original article used the 5′
  flanking sequence of the MOB1 promoter rather than a target-gene 3′UTR.
- Explicitly rejected the unrelated, retracted circ_0008360/CCND2 article as
  evidence for AKAP12.
- Kept the evaluation gate closed with 10 of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.13 — 2026-09-09

- Completed the miR-124-3p–ROCK2 original-paper review.
- Recorded the article's 2021 correction and excluded the pair fail-closed
  under `FULL_TEXT_UNAVAILABLE`, because the corrected Figure 6 and the
  WT/cognate-mutant construct details could not be inspected reliably.
- Preserved the distinction between evidence inaccessibility and biological
  falsification: the accessible abstract still reports ROCK2 3′UTR targeting
  and endogenous mRNA/protein suppression.
- Kept the evaluation gate closed with 13 of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.12 — 2026-09-09

- Completed original-paper screening for miR-17-5p–ITGB8 and STAT3.
- Included STAT3 on the basis of an exact-miR-17, full-length human 3′UTR
  reporter with single- and double-site mutant rescue plus endogenous protein
  support.
- Excluded ITGB8 under `OTHER_WITH_NOTE`: the database-linked article tested
  WT/mutant reporters for ADAR1 only, while ITGB8 received qPCR-only follow-up
  that the authors described as questionable.
- Kept the evaluation gate closed with 14 of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.11 — 2026-09-09

- Completed original-paper screening for miR-17-5p–JAK1, RUNX1, and CDKN1A;
  all three satisfy the human-target, cognate WT/mutant 3′UTR criterion and
  have orthogonal endogenous support.
- Anchored RUNX1 to the original 2007 AML1 study and documented why its
  individual anti-miR-17-5p experiment supports exact-arm attribution despite
  the shared miR-17-family seed site.
- Anchored CDKN1A to the stronger two-site-mutant and mRNA-protector study
  rather than a WT-only reporter paper.
- Kept the evaluation gate closed with 16 of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.10 — 2026-09-09

- Completed original-paper screening for miR-16-5p–WEE1, YAP1, KRAS, and
  PPM1D; all four satisfy the cognate WT/mutant-or-deletion 3′UTR criterion
  and have orthogonal endogenous support.
- Anchored YAP1 to the explicit bidirectional psiCHECK-2 experiment in PMID
  34632051 and recorded its exact cognate-site substitution.
- Recorded the 2021 KRAS author correction and verified that it changes Figure
  4D invasion images rather than the Figure 3E reporter evidence used here.
- Kept the evaluation gate closed with 19 of 33 decisions pending; no
  APAmiRank rank, score, percentile, or outcome table was opened or joined.
- Re-ran structural validation and the complete regression suite.

## 0.2.9 — 2026-09-09

- Completed original-paper screening for miR-155-5p–APC, HBP1, and CLDN1,
  each with a cognate WT/mutant 3′UTR reporter comparison and orthogonal
  endogenous support.
- Completed original-paper screening for miR-16-5p–CDS2 using the entire human
  CDS2 3′UTR, a seed-complementary mutant, and near-complete mutant rescue.
- Preserved publication-integrity provenance and explicitly excluded an
  unrelated retracted HBP1 paper from the evidence base.
- Kept the evaluation gate closed with 23 of 33 decisions pending and did not
  access or join APAmiRank rank, score, percentile, or outcome tables.
- Re-ran the complete regression suite; all 28 tests passed.

## 0.2.8 — 2026-09-09

- Completed original-paper screening for miR-124-3p–MAPK14 using the
  full-length 3′UTR pCI-FLuc reporter, cognate seed mutant, mutant rescue, and
  endogenous p38α evidence in human cells.
- Completed original-paper screening for miR-124-3p–AHR using an independent
  human intestinal-cell study with WT and seed-mutant AHR 3′UTR reporters and
  orthogonal endogenous-protein and patient-tissue support.
- Recorded the 2021 ROCK2 source-article correction and retained ROCK2 as
  pending until its construct-level impact is fully resolved.
- Kept the evaluation gate closed with 29 of 33 decisions pending and did not
  access or join APAmiRank rank, score, percentile, or outcome tables.
- Re-ran the complete regression suite; all 28 tests passed.

## 0.2.7 — 2026-09-08

- Added reproducible PubMed candidate discovery for every locked miRNA-target
  pair, yielding 298 pair-article candidates across 299 unique PMIDs.
- Added a non-decisional PMC full-text passage extractor and parser tests to
  accelerate manual WT/mutant 3′UTR screening without weakening eligibility.
- Completed original-paper review for the first two locked pairs: miR-124-3p–
  GRB2 and miR-124-3p–PLEC both satisfy the strict inclusion criteria.
- Inspected the linked PLEC correction and documented that it replaces a
  mislabeled Fig. 4A image rather than altering the Fig. 3 reporter evidence.
- Kept the evaluation gate closed with 31 of 33 screening decisions pending;
  no APAmiRank rank or score was joined during literature review.
- Expanded the regression suite to 28 tests.

## 0.2.6 — 2026-09-08

- Locked a deterministic 33-pair luciferase holdout selected without reading
  APAmiRank score or rank columns and excluding every v0.2.5 pilot-audit pair.
- Added a checksum-protected candidate manifest, rank-free screening form,
  original-paper eligibility protocol, and explicit correction/retraction checks.
- Added a fail-closed evaluation gate that blocks rank joins until all screening
  decisions are complete, provenance fields are populated, and the completed
  form's SHA-256 matches its lock.
- Expanded the regression suite to 23 tests.

## 0.2.5 — 2026-09-08

- Added a transparent manual audit for exact-miRNA WT/mutant 3′UTR luciferase
  evidence, retaining included, excluded, and unresolved records.
- Added separate full-text-verified and construct-corroborated sensitivity
  scopes and deterministic joins to the frozen reporter-panel ranks.
- Added deliberately descriptive small-sample summaries, regression tests, and
  explicit post-selection-bias safeguards; no inferential test is emitted.

## 0.2.4 — 2026-09-08

- Added a dated official-miRTarBase Reporter Assay snapshot workflow and a
  seven-miRNA positive--unlabeled reporter-validation benchmark.
- Added a primary Reporter Assay plus Western blot/qPCR evidence definition,
  four-context consensus ranks, strict-universe and evidence-tier sensitivities,
  and per-cell-context summaries.
- Added paired stratified bootstrap estimates of APA-exposure AUC increment
  over unweighted-additive ranking, with minimum positive-count safeguards.
- Added publication figures, reporter benchmark tests, source checksums, and
  explicit construct/cell-context/publication-bias limitations.

## 0.2.3 — 2026-09-08

- Added a positive--unlabeled direct-interaction benchmark using GSE73057
  CLEAR-CLIP and GSE50452 CLASH.
- Added a prespecified five-miRNA panel, transcript-coordinate 3′UTR mapping,
  candidate-coverage reporting, rank-distribution tests, top-decile enrichment,
  and Benjamini--Hochberg correction.
- Added reproducible download, harmonization, panel-run, evaluation, and figure
  scripts plus direct-interaction regression tests.
- Expanded benchmark dependencies to include pandas and openpyxl.

## 0.2.2 — 2026-09-08

- Added prespecified covariate-adjusted GSE52531 models controlling baseline
  expression, benchmark 3′UTR length, canonical site count, and mean local AU.
- Added deterministic repeated 10-fold cross-validation, nested-model tests,
  standardized coefficients, and numerical-condition diagnostics.
- Added gene-level covariate input exports and a documented checkpoint report.
- Added a `benchmark` optional dependency group for NumPy and SciPy.
- Added covariate-model unit tests.

## 0.2.1 - 2026-09-07

- Added conservative, strand-aware mapping of GSE52527 3P-tag clusters to
  contiguous RefSeq terminal-exon 3′UTRs on hg19.
- Added per-cell-line and harmonized four-cell-line PAS/usage preparation.
- Added coordinate-exact strand-correct extraction from the UCSC hg19 2bit
  reference.
- Added reproducible architecture benchmark runners and outcome evaluation
  against matched GSE52530 miRNA perturbation data.
- Added paired bootstrap uncertainty and cognate-versus-non-cognate AIR
  comparisons.
- Expanded the test suite to 11 tests.

## 0.2.0 - 2026-09-07

- Added optional condition-specific isoform usage input.
- Added strict locus/condition completeness and fraction-sum validation.
- Added site exposure and exposure-weighted binding-evidence scores.
- Added additive condition-specific gene rankings and pairwise contrasts.
- Added synthetic proximal-versus-distal usage example and validation tests.
- Added a reproducible GSE52527/GSE52530 downloader and a header-driven,
  replicate-preserving expression converter for the proof-of-concept benchmark.
- Retained separation of binding evidence, APA exposure, and post hoc outcome
  annotations.

## 0.1.0 - 2026-09-01

- First generic APAmiRank release candidate.
- Generates canonical seed motifs from any supplied miRNA sequence.
- Represents common, intermediate, and distal-extension site availability
  across ordered alternative 3′UTR isoforms.
- Separates primary binding ranking from external site annotations and post hoc
  expression/rescue annotations.
- Provides architecture-only and full thermodynamic scoring modes.
- Removes miR-548ah-3p-, VIM-, cell-line-, and contrast-specific assumptions.
- Adds configuration, synthetic example, unit tests, CI, citation metadata, and
  release documentation.
