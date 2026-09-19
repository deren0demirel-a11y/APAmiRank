# Verification record for v0.2.18

Verification date: 2026-09-10

The following checks passed in the packaging environment:

- mature-miRNA normalization and dynamic canonical seed construction;
- different miRNA sequences generate different target motifs;
- canonical-site enumeration in plus- and minus-strand synthetic loci;
- common versus distal-extension site classification;
- complete architecture-mode run from YAML configuration;
- strict validation of complete condition-specific isoform fractions;
- expected common-site exposure across proximal and distal contexts;
- expected condition dependence of distal-extension site exposure;
- condition-specific additive gene scoring and pairwise contrasts;
- preservation of distinct matched mock groups in GSE52530 expression parsing;
- strand-aware conversion of 3P-tag cluster boundaries to PAS coordinates;
- bounded-span cross-cell-line PAS harmonization;
- exact sequence length validation after UCSC twoBit extraction;
- deterministic site and gene ranking outputs;
- preservation of primary ranks after post hoc gene annotation;
- seven-feature thermodynamic scoring schema using synthetic evidence;
- RNAup energy/coordinate parser test;
- Python byte-code compilation;
- TOML, YAML, and CFF parsing.
- CLASH transcript-coordinate classification against hg19 `ensGene`;
- positive--unlabeled CLEAR-CLIP/CLASH rank and enrichment evaluation;
- Benjamini--Hochberg adjustment within prespecified five-miRNA test families;
- publication-figure generation in PNG, SVG, and PDF formats.
- validation of strict luciferase eligibility fields and separation of the
  full-text and construct-corroborated scopes;
- deterministic strict-luciferase rank joins and descriptive summaries.
- deterministic rank-independent holdout selection and exclusion of prior
  pilot-audit pairs;
- fail-closed rejection of pending or checksum-mismatched screening forms.
- fielded PubMed query construction, pair-specific article mapping, and
  assay-language candidate prioritization;
- extraction of reporter-plus-design passages from representative PMC XML.
- deterministic miRNA-matched rank-uniform null simulation with sampling
  without replacement inside each ranking universe;
- Holm adjustment of the four secondary ranking-enrichment metrics;
- complete-miRNA-cluster bootstrap and publication-figure generation.

Thirty-three pytest tests passed in the v0.2.18 development environment. The packaging
environment did not contain RNAhybrid, RNAup, or Snakemake executables.
Consequently, the external-binary thermodynamic execution and Snakemake entry
point were not run here. Their parsers, scoring path, pinned Conda dependencies,
and GitHub Actions workflow are included. Before a public release, create the
provided Conda environment and run:

```bash
pytest -q
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.example.yaml \
  --cores 4
```

The included example intentionally uses a synthetic miRNA and synthetic loci;
it is a software test and has no biological interpretation.

The GSE52527/GSE52530 processed files and hg19 reference were downloaded again
from their recorded official URLs and local checksums were regenerated. The
architecture benchmark counts and R² values were reproduced. The covariate
evaluation completed 100 repeated 10-fold splits for four expression thresholds
and eight cell-line/miRNA combinations. Full-model condition numbers were
4.07--5.40 at the primary RPKM >=1 threshold.

TargetScanHuman Release 8.0 was not included in the executed comparison because
the official bulk endpoint returned HTTP 502 and the documented Figshare mirror
returned HTTP 403 during verification. No substitute score source was used.

The direct-interaction benchmark downloaded GSE50452 CLASH, GSE73057
CLEAR-CLIP, UCSC hg19 `ensGene`, and miRBase mature sequences from the URLs
recorded in `benchmarks/direct_interactions/README.md`; SHA-256 checksums were
recorded. Ten architecture-mode panel runs completed. The primary CLEAR-CLIP
3′UTR benchmark produced top-decile enrichment above twofold for all five
miRNAs; four of five APA-exposure rank tests passed BH q<0.05. The independent
CLASH 3′UTR analysis was weaker and no APA-exposure test passed BH q<0.05.

The reporter benchmark downloaded a dated human Reporter Assay search export
from the official miRTarBase portal and recorded its checksum. Twenty-eight
architecture-mode runs completed across seven miRNAs and four APA contexts.
In the primary Reporter Assay plus Western blot/qPCR set, six miRNAs with at
least five observable positives had consensus APA-exposure PU-AUCs of
0.722--0.899 and all passed BH q<0.001. APA-exposure AUC exceeded unweighted
additive AUC for all seven panel miRNAs, but paired bootstrap confidence
intervals excluded zero after BH correction only for miR-16, miR-25, and
miR-124. miR-186 had only two observable primary positives and was suppressed
from incremental inference.

The construct-level luciferase audit retained four full-text-verified exact
miRNA–target pairs with WT repression and loss/weakening after 3′UTR-site
mutation. All four were observable and fell in the top APA-exposure decile
(mean percentile 0.9673; median 0.9733). Adding one separately labelled
construct-corroborated ETS1 record yielded five of five in the top decile.
Because literature candidates were inspected after the broader reporter ranks
were available and the sample is very small, these are explicitly descriptive
pilot results; no inferential test was run.

The locked follow-up manifest contains 33 previously unaudited pairs selected
without any APAmiRank score or rank field: five for six panel miRNAs and all
three eligible miR-186-5p pairs. Regeneration reproduced the byte-identical
manifest and SHA-256
`e5c28f59a84ddfaf63d9f1d18c69c24ba730cc956c5a0628c82131975410b125`.
The evaluation gate was first tested against the initialized form and correctly
refused to join ranks while decisions remained pending. PubMed discovery
returned 298 pair-article candidates across 299 unique PMIDs. Original-paper
review is now complete for all 33 pairs: 21 included and 12 excluded. The final
33-row, 25-column screening form was frozen before rank access with SHA-256
`8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5`.
The gate then accepted the checksum and joined only the 21 strictly included
pairs to the prespecified consensus-rank table. Twelve pairs were observable
(coverage 0.5714). Their mean APA-exposure percentile was 0.7847 and median was
0.8185; four of 12 were in the top decile and seven of 12 in the top quintile.
The median paired APA-exposure minus unweighted-additive percentile was 0.0651.
The pair-level and summary outputs have SHA-256 values
`7a7459ae95fe210488a6d63e80cdcee625ba6eed9b80ab6d083f8860d22477bf`
and `a16806c7712ea8192010acfec32376d760fc83c5b0de19bce1c7e0c90db46a1a`,
respectively. These are descriptive holdout metrics; the 57.1% candidate
coverage is an explicit limitation and missing pairs were not assigned scores.

After the primary descriptive result was available, a secondary exploratory
uncertainty analysis ran 50,000 miRNA-matched rank-uniform null iterations and
50,000 complete-cluster bootstrap iterations. It is explicitly labelled
post-unblinding and not prespecified. The pooled mean APA percentile was 0.7847
versus a simulated null mean of 0.4992 (one-sided Monte Carlo P=0.00020;
Holm-adjusted P=0.00080). Holm-adjusted P values were 0.01164 for the median,
0.02534 for the top-decile fraction, and 0.01164 for the top-quintile fraction.
The miRNA-equal mean percentile was 0.7915 (six-cluster bootstrap interval
0.7030–0.8847). The paired APA-minus-unweighted cluster-equal estimate was
0.0423 with interval −0.0136 to 0.0899, so incremental APA benefit was not
independently established. The null is conditional on 12 observable pairs and
does not reproduce score-tie multiplicities; only six miRNA clusters contribute.

Secondary output SHA-256 values:

- `holdout_uncertainty_summary.tsv`: `40fd9fb563efd052674e1c90e4a87416b0c8e27385c8b95062e21f3905d005df`
- `holdout_mirna_cluster_summary.tsv`: `d11a4c2eafd6939ddd34f784470a0c97c38c17060272a668bb9d36431fbf4d86`
- `holdout_rank_uniform_null.tsv.gz`: `c9214c36df8503f8e2d5f8215261c5bb347de3364103fef791c435dafe2e17f0`
- `holdout_cluster_bootstrap.tsv.gz`: `446dcde1bcfe147101c852be451de8fa3ac06e42242b8d91caf30dccd5ad884d`
- figure PDF/PNG/SVG: `1a8e3edd1bd8f0e30af6601867fe067ec65c1a99d8fe3caf0ca6e9788b68fc37`,
  `6e3520f84f99aa46bfb84c7e6259b73aab06da9062bdafc0bb85885da228a1a0`, and
  `985307774c335fed79f81ea342413d0fc24e2d411064c4f65bc1aa2111c90727`.

The manuscript-ready validation narrative was checked against the benchmark
result records. Its Word edition was rendered to eight PNG pages and every
page was inspected for clipping, overlap, table wrapping, figure readability,
panel/legend correspondence, page breaks, and reference numbering. Final
manuscript SHA-256 values:

- `manuscript/VALIDATION_METHODS_RESULTS.md`: `81306f293986f312b32ab332504da5a1b6a3c69be23551009e80f58e5cdebbcd`
- `manuscript/APAmiRank_Validation_Methods_Results_v0.2.18.docx`: `2b8df258cc064423abc746c22cea51a5f66dcfa98bccc5e1e56e160e875f35e7`
