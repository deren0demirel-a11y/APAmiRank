# Direct-interaction benchmark

This benchmark asks whether APAmiRank prioritizes genes supported by ligation-based,
direct miRNA--mRNA interaction assays. It is a **positive--unlabeled** evaluation:
absence from CLASH or CLEAR-CLIP is not interpreted as a negative interaction.

## Sources and scope

| Dataset | Assay / cell line | Primary use | Source |
|---|---|---|---|
| GSE73057 | CLEAR-CLIP, Huh-7.5 | Primary gene-level benchmark restricted to annotated 3'UTR chimeras | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73057 |
| GSE50452 | CLASH, HEK293 | Independent replication; transcript coordinates are annotated as 5'UTR/CDS/3'UTR with hg19 `ensGene` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50452 |
| hg19 `ensGene` | UCSC transcript models | CLASH transcript-coordinate region assignment | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz |
| miRBase mature sequences | miRNA sequence reference | Panel sequences and identifier audit | https://www.mirbase.org/download/mature.fa |

Primary publications: Helwak et al., *Cell* 2013,
https://doi.org/10.1016/j.cell.2013.03.043; Moore et al., *Nature
Communications* 2015, https://doi.org/10.1038/ncomms9864.

The predeclared panel contains miR-16-5p, miR-17-5p, miR-25-3p,
miR-92a-3p, and miR-186-5p because all five have usable direct-interaction
support in both datasets. The earlier miR-124/miR-155 perturbation benchmark
cannot be reused as a well-powered direct benchmark: GSE50452 contains only two
miR-124 chimeras (both mapped to CDS), while GSE73057 contains one annotated
3'UTR miR-155 chimera and that gene lacks a canonical site in the APAmiRank
candidate universe.

Huh-7.5 (CLEAR-CLIP) and Huh7 (APA profiles) are related but not identical cell
lines; this comparison is therefore near-matched rather than exact-matched.
HEK293 is cell-line matched between CLASH and the APA profiles.

## Reproduction

From this directory:

```bash
bash download_data.sh raw
python prepare_direct_interactions.py \
  --clash-dir raw/GSE50452 \
  --clear-workbook raw/GSE73057_Huh7_miRNA_chimera_interactions.xlsx \
  --ensgene raw/hg19_ensGene.txt.gz \
  --output derived/direct_interactions.tsv

PYTHONPATH=../../src python run_panel.py \
  --derived-dir ../gse52531/derived \
  --output-root runs

PYTHONPATH=../../src python evaluate_direct_interactions.py \
  --interactions derived/direct_interactions.tsv \
  --runs runs \
  --output-dir evaluation
```

Install the benchmark extra (`pip install -e '.[benchmark]'`) before running.

## End points

- Percentile rank among genes containing at least one canonical 3'UTR site.
- Positive--unlabeled AUC (direct positives versus the unlabeled candidate set),
  reported descriptively and never called specificity.
- One-sided Mann--Whitney test, with Benjamini--Hochberg correction across the
  five miRNAs within each dataset / region / rank-metric family.
- Top-10% and top-20% enrichment with hypergeometric P values.
- `source_recall` uses all unique direct-interaction genes as denominator;
  `observable_recall` uses only direct genes that enter APAmiRank's canonical-site
  candidate universe. Candidate coverage is reported explicitly.

Four ranks are retained: best-site, site-count opportunity-adjusted, unweighted
additive, and APA-exposure-weighted additive. Comparing the last two is a paired
descriptive sensitivity analysis, not an isolated causal estimate of APA benefit.

## Interpretation boundary

The CLEAR-CLIP 3'UTR analysis is the primary result. CLASH 3'UTR and all-mRNA
analyses are independent secondary checks. Differences in protocol, genome
annotation, cell state, noncanonical binding, transcript expression, and the
canonical-site eligibility rule can all reduce overlap. No absent interaction is
treated as a confirmed negative.
