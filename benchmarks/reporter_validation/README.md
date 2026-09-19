# Reporter-assay validation benchmark

This benchmark tests whether APAmiRank prioritizes human miRNA--target pairs
annotated by miRTarBase as supported by a `Reporter Assay`. It is a
positive--unlabeled analysis: all other candidate genes remain unlabeled and
are not called true negatives.

## Evidence snapshot

The source CSV was exported on 2026-09-08 from the official miRTarBase search
interface with `mode=method`, `species=hsa`, and `methods=reporter_assay`:

https://awi.cuhk.edu.cn/miRTarBase/search/results/download/?mode=method&species=hsa&methods=reporter_assay

The download endpoint is not version-tagged. It is therefore described as a
dated portal snapshot rather than asserted to be an immutable v10 file. The
latest peer-reviewed release at the time of analysis is miRTarBase 10.0:
https://doi.org/10.1093/nar/gkae1072.

The snapshot contains 7,790 rows, of which 7,669 are exact human-miRNA to
human-target records. The seven-miRNA panel contributes 513 unique human--human
pairs after harmonization. miRTarBase reported 24,530 reporter-assay-supported
MTIs across its release-wide collection; that database-wide count and this
human filtered export are not interchangeable.

`Reporter Assay` is retained as the database's evidence category. The export
does not establish that every record used a luciferase construct, a mutant
site, the complete 3'UTR, or the same cell line as the APA data.

## Frozen panel and evidence sets

The panel contains miR-16-5p, miR-17-5p, miR-25-3p, miR-92a-3p, miR-124-3p,
miR-155-5p, and miR-186-5p. The primary positive set requires Reporter Assay
plus at least one orthogonal low-throughput assay (`Western Blot` or `qPCR`).
Sensitivity sets retain all reporter records, Tier A/B records, or Tier A only.

Because the reporter literature spans many biological contexts and the CSV
does not provide a usable cell line for every record, no interaction is called
cell-matched. The primary consensus universe requires a gene to be rankable in
at least two of the four APA contexts (HeLa, HEK293, Huh7, IMR90), then averages
its available percentile ranks and re-ranks the consensus score. Requiring all
four contexts is retained as a strict sensitivity analysis. Individual cell
contexts are reported separately.

## Reproduction

```bash
bash download_reporter_snapshot.sh raw
python prepare_reporter_evidence.py \
  --input raw/miRTarBase_human_reporter_2026-09-08.csv \
  --output derived/reporter_evidence.tsv

PYTHONPATH=../../src python run_panel.py \
  --derived-dir ../gse52531/derived \
  --output-root runs

PYTHONPATH=../../src python evaluate_reporter_validation.py \
  --evidence derived/reporter_evidence.tsv \
  --runs runs \
  --output-dir evaluation
```

Install `.[benchmark]` first. The evaluation reports candidate coverage,
positive--unlabeled AUC, one-sided Mann--Whitney tests, BH correction across the
seven-miRNA family, and top-decile/top-quintile hypergeometric enrichment.

The incremental APA analysis compares consensus APA-exposure AUC with the
unweighted-additive AUC on the same genes. A paired stratified bootstrap samples
positive and unlabeled genes with replacement for 5,000 iterations. Results
with fewer than five observable positive genes are suppressed.

## Interpretation limits

Reporter-positive interactions are affected by literature-selection and
construct-selection bias. Many reporters were chosen after canonical seed
prediction, which can favor a canonical-site ranking method. Gene symbols and
miRNA names are database annotations rather than independent remapping of every
paper. These results support prioritization of previously validated targets;
they do not estimate specificity, false-positive rate, or performance on novel
targets.
