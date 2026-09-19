# GSE52531 proof-of-concept benchmark

This benchmark combines the two relevant subseries of GSE52531:

- **GSE52527:** one processed human 3P-seq poly(A)-site BED file for each of
  HeLa, HEK293, Huh7, and IMR90 (genome build hg19).
- **GSE52530:** RNA-seq after miR-124, miR-155, or mock transfection in the
  same four cell-line contexts. The series contains two biological replicates
  per transfected miRNA and corresponding mock samples; Huh7 has four mock
  libraries in the series.

Run `bash download_processed.sh` to retrieve the four 3P-seq files and both
the unnormalized and quantile-normalized expression matrices for each cell
line, then write `raw/SHA256SUMS`. Public source data are intentionally excluded
from the repository.

Use the quantile-normalized (`*.expData.qn.txt.gz`) matrices for the primary
benchmark. Convert one matrix to explicit replicate-level long format with:

```bash
python benchmarks/gse52531/prepare_expression.py \
  benchmarks/gse52531/raw/GSE52530_HeLa.expData.qn.txt.gz \
  benchmarks/gse52531/derived/GSE52530_HeLa.expression_long.tsv \
  --cell-line HeLa
```

The converter deliberately does not calculate fold changes or add a
pseudocount. It also preserves the separate miR-124 and miR-155 mock pairs in
Huh7 (`mock_miR-124` and `mock_miR-155`) instead of pooling them. Those analysis
choices and comparison pairs must be prespecified in the benchmark.

## Planned benchmark comparisons

1. Longest annotated 3′UTR with no usage weighting.
2. Dominant 3P-seq isoform only.
3. APAmiRank architecture score with annotation-only site availability.
4. APAmiRank condition-specific exposure weighting.
5. TargetScan weighted context++ score, using an explicitly recorded release.

Primary outcomes will be the association between target score and miRNA-induced
RNA log2 fold change, top-k repression, and the incremental explanatory value
of isoform exposure beyond seed architecture and local AU content.

## Reproducible preparation and first-pass evaluation

After `download_processed.sh`, run `download_hg19_reference.sh` (approximately
0.8 GB for the hg19 2bit reference). The following scripts implement the
benchmark stages:

- `prepare_3p_usage.py`: conservative per-cell-line RefSeq PAS mapping;
- `prepare_joint_3p_usage.py`: common PAS atlas and four-cell-line usage matrix;
- `build_apamirank_inputs.py`: exact strand-correct 2bit sequence extraction;
- `run_architecture_benchmark.py`: per-cell-line miR-124/miR-155 runs;
- `run_joint_architecture.py`: common-atlas runs;
- `evaluate_architecture.py`: matched perturbation evaluation and bootstrap CI;
- `evaluate_joint_architecture.py`: cognate/non-cognate AIR comparison.
- `evaluate_covariate_models.py`: prespecified covariate-adjusted and repeated
  cross-validated comparisons of unweighted and exposure-weighted scores.

The first development checkpoint and its limitations are recorded in
`PRELIMINARY_RESULTS.md`; the adjusted-model checkpoint is recorded in
`COVARIATE_RESULTS.md`.

Run the adjusted analysis after `evaluate_architecture.py` with:

```bash
python benchmarks/gse52531/evaluate_covariate_models.py \
  --evaluation-dir benchmarks/gse52531/evaluation \
  --runs-root benchmarks/gse52531/runs \
  --derived-dir benchmarks/gse52531/derived \
  --output-dir benchmarks/gse52531/covariate_evaluation
```

## Required harmonization before analysis

- Keep the original hg19 coordinates for the first benchmark implementation or
  perform a documented liftOver of every coordinate-bearing input.
- Map historical miRNA names to an explicitly versioned mature-sequence source.
- Derive isoform usage only after excluding internal-priming candidates using
  documented criteria compatible with the original 3P-seq processing.
- Treat the single 3P-seq library per cell type as a limitation; it does not
  provide biological-replicate uncertainty for APA usage.
- Do not use absence from a CLIP or curated database as a confirmed negative.
