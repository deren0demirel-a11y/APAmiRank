# APAmiRank

**APAmiRank** is an alternative polyadenylation-aware miRNA target-ranking
workflow. It enumerates canonical miRNA target sites across ordered 3′UTR
isoforms, records whether each site is shared by proximal isoforms or becomes
available only in intermediate/distal extensions, and produces site- and
gene-level descriptive rankings. When condition-specific isoform fractions are
provided, it also quantifies the fraction of expressed isoforms that expose
each site and reports exposure-weighted condition contrasts.

The repository is intentionally generic: it contains no miR-548ah-3p, VIM,
cell-line, treatment, or predefined rescue assumptions. Any mature miRNA
sequence and compatible APA isoform annotation can be supplied.

## Evidence layers

APAmiRank keeps conceptually distinct evidence layers separate:

1. **Primary binding ranking** from canonical seed architecture, local AU
   context, and optionally RNAhybrid/RNAup evidence.
2. **APA availability annotation** describing common versus extension-specific
   target-site exposure; APA status is not silently treated as proof of binding.
3. **External site annotations** joined after ranking.
4. **Post hoc gene annotations**, including differential-expression or rescue
   results, joined after ranking without changing the primary ranks.
5. **Condition-specific exposure**, calculated from supplied isoform usage
   fractions after the binding-evidence score has been fixed.

The outputs prioritize candidates. They are not calibrated probabilities of
direct targeting and do not replace mutant reporters, AGO occupancy assays, or
endogenous functional perturbation.

The additive condition-level gene score may exceed one. It summarizes the
exposure-weighted evidence of all enumerated sites and is not a probability.

For sites shared by multiple isoforms, local sequence context is calculated
from the longest isoform carrying the site so that the downstream AU window is
not truncated at a proximal PAS.

## Installation

### Lightweight architecture-only mode

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test,workflow,benchmark]"
pytest -q
```

### Full thermodynamic mode

```bash
conda env create -f environment.yml
conda activate apamirank
pytest -q
```

The full environment installs RNAhybrid 2.1.2 and ViennaRNA/RNAup 2.7.2 from
Bioconda/conda-forge. Pinning these versions is recommended for reproducible
comparisons.

## Quick start

Run the included synthetic example:

```bash
apamirank run --config config/config.example.yaml
```

Or through Snakemake:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.example.yaml \
  --cores 4
```

The example uses an artificial miRNA rather than miR-548ah-3p and includes
plus- and minus-strand loci, demonstrating that seed generation and coordinate
mapping are parameterized.

## Preparing QAPA inputs

APAmiRank's stable main interface is a metadata TSV plus strand-correct FASTA.
For QAPA BED7 annotations, these can be generated directly:

```bash
apamirank prepare-qapa \
  --qapa-bed qapa_3utrs.bed \
  --genome-fasta GRCh38.fa \
  --metadata-out input/isoforms.tsv \
  --fasta-out input/isoforms.fa
```

QAPA and genome resources are not redistributed. Record their release,
reference build, and checksum in the analysis provenance.

## Configuration

Copy the example before editing:

```bash
cp config/config.example.yaml config/config.local.yaml
```

Set the mature miRNA identifier and 5′→3′ sequence, input paths, output
directory, scoring mode, and deterministic random seed. Paths are resolved
relative to the configuration file. See [`docs/INPUT_SCHEMA.md`](docs/INPUT_SCHEMA.md).

To enable condition-specific APA analysis, supply `input.isoform_usage` as a
long-form TSV with `isoform_id`, `condition`, and `usage_fraction`. Fractions
must be complete and sum to approximately one within each analyzed
gene/locus/condition. Missing isoforms are rejected rather than silently
treated as zero.

### Architecture-only scoring

`scoring.mode: architecture` ranks sites using canonical seed strength and AU
context. It is fast and useful for screening or CI, but is not equivalent to
the seven-feature thermodynamic analysis.

### Thermodynamic scoring

`scoring.mode: thermodynamic` runs RNAhybrid and RNAup when
`input.thermodynamic_evidence` is null. Alternatively, a complete precomputed
evidence table may be supplied. The seven-feature model uses seed strength, AU
context, RNAhybrid concordance/MFE, and RNAup concordance/interaction/opening
energies. Parsing is strict by default.

## Outputs

| File | Purpose |
|---|---|
| `01_canonical_sites.tsv` | Collapsed canonical sites with APA availability |
| `01b_thermodynamic_evidence.tsv` | Optional RNAhybrid/RNAup evidence |
| `02_site_ranking.tsv` | Primary site-level ranking and robustness |
| `02_leave_one_feature_out.tsv` | Feature-sensitivity summary |
| `03_gene_ranking.tsv` | Best-site gene summary and opportunity adjustment |
| `03b_condition_site_ranking.tsv` | Site exposure and exposure-weighted score in each condition |
| `03c_condition_gene_ranking.tsv` | Additive exposure-weighted gene ranking in each condition |
| `03d_condition_gene_contrasts.tsv` | Pairwise changes in gene evidence score and rank |
| `04_site_ranking_annotated.tsv` | Site ranking plus optional external/post hoc annotations |
| `04_gene_ranking_annotated.tsv` | Gene ranking plus optional post hoc annotations |
| `run_manifest.json` | Output paths and SHA-256 checksums |

The empirical opportunity statistic is descriptive and should not be reported
as a target probability or inferential P value.

## Tests

```bash
pytest -q
```

Tests verify generic seed construction, plus/minus-strand site mapping,
common/distal availability, strict isoform-usage validation, condition-specific
site exposure, post hoc rank preservation, and a complete architecture-mode
run.

The checks performed for this packaged release candidate are recorded in
[`docs/VERIFICATION.md`](docs/VERIFICATION.md).

The prespecified analysis logic for the public-data validation is in
[`docs/BENCHMARK_PROTOCOL.md`](docs/BENCHMARK_PROTOCOL.md). It should be frozen
before outcome labels are joined to rankings.

The orthogonal positive--unlabeled CLASH/CLEAR-CLIP benchmark, its exact source
files, commands, metrics, and interpretation limits are documented in
[`benchmarks/direct_interactions/README.md`](benchmarks/direct_interactions/README.md).

The low-throughput validation layer based on the dated official-miRTarBase
Reporter Assay export is documented in
[`benchmarks/reporter_validation/README.md`](benchmarks/reporter_validation/README.md).

The manual construct-level audit requiring a WT/mutant 3′UTR reporter is
documented in
[`benchmarks/strict_luciferase/README.md`](benchmarks/strict_luciferase/README.md).
Its pilot summaries are descriptive only because the small candidate set was
reviewed after the broader reporter ranks had been inspected.

The next-stage, rank-independent pair-level holdout is locked in
[`benchmarks/locked_luciferase_holdout/`](benchmarks/locked_luciferase_holdout/README.md).
Its 33-pair manifest was selected deterministically without score columns; a
fail-closed gate prevents evaluation until every original-paper screening
decision is complete and the final form is checksummed. Version 0.2.18 completes
all 33 original-paper reviews: 21 included and 12 excluded. The final screening
form was frozen by SHA-256 before the prespecified rank file was joined. Of 21
strictly included pairs, 12 were observable in the ranking universe; their mean
and median APA-exposure percentiles were 0.7847 and 0.8185, respectively. Four
of 12 were in the top decile and seven of 12 in the top quintile. The complete
pair-level results and coverage limitation are documented in
[`benchmarks/locked_luciferase_holdout/RESULTS.md`](benchmarks/locked_luciferase_holdout/RESULTS.md).
An explicitly post-unblinding secondary analysis compares the observable pairs
with a miRNA-matched rank-uniform null and adds miRNA-cluster bootstrap
intervals. Conditional ranking enrichment remains detectable after Holm
adjustment, whereas the cluster-bootstrap interval for the incremental
APA-minus-unweighted effect crosses zero. These exploratory results do not
resolve the 57.1% coverage limitation.

A manuscript-ready synthesis of the validation methods, results, figure
legend, limitations, and permitted claim boundary is provided in
[`manuscript/VALIDATION_METHODS_RESULTS.md`](manuscript/VALIDATION_METHODS_RESULTS.md),
with a styled Word edition in the same directory.

## Citation and release status

This is version 0.2.18, a research-validation release candidate. Complete the items in
[`docs/PUBLIC_RELEASE_CHECKLIST.md`](docs/PUBLIC_RELEASE_CHECKLIST.md) before a
public release. In particular, an open-source license has not been selected;
no license is granted by this draft package.
