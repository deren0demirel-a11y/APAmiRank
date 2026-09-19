# Secondary uncertainty-analysis specification

Specification date: 2026-09-10  
Status: post-unblinding, exploratory secondary analysis

This analysis was designed after the initial locked-holdout descriptive ranks
had been calculated. It was not part of the 2026-09-08 screening protocol and
must not be described as preregistered or as a prespecified primary inferential
test. The original candidate selection, literature decisions, screening-form
checksum and rank join remain unchanged.

## Rank-uniform matched null

For each of the six miRNAs with at least one observable strict inclusion, the
simulation preserves the observed number of holdout pairs and the size of that
miRNA's full `at_least_2_of_4_contexts` ranking universe. In each of 50,000
iterations, ordinal rank positions are sampled without replacement within each
miRNA and converted to percentiles using the same rank-to-percentile formula.
The pooled mean, median, top-decile fraction and top-quintile fraction are then
recalculated. One-sided Monte Carlo P values use the plus-one correction, and
the four correlated tests are adjusted by Holm's method.

This is a theoretical rank-uniform null, not an empirical negative-gene set. It
does not preserve score-tie multiplicities because the complete per-gene rank
tables were not retained in the packaged benchmark. It is conditional on the
12 observable pairs and therefore does not correct potential bias from the nine
strict inclusions absent from the ranking universe.

## Cluster-aware bootstrap

Each miRNA is treated as a cluster. Cluster-specific statistics are calculated
first, after which the six complete clusters are sampled with replacement for
50,000 iterations and averaged with equal miRNA weight. Percentile 2.5% and
97.5% quantiles are reported. For the paired APA-minus-unweighted comparison,
the cluster statistic is the within-miRNA median difference.

With only six clusters, these intervals are unstable estimates of uncertainty
and should be interpreted as a sensitivity analysis rather than definitive
population-level inference. No P value is assigned to the paired APA increment.

Random seed: `20260910` for the null and `20260911` for the cluster bootstrap.
