# Methods template

Canonical target-site architecture was evaluated using APAmiRank v0.2.2
(Alternative Polyadenylation-aware miRNA target ranking). Strand-correct 3′UTR
isoforms were ordered from proximal to distal polyadenylation-site usage within
each locus. Canonical 8mer, 7mer-m8, 7mer-A1, and 6mer target motifs were
generated from the supplied mature miRNA sequence and enumerated across all
eligible isoforms. Genomically identical sites present in multiple isoforms
were collapsed while retaining their first and last available PAS ranks and
their classification as common or extension-specific sites.

[Select the applicable scoring sentence.]

**Architecture-only mode:** Sites were ranked using an equal-weight descriptive
score comprising canonical seed strength and the percentile of local AU
content. This mode does not incorporate thermodynamic predictions.

**Thermodynamic mode:** Sites were ranked using an equal-weight seven-feature
descriptive score comprising canonical seed strength, local AU context,
RNAhybrid coordinate concordance and minimum free energy, and RNAup coordinate
concordance, total interaction energy, and target-opening energy. Rank
robustness was assessed by leave-one-feature-out analysis and random positive
feature-weight resampling. A Monte Carlo opportunity adjustment summarized the
best site per gene relative to genes with the same number of enumerated sites.

External prediction scores and expression/rescue information were added only
after the primary binding-evidence ranks were fixed. APAmiRank outputs are
candidate-prioritization catalogues and do not establish direct endogenous
miRNA targeting or functional mediation.

**Condition-specific APA exposure:** For each site and condition, site exposure
was calculated as the sum of the supplied usage fractions of all isoforms
carrying that genomic site. The fixed binding-evidence score was multiplied by
site exposure to obtain an exposure-weighted descriptive score. Gene-level
condition scores were calculated as the sum of exposure-weighted evidence
scores across enumerated sites, and pairwise condition contrasts were reported.
Input fractions were required to be complete and to sum to 1 within each
locus/condition (tolerance: [REPORT VALUE]). These additive scores were treated
as relative prioritization statistics rather than calibrated probabilities of
targeting.
