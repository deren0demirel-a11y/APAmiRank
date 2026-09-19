# Locked-holdout screening progress

Checkpoint date: 2026-09-10  
Protocol: `APAmiRank-locked-luciferase-holdout-v1`

## Current accounting

| State | Pairs |
|---|---:|
| Locked | 33 |
| Complete — include | 21 |
| Complete — exclude | 12 |
| Pending | 0 |

All decisions were completed and the screening form was frozen at SHA-256
`8f7ecc1f89a3819365fdfc5ccbf7799d6d13ac7042874c0023d93fe779160ae5`
before any APAmiRank rank, score, percentile, or outcome table was joined. The
evaluation gate was then opened with that exact checksum.

## Discovery checkpoint

The fielded PubMed discovery script returned 298 pair-article candidates and
299 unique PMIDs. Five pairs had no title/abstract co-mention candidate:
miR-186-5p–CSNK2A1, miR-186-5p–MOB1A, miR-25-3p–TCEAL1,
miR-92a-3p–CPEB2, and miR-92a-3p–MTO1. This is a retrieval limitation, not an
exclusion decision; abstracts can omit either interaction partner.

The official miRTarBase record interface returned HTTP 429 on the detail pages
after one controlled retry. The dated human Reporter Assay CSV remained
available and its SHA-256 matched the frozen benchmark snapshot:
`5892eddbf8c51d35bf8379254d87c5c084aa34a9d8c4527cb1e16d3d7cd21c08`.

## Completed original-paper reviews

| Pair | PMID | Decision | Construct-level basis |
|---|---:|---|---|
| hsa-miR-124-3p–GRB2 | 28496318 | Include | Human GRB2 3′UTR in pISO; cognate site mutant; WT repression by miR-124 lost in mutant; qPCR and Western blot support in PASMCs. |
| hsa-miR-124-3p–PLEC | 35024318 | Include | PLEC 0–1024-nt 3′UTR in pEZX-MT06; three site-mutant constructs; >40% WT repression; response lost when distal MRE2 was mutated; endogenous PLEC protein support. |
| hsa-miR-124-3p–MAPK14 | 23109423 | Include | Full-length human p38α/MAPK14 3′UTR in pCI-FLuc; cognate miR-124 seed mutant; WT repression by miR-124 lost in the mutant; endogenous p38α protein and neuronal antagomir support. |
| hsa-miR-124-3p–AHR | 26802080 | Include | Human AHR 3′UTR reporter in Caco-2 and HT-29 cells; seed-matching nucleotides mutated; pre-miR-124 repression and anti-miR-124 activation were both largely abolished by mutation; endogenous AHR protein and patient-tissue support. |
| hsa-miR-124-3p–ROCK2 | 21672940 | Exclude — `FULL_TEXT_UNAVAILABLE` | The abstract reports direct ROCK2 3′UTR targeting and endogenous mRNA/protein suppression, but the WT/cognate-mutant design could not be verified. A 2021 correction republishes corrected figures, including Figure 6, so the record is excluded fail-closed rather than inferred from the abstract. |
| hsa-miR-186-5p–AKAP12 | 20979053 | Exclude — `FULL_TEXT_UNAVAILABLE` | The original Hepatology study is cited for an AKAP12 3′UTR luciferase assay, but its construct methods and figures were inaccessible; the open 2018 follow-up confirms exact-miR-186-5p endogenous regulation but does not repeat the reporter assay. |
| hsa-miR-186-5p–CSNK2A1 | 23137536 | Exclude — `FULL_TEXT_UNAVAILABLE` | The abstract reports CKIIα/CSNK2A1 3′UTR reporter analysis after a joint four-miRNA perturbation, but an exact miR-186-only WT response and cognate-mutant rescue could not be verified from the inaccessible full text. |
| hsa-miR-186-5p–MOB1A | 33413543 | Exclude — `NOT_3UTR` | The original methods explicitly clone the 5′ flanking sequence of the MOB1 promoter into pGL3-basic. WT/mutant reporter response, RNA pull-down, and qPCR support are reported, but the construct is not a target-gene 3′UTR. |
| hsa-miR-155-5p–CKAP5 | 21799781 | Exclude — `NO_MUTANT` | Full-length CKAP5 3′UTR in pMIR-REPORT was repressed by miR-155 in HEK293T cells, with SILAC and Western-blot support, but no cognate site-mutant or site-deletion reporter was tested. |
| hsa-miR-155-5p–SEL1L | 23661430 | Exclude — `NO_MUTANT` | A SEL1L 3′UTR segment containing the predicted site and >100-bp flanks was tested in PL18 and MDA-Panc3 cells; endogenous and precursor-overexpression evidence supported regulation, but no cognate mutant reporter was tested. |
| hsa-miR-155-5p–APC | 22610915 | Include | Approximately 300-bp APC 3′UTR insert in pGL3; four-nucleotide mutation of the miR-155 site; WT repression in miR-155-overexpressing HEK293 cells was completely abrogated by mutation; endogenous APC protein and patient-sample support. |
| hsa-miR-155-5p–HBP1 | 33769664 | Include | Human HBP1 3′UTR in pYr-MirTarget; cognate site mutant; miR-155 inhibited WT activity but had no inhibitory effect on the mutant; endogenous HBP1 and human-tissue support. Independently replicated with a conserved site-B mutant in PMID 24675724. |
| hsa-miR-155-5p–CLDN1 | 36819593 | Include | CLDN1 3′UTR in pGL3; four-nucleotide cognate-site mutant; miR-155-5p markedly repressed WT in NCM460 cells and mutation completely abolished suppression; endogenous CLDN1 and human-tissue support. |
| hsa-miR-16-5p–CDS2 | 23941513 | Include | Entire human CDS2 3′UTR in p-MIR-report; miR-16 repressed WT and cognate seed-complementary mutations almost fully rescued activity in A549 cells; multi-cell-line transcriptomic support. |
| hsa-miR-16-5p–WEE1 | 33583218 | Include | Human WEE1 3′UTR in pGL3-basic; cognate seed-site mutant; miR-16-5p repressed WT in H838 and A549 cells but did not significantly affect the mutant; endogenous WEE1 mRNA/protein and human-tissue support. |
| hsa-miR-16-5p–YAP1 | 34632051 | Include | A 500-bp human YAP1 3′UTR fragment in psiCHECK-2; cognate seed-site mutant; miR-16-5p overexpression repressed WT but not mutant in HuCCT1 cells, while inhibitor reciprocally increased only WT activity in QBC939 cells; endogenous protein, Ago2-RIP, and rescue support. |
| hsa-miR-16-5p–KRAS | 27857191 | Include | Entire human KRAS 3′UTR in pMIR-REPORT; cognate seed-site mutant; pre-miR-16 repressed only WT in SW480 cells and anti-miR-16 reciprocally increased only WT in Caco2 cells; endogenous protein and phenotypic-rescue support. |
| hsa-miR-16-5p–PPM1D | 20668064 | Include | A 1.1-kbp human WIP1/PPM1D 3′UTR in pRL; cognate six-nucleotide site deletion; pre-miR-16 repressed WT and antagomir-16 increased WT, while deletion almost completely abolished both effects; endogenous protein and rescue support. |
| hsa-miR-17-5p–JAK1 | 38233760 | Include | WT and cognate-site-mutant human JAK1 3′UTR reporters in HEK293T cells; miR-17-5p repressed WT but not mutant; Ago2 RIP and reciprocal endogenous JAK1 mRNA/protein responses in THP-1 and HL-60 cells. |
| hsa-miR-17-5p–RUNX1 | 17589498 | Include | Human AML1/RUNX1 target 3′UTR and 60-bp AML1short WT/mutant reporter; mutant relieved repression and individual anti-miR-17-5p selectively increased WT/MUT activity; endogenous AML1 protein and monocytopoiesis support. |
| hsa-miR-17-5p–CDKN1A | 26482648 | Include | A 1531-bp human CDKN1A 3′UTR in psiCHECK-2 with both miR-17-family sites mutated; miR-17 repressed WT but not mutant; endogenous CDKN1A and independent mRNA-protector rescue support. |
| hsa-miR-17-5p–ITGB8 | 24920276 | Exclude — `OTHER_WITH_NOTE` | The linked article tested WT/mutant ADAR1 reporters, not ITGB8. ITGB8 received qPCR-only follow-up in C81-61 cells, and the authors described its regulation as questionable. |
| hsa-miR-17-5p–STAT3 | 23059786 | Include | Full-length human STAT3 3′UTR in pMIR-REPORT; miR-17 repressed WT, while mutation of the second site and especially both conserved sites abrogated repression; endogenous STAT3 protein support in U937T cells. |
| hsa-miR-25-3p–FBXW7 | 25512615 | Include | WT FBXW7 3′UTR reporter was reduced by approximately 50% after miR-25 expression, while the cognate-site-mutant reporter showed no obvious response; endogenous FBXW7 and AURKA-pathway support. |
| hsa-miR-25-3p–MDM2 | 22431589 | Include | WT and seed-region-mutant MDM2 3′UTR reporters in U87 cells; miR-25 repressed WT but not the responsive-region mutant; endogenous MDM2, p53 and p21 support. |
| hsa-miR-25-3p–CCL26 | 22431589 | Exclude — `WRONG_TARGET` | The linked miRTarBase description and source paper test TSC1, not CCL26; no CCL26 3′UTR reporter is reported in that source. |
| hsa-miR-25-3p–TP53 | 20935678 | Include | Human TP53 3′UTR in pRL-TK; positions 92–98 changed from `TGCAATA` to `ACGTTAT`; miR-25 repression was abolished by the cognate mutation; endogenous p53 support. |
| hsa-miR-25-3p–TCEAL1 | 23028803 | Exclude — `WRONG_TARGET` | The linked miRTarBase description and source paper test BIM as the miR-25 target, not TCEAL1; no TCEAL1 3′UTR reporter is reported in that source. |
| hsa-miR-92a-3p–NRF1 | 33505436 | Exclude — `OTHER_WITH_NOTE` | Exact-miR-92a-3p WT/mutant Nrf1 reporter evidence is reported, but the biological model is mouse and the insert species is not identified; a human NRF1 target record therefore cannot be admitted under the locked source-population rule. |
| hsa-miR-92a-3p–CPEB2 | 20660482 | Include | Human CPEB2 3′UTR fragment in pGL3-Control; Mut92 cognate-site mutant; miR-92 overexpression repressed WT but not mutant, while LNA-mediated depletion reciprocally raised WT to mutant levels; endogenous CPEB2 mRNA support. |
| hsa-miR-92a-3p–CDH1 | 27801803 | Include | WT and cognate-site-mutant CDH1 3′UTR reporters in human glioma cells; miR-92a-3p inhibition de-repressed WT but not mutant; endogenous CDH1 mRNA/protein support. |
| hsa-miR-92a-3p–STAT3 | 23820254 | Exclude — `WRONG_TARGET` | The linked paper places STAT3 upstream of miR-92a and tests a RECK 3′UTR reporter; it reports no STAT3 3′UTR reporter. |
| hsa-miR-92a-3p–MTO1 | — | Exclude — `FULL_TEXT_UNAVAILABLE` | The frozen database row labels reporter and orthogonal support, but the original reporter article and construct could not be resolved from accessible sources; no evidence was inferred. |

## Locked evaluation result

The gate joined the 21 strict inclusions to the prespecified consensus-rank
table. Twelve were observable (57.1% candidate coverage). The observable set
had mean and median APA-exposure percentiles of 0.7847 and 0.8185; four of 12
were in the top decile and seven of 12 in the top quintile. The median paired
APA-exposure minus unweighted-additive percentile was 0.0651. Missing candidates
were retained as missing and were not imputed. Detailed results are in
`RESULTS.md` and `evaluation/`.

A secondary rank-uniform null and cluster-bootstrap analysis was added only
after this descriptive result was available. It is therefore labelled
post-unblinding and exploratory in `SECONDARY_ANALYSIS_PROTOCOL.md`; it does not
change the screening checksum, candidate accounting, or primary evaluation.

The PLEC article has an author correction (PMID 37655312; DOI
10.1016/j.apsb.2023.07.002). It replaces a mislabeled representative image in
Fig. 4A and states that the results and conclusions are unaffected. The strict
reporter evidence is in Fig. 3, so the correction does not alter the inclusion
decision.

The ROCK2 exclusion is an evidence-access decision, not evidence against the
biological interaction. PMID 21672940 reports direct miR-124 targeting of the
ROCK2 3′UTR and suppression of endogenous ROCK2 mRNA and protein. However, its
2021 correction (PMID 34497120) identifies incorrect representative images and
republishes corrected Figure 6 and Figure S5. The corrected figure and original
construct methods could not be inspected reliably, so the required WT versus
cognate-mutant reporter response was not assumed from the abstract.

The AKAP12 exclusion is also fail-closed. PMID 20979053 is the original source
later cited for the AKAP12 3′UTR luciferase experiment, but the publisher's
construct methods and figures could not be inspected reliably. The open
follow-up PMID 29653561 verifies reciprocal endogenous AKAP12 mRNA/protein
responses after exact miR-186-5p inhibition or overexpression, but it does not
repeat a reporter assay. PMID 34223793 is a retracted circRNA/CCND2 article in
which AKAP12 appears only in background text and is not used as evidence.

The CSNK2A1 exclusion uses the original PMID 23137536. Its abstract reports
reporter analysis and RT-PCR consistent with targeting of the CKIIα 3′UTR, but
describes miR-186, miR-216b, miR-337-3p, and miR-760 as a joint perturbation.
Without accessible construct-level text, neither a miR-186-only WT response nor
cognate-mutant rescue was inferred.

The MOB1A exclusion is based on the open original article PMID 33413543. The
paper reports miR-186 mimic/inhibitor reporter effects, a mutant lacking that
response, RNA pull-down, and reciprocal endogenous qPCR evidence. Nevertheless,
its methods state that the cloned insert was the 5′ flanking sequence of the
MOB1 promoter. It therefore fails the prespecified 3′UTR criterion under
`NOT_3UTR` rather than on biological plausibility.

The CKAP5 article supports a direct miR-155 interaction biologically, including
a 0.65-fold full-length 3′UTR reporter signal and protein-level confirmation.
It is excluded only from this stricter mutant-rescue holdout because the paper
does not report a cognate site-mutant or site-deletion construct.

The SEL1L article likewise supports miR-155-dependent regulation, but its
reporter comparison used high- versus low-endogenous-miR-155 PDA cell contexts
and an empty reporter control, without a cognate binding-site mutant. It is
therefore excluded from the strict holdout under `NO_MUTANT`.

The ITGB8 exclusion records a source-evidence mismatch rather than negative
reporter biology. PMID 24920276 used WT/mutant psiCheck2 constructs only for
ADAR1. Its ITGB8 experiment was qPCR-only, and the authors characterized the
observed regulation as questionable. Because no ITGB8 3′UTR reporter was tested,
the pair is coded `OTHER_WITH_NOTE` rather than `NO_MUTANT`.

The STAT3 inclusion is anchored to PMID 23059786. The full-length human STAT3
3′UTR contained two reported miR-17/20a sites. Exact miR-17 mimic repressed the
WT reporter in HEK293T cells; the M2 single-site mutant and especially the
M1+M2 double mutant clearly abrogated repression. Ectopic miR-17 also reduced
endogenous STAT3 protein in U937T human AML cells.

The APC inclusion is anchored to the open full text of PMID 22610915 rather
than the paywalled PTC article PMID 23796566. The Hepatology study directly
reports WT repression and complete loss of suppression after a four-nucleotide
mutation in the cognate APC 3′UTR site, satisfying every locked criterion.

The HBP1 inclusion is anchored to the human-3′UTR experiment in PMID 33769664.
That study used a pYr-MirTarget WT construct and a QuikChange cognate-site
mutant in HEK293T cells; WT activity was inhibited by miR-155, while the mutant
showed no inhibitory response. The construct-level conclusion is independently
supported by PMID 24675724, where mutation of conserved HBP1 site B almost
rescued reporter suppression. The unrelated 2016 HBP1/MMR article retracted in
2024 is not used as evidence.

The CLDN1 inclusion uses the explicit mature-arm experiment in PMID 36819593.
The original article reports a pGL3 CLDN1 3′UTR WT/mutant comparison in human
NCM460 cells: 40 nM miR-155-5p repressed WT after 48 hours, whereas a
four-nucleotide mutation in the cognate site completely eliminated repression.

The CDS2 inclusion is based on PMID 23941513. The entire human CDS2 3′UTR was
placed downstream of firefly luciferase, and bases complementary to positions
2–8 of mature miR-16 were mutated. Pre-miR-16 repressed WT, anti-miR-16 raised
reporter activity, and the cognate mutation almost fully rescued repression in
three independent A549 experiments. The paper does not report CDS2-specific
protein validation; that limitation is retained explicitly in the form.

The WEE1 inclusion is anchored to the explicit mature-arm WT/mutant experiment
in PMID 33583218. In both H838 and A549 cells, miR-16-5p reduced the pGL3-basic
WT reporter to approximately 0.4 of control, while the cognate seed-site mutant
remained near control and was not significantly affected. The original figure,
methods, and caption were inspected. PMID 24317448 independently supports WT
WEE1 3′UTR regulation but is not used as the mutant-rescue evidence.

The YAP1 inclusion is anchored to PMID 34632051 rather than the alternate
LINC00649 study. The Molecular Therapy – Oncolytics article explicitly cloned a
500-bp YAP1 3′UTR fragment into psiCHECK-2 and changed the cognate WT
`UGCUGCUG` site to `GGATCCG`. Mature miR-16-5p overexpression repressed WT but
not mutant activity in HuCCT1 cells; miR-16-5p inhibition produced the
reciprocal WT-only increase in QBC939 cells. The original Figure 5, methods,
PubMed record, and publication-notice status were inspected.

The KRAS inclusion is based on PMID 27857191. The entire human KRAS 3′UTR was
cloned into pMIR-REPORT, and the cognate `TGCTGCT` seed-recognizing sequence was
changed to `ACGACGA`. Pre-miR-16 repressed only WT activity in SW480 cells;
anti-miR-16 reciprocally increased only WT activity in Caco2 cells. A 2021
author correction (DOI `10.1038/s41598-021-99119-w`) replaces partially
duplicated cell-invasion images in Figure 4D. It does not alter the Figure 3E
luciferase experiment used for this holdout decision.

The PPM1D inclusion uses PMID 20668064. A 1.1-kbp human WIP1/PPM1D 3′UTR was
cloned into a pRL Renilla reporter, and six nucleotides recognized by miR-16
were deleted. Pre-miR-16 reduced WT activity by approximately 75%, whereas
antagomir-16 increased it by approximately 50%; the cognate deletion almost
completely abolished both effects. Endogenous Wip1 protein regulation in U2OS
and MCF-7 cells and a Wip1 phenotypic-rescue experiment provide orthogonal
support.

The JAK1 inclusion is anchored to PMID 38233760. The methods explicitly name
JAK1-3UTR (miR-17-5p)-WT and -MUT constructs, and the original Figure 7 shows
the cognate-site alignment. Mature miR-17-5p repressed WT but not MUT in
HEK293T cells; Ago2 RIP and reciprocal JAK1 mRNA/protein responses in two human
AML lines provide orthogonal support.

The RUNX1 inclusion is anchored to the original 2007 study (PMID 17589498),
not the later 2015 mechanistic follow-up. Figure 3 tests a 60-bp human
AML1short 3′UTR site and a four-base cognate mutant. Because miR-17-5p, miR-20a,
and miR-106a share the tested seed site, exact miR-17-5p attribution is retained
only because the paper also tests anti-miR-17-5p individually and shows the
WT/MUT-selective reporter response.

The CDKN1A inclusion uses PMID 26482648 rather than the 2014 synovial-sarcoma
paper, which reports only a WT p21 3′UTR reporter. The selected study tests the
1531-bp human CDKN1A 3′UTR with both miR-17-family sites mutated, shows
miR-17-dependent repression only for WT, and independently blocks the same
interaction with a CDKN1A mRNA protector.

The three miR-25-3p inclusions satisfy the locked WT/cognate-mutant rule in
their original articles. PMID 25512615 reports approximately 50% repression of
the WT FBXW7 3′UTR reporter with no obvious response from its cognate mutant.
PMID 22431589 reports miR-25-dependent repression of WT MDM2 3′UTR and loss of
that response after mutation of the responsive seed region. PMID 20935678
specifies the TP53-site change from `TGCAATA` to `ACGTTAT` and reports that it
abolished miR-25 regulation while preserving regulation by another miRNA.

The CCL26 and TCEAL1 exclusions are provenance failures, not negative
biological results. The CCL26-labelled miRTarBase record links to PMID 22431589,
whose original description and article test TSC1. The TCEAL1-labelled record
links to PMID 23028803, whose original description and article test BIM for
miR-25 and p21 for miR-106b. Neither linked source reports a reporter construct
for the gene named by the locked pair, so both are coded `WRONG_TARGET`.

The final miR-92a-3p decisions preserve the same distinction. NRF1 is excluded
because a human reporter insert cannot be established from a mouse-model paper,
not because its WT/mutant response failed. STAT3 is a target-role mismatch: the
linked paper's direct target and reporter insert are RECK. MTO1 is an
evidence-access exclusion because the database-labelled construct could not be
verified. CPEB2 and CDH1 independently satisfy the exact mature-miRNA, human
3′UTR, WT-repression and cognate-mutant-rescue criteria.
