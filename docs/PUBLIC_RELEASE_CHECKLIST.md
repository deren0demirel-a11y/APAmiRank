# Public-release checklist

- [ ] Confirm repository authorship and contribution statements.
- [ ] Obtain institutional/co-author approval for an open-source license.
- [ ] Replace `LICENSE_PENDING.md` with the approved license.
- [ ] Replace the placeholder GitHub URL in `CITATION.cff`.
- [ ] Add the manuscript DOI or preprint identifier when available.
- [ ] Confirm that no patient, unpublished collaborator, or licensed reference
      data are committed.
- [ ] Keep QAPA, reference-genome, TargetNet, RNAhybrid, and ViennaRNA resources
      external unless their licenses explicitly allow redistribution.
- [ ] Run `pytest -q` and the synthetic example.
- [ ] Verify the release archive checksum.
- [ ] Create a versioned GitHub release and archive it with Zenodo if desired.
