from __future__ import annotations

from .io import read_tsv, write_tsv


def _unique_map(rows, key, label):
    result = {}
    for row in rows:
        value = row.get(key, "")
        if not value:
            raise ValueError(f"{label} is missing key column/value: {key}")
        if value in result:
            raise ValueError(f"Duplicate {label} key: {value}")
        result[value] = row
    return result


def annotate_rankings(site_ranking, gene_ranking, output_sites, output_genes, external_sites=None, posthoc_genes=None):
    sites, genes = read_tsv(site_ranking), read_tsv(gene_ranking)
    external = _unique_map(read_tsv(external_sites), "site_id", "external site annotation") if external_sites else {}
    posthoc = _unique_map(read_tsv(posthoc_genes), "gene", "post hoc gene annotation") if posthoc_genes else {}

    def add_external(row):
        result = dict(row)
        hit = external.get(row["site_id"])
        result["external_site_annotation_present"] = "Yes" if hit else "No"
        if hit:
            for key, value in hit.items():
                if key != "site_id":
                    result[f"external_{key}"] = value
        return result

    def add_posthoc(row):
        result = dict(row)
        hit = posthoc.get(row["gene"])
        result["posthoc_gene_annotation_present"] = "Yes" if hit else "No"
        if hit:
            for key, value in hit.items():
                if key != "gene":
                    result[f"posthoc_{key}"] = value
        return result

    annotated_sites = [add_posthoc(add_external(row)) for row in sites]
    annotated_genes = [add_posthoc(row) for row in genes]
    if [row["equal_weight_rank"] for row in annotated_sites] != [row["equal_weight_rank"] for row in sites]:
        raise AssertionError("Primary site ranks changed during annotation")
    if [row["opportunity_adjusted_descriptive_rank"] for row in annotated_genes] != [row["opportunity_adjusted_descriptive_rank"] for row in genes]:
        raise AssertionError("Primary gene ranks changed during annotation")
    write_tsv(output_sites, annotated_sites)
    write_tsv(output_genes, annotated_genes)
    return annotated_sites, annotated_genes
