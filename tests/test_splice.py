from transcript_isoform_plotter.splice import (
    matching_annotations,
    read_gene_table,
    read_splice_annotations,
)


def test_splice_parser_and_tolerant_matching(tmp_path):
    path = tmp_path / "splice.txt"
    path.write_text(
        "# affected_exon\tsplicing_boundary\n"
        "98, 202\t200, 400\n"
    )
    annotations = read_splice_annotations(path)
    assert matching_annotations((100, 200), annotations, tolerance=2) == annotations
    assert matching_annotations((100, 200), annotations, tolerance=1) == []


def test_gene_table_supports_optional_columns_and_repeated_genes(tmp_path):
    path = tmp_path / "genes.tsv"
    path.write_text(
        "gene\taffected_exon\tsplicing_boundary\n"
        "GENE1\t100, 200\t200, 400\n"
        "GENE1\t400, 500\n"
        "GENE2\n"
    )
    genes, annotations = read_gene_table(path)
    assert genes == ["GENE1", "GENE2"]
    assert annotations["GENE1"] == [
        {"exon": (100, 200), "boundary": (200, 400)},
        {"exon": (400, 500), "boundary": None},
    ]
    assert annotations["GENE2"] == []
