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
    # The affected interval still matches through its >80% exon overlap.
    assert matching_annotations((100, 200), annotations, tolerance=1) == annotations


def test_partial_affected_intervals_match_gtf_exons():
    cases = [
        # SH3YL1: affected interval overhangs the exon start by 5 bp.
        ((224868, 224920), (224863, 224920)),
        # PLEKHA1: affected interval is a subregion of a larger exon.
        ((122428276, 122428420), (122428275, 122428316)),
        # PAK6: affected interval is internal to a much larger exon.
        ((40275927, 40277487), (40276105, 40277115)),
    ]
    for exon, affected_interval in cases:
        annotation = [{"exon": affected_interval, "boundary": None}]
        assert matching_annotations(exon, annotation, tolerance=2) == annotation

    mostly_nonoverlapping = [{"exon": (150, 300), "boundary": None}]
    assert matching_annotations((100, 175), mostly_nonoverlapping, tolerance=2) == []


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
