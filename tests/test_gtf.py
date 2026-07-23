import gzip

import pytest

from transcript_isoform_plotter.gtf import (
    filter_gene_transcripts,
    normalize_reference_gene_types,
    normalize_transcript_filters,
    parse_region,
    read_genes,
    read_reference_region,
)


GTF = """##format: gtf
chr22\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1"; transcript_name "TX1.1"; transcript_type "protein_coding"; tag "basic"; tag "MANE_Select"; tag "appris_principal_1";
chr22\ttest\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\tCDS\t150\t200\t.\t+\t0\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\texon\t400\t500\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\tCDS\t400\t450\t.\t+\t0\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\ttranscript\t120\t450\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2"; transcript_type "retained_intron"; tag "Ensembl_canonical";
chr22\ttest\texon\t120\t220\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2";
chr22\ttest\texon\t350\t450\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2";
chr1\ttest\ttranscript\t1000\t1300\t.\t+\t.\tgene_id "G2"; gene_name "GENE2"; transcript_id "TX3"; transcript_type "processed_transcript";
"""


@pytest.mark.parametrize("compressed", [False, True])
def test_reads_plain_and_gzipped_gtf(tmp_path, compressed):
    path = tmp_path / ("test.gtf.gz" if compressed else "test.gtf")
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(GTF)
    else:
        path.write_text(GTF)
    genes = read_genes(path, ["GENE1"])
    assert list(genes["GENE1"]["transcripts"]) == ["TX1", "TX2"]
    assert genes["GENE1"]["transcripts"]["TX1"]["label"] == "TX1"
    assert genes["GENE1"]["transcripts"]["TX1"]["exons"] == [(100, 200), (400, 500)]
    assert genes["GENE1"]["transcripts"]["TX1"]["cds"] == [(150, 200), (400, 450)]


def test_matches_gene_id_and_uses_transcript_span_fallback(tmp_path):
    path = tmp_path / "test.gtf"
    path.write_text(GTF)
    genes = read_genes(path, ["G2"])
    assert genes["G2"]["transcripts"]["TX3"]["exons"] == [(1000, 1300)]
    assert genes["G2"]["transcripts"]["TX3"]["fallback_span"] is True


def test_missing_gene_is_reported(tmp_path):
    path = tmp_path / "test.gtf"
    path.write_text(GTF)
    with pytest.raises(ValueError, match="not found"):
        read_genes(path, ["NO_SUCH_GENE"])


def test_transcript_filters_use_tags_types_and_cds(tmp_path):
    path = tmp_path / "test.gtf"
    path.write_text(GTF)
    genes = read_genes(path, ["GENE1"])
    expected = {
        "basic": ["TX1"],
        "mane": ["TX1"],
        "appris_principal": ["TX1"],
        "canonical": ["TX2"],
        "protein_coding": ["TX1"],
        "coding": ["TX1"],
        "noncoding": ["TX2"],
        "type:retained_intron": ["TX2"],
        "basic,canonical": ["TX1", "TX2"],
    }
    for expression, transcript_ids in expected.items():
        filters = normalize_transcript_filters(expression)
        filtered = filter_gene_transcripts(genes, filters)
        assert list(filtered["GENE1"]["transcripts"]) == transcript_ids


def test_transcript_filter_errors_are_clear(tmp_path):
    with pytest.raises(ValueError, match="Unknown transcript filter"):
        normalize_transcript_filters("not_a_filter")

    path = tmp_path / "test.gtf"
    path.write_text(GTF)
    genes = read_genes(path, ["GENE2"])
    with pytest.raises(ValueError, match="No transcripts passed"):
        filter_gene_transcripts(genes, ["mane"])


def test_region_parser_allows_commas():
    assert parse_region("chr22:41,900,000-42,100,000") == (
        "chr22",
        41_900_000,
        42_100_000,
    )


def test_reference_region_prefers_curated_refseq_accession(tmp_path):
    path = tmp_path / "reference.gtf"
    path.write_text(
        'chr22\ttest\ttranscript\t100\t900\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "ENST1"; tag "MANE_Select";\n'
        'chr22\ttest\texon\t100\t300\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "ENST1"; tag "MANE_Select";\n'
        'chr22\ttest\texon\t700\t900\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "ENST1"; tag "MANE_Select";\n'
        'chr22\ttest\ttranscript\t200\t800\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "NM_000001.1";\n'
        'chr22\ttest\texon\t200\t350\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "NM_000001.1";\n'
        'chr22\ttest\texon\t650\t800\t.\t+\t.\t'
        'gene_id "G"; gene_name "GENE"; transcript_id "NM_000001.1";\n'
    )
    reference = read_reference_region(path, ("chr22", 1, 1000))
    assert reference["genes"]["GENE"]["transcript_id"] == "NM_000001.1"


def test_reference_gene_type_filters(tmp_path):
    path = tmp_path / "reference_types.gtf"
    path.write_text(
        'chr22\ttest\ttranscript\t100\t300\t.\t+\t.\t'
        'gene_id "G1"; gene_name "CODING"; transcript_id "TX1"; '
        'gene_type "protein_coding";\n'
        'chr22\ttest\texon\t100\t300\t.\t+\t.\t'
        'gene_id "G1"; gene_name "CODING"; transcript_id "TX1"; '
        'gene_type "protein_coding";\n'
        'chr22\ttest\ttranscript\t400\t600\t.\t+\t.\t'
        'gene_id "G2"; gene_name "MIR"; transcript_id "TX2"; '
        'gene_type "miRNA";\n'
        'chr22\ttest\texon\t400\t600\t.\t+\t.\t'
        'gene_id "G2"; gene_name "MIR"; transcript_id "TX2"; '
        'gene_type "miRNA";\n'
        'chr22\ttest\ttranscript\t700\t900\t.\t+\t.\t'
        'gene_id "G3"; gene_name "PSEUDO"; transcript_id "TX3"; '
        'gene_biotype "processed_pseudogene";\n'
        'chr22\ttest\texon\t700\t900\t.\t+\t.\t'
        'gene_id "G3"; gene_name "PSEUDO"; transcript_id "TX3"; '
        'gene_biotype "processed_pseudogene";\n'
    )
    region = ("chr22", 1, 1000)
    expected = {
        "protein_coding": ["CODING"],
        "mirna": ["MIR"],
        "pseudogene": ["PSEUDO"],
        "noncoding": ["MIR", "PSEUDO"],
        "all": ["CODING", "MIR", "PSEUDO"],
        "type:processed_pseudogene": ["PSEUDO"],
        "protein_coding,mirna": ["CODING", "MIR"],
    }
    for expression, gene_names in expected.items():
        filters = normalize_reference_gene_types(expression)
        reference = read_reference_region(path, region, filters)
        assert list(reference["genes"]) == gene_names


def test_reference_gene_type_filter_errors_are_clear(tmp_path):
    with pytest.raises(ValueError, match="Unknown reference gene type"):
        normalize_reference_gene_types("not_a_type")
