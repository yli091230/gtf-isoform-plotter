import gzip

import pytest

from transcript_isoform_plotter.gtf import parse_region, read_genes, read_reference_region


GTF = """##format: gtf
chr22\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1"; transcript_name "TX1.1";
chr22\ttest\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\tCDS\t150\t200\t.\t+\t0\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\texon\t400\t500\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\tCDS\t400\t450\t.\t+\t0\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX1";
chr22\ttest\ttranscript\t120\t450\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2";
chr22\ttest\texon\t120\t220\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2";
chr22\ttest\texon\t350\t450\t.\t-\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "TX2";
chr1\ttest\ttranscript\t1000\t1300\t.\t+\t.\tgene_id "G2"; gene_name "GENE2"; transcript_id "TX3";
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
