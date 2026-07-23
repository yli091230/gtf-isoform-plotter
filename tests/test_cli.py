from transcript_isoform_plotter.cli import main

from .test_gtf import GTF


def test_cli_writes_editable_pdf_with_optional_splice_omitted(tmp_path):
    gtf = tmp_path / "test.gtf"
    output = tmp_path / "genes.pdf"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--genes", "GENE1,GENE2",
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")
    assert output.stat().st_size > 1_000


def test_cli_reference_mode_writes_region_track(tmp_path):
    gtf = tmp_path / "test.gtf"
    output = tmp_path / "reference.pdf"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--plot-type", "reference",
        "--region", "chr22:50-600",
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")


def test_cli_accepts_gene_table_in_isoform_mode(tmp_path):
    gtf = tmp_path / "test.gtf"
    gene_file = tmp_path / "genes.tsv"
    output = tmp_path / "gene_table.pdf"
    gtf.write_text(GTF)
    gene_file.write_text(
        "GENE1\t100, 200\t200, 400\n"
        "GENE2\n"
    )
    result = main([
        "--gtf-file", str(gtf),
        "--plot-type", "isoforms",
        "--gene-file", str(gene_file),
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")
