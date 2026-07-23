from transcript_isoform_plotter.cli import main

from .test_gtf import GTF


def test_cli_writes_editable_pdf_with_optional_splice_omitted(tmp_path):
    gtf = tmp_path / "test.gtf"
    output_dir = tmp_path / "plots"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--genes", "GENE1,GENE2",
        "--output-dir", str(output_dir),
    ])
    assert result == 0
    assert (output_dir / "GENE1.pdf").read_bytes().startswith(b"%PDF")
    assert (output_dir / "GENE2.pdf").read_bytes().startswith(b"%PDF")
    assert (output_dir / "GENE1.pdf").stat().st_size > 1_000


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


def test_cli_reference_mode_accepts_all_gene_types(tmp_path):
    gtf = tmp_path / "test.gtf"
    output = tmp_path / "reference_all.pdf"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--plot-type", "reference",
        "--region", "chr22:50-600",
        "--reference-gene-type", "all",
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")


def test_cli_accepts_gene_table_in_isoform_mode(tmp_path):
    gtf = tmp_path / "test.gtf"
    gene_file = tmp_path / "genes.tsv"
    output_dir = tmp_path / "gene_table"
    gtf.write_text(GTF)
    gene_file.write_text(
        "GENE1\t100, 200\t200, 400\n"
        "GENE2\n"
    )
    result = main([
        "--gtf-file", str(gtf),
        "--plot-type", "isoforms",
        "--gene-file", str(gene_file),
        "--output-dir", str(output_dir),
    ])
    assert result == 0
    assert (output_dir / "GENE1.pdf").read_bytes().startswith(b"%PDF")
    assert (output_dir / "GENE2.pdf").read_bytes().startswith(b"%PDF")


def test_single_isoform_gene_keeps_explicit_output_filename(tmp_path):
    gtf = tmp_path / "test.gtf"
    output = tmp_path / "custom.pdf"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--genes", "GENE1",
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")


def test_cli_filters_isoforms_by_annotation_tag(tmp_path):
    gtf = tmp_path / "test.gtf"
    output = tmp_path / "basic.pdf"
    gtf.write_text(GTF)
    result = main([
        "--gtf-file", str(gtf),
        "--genes", "GENE1",
        "--transcript-filter", "basic",
        "--output", str(output),
    ])
    assert result == 0
    assert output.read_bytes().startswith(b"%PDF")
