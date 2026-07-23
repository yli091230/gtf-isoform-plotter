"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .gtf import (
    filter_gene_transcripts,
    normalize_reference_gene_types,
    normalize_transcript_filters,
    parse_region,
    read_genes,
    read_reference_region,
)
from .plotting import write_gene_pdfs, write_pdf, write_reference_pdf
from .splice import read_gene_table, read_splice_annotations


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create editable transcript isoform diagrams from a GTF file."
    )
    parser.add_argument("--gtf-file", required=True, type=Path, help="GTF or GTF.GZ file")
    parser.add_argument(
        "--plot-type",
        choices=("isoforms", "reference"),
        default="isoforms",
        help="Plot all isoforms by gene, or one representative gene track by region",
    )
    parser.add_argument(
        "--genes",
        help="For --plot-type isoforms: gene name or comma-separated names/IDs",
    )
    parser.add_argument(
        "--gene-file",
        type=Path,
        help=(
            "For --plot-type isoforms: TSV with gene in column 1 and optional "
            "affected-exon and splice-boundary pairs in columns 2 and 3"
        ),
    )
    parser.add_argument(
        "--region",
        help="For --plot-type reference: chromosome range, e.g. chr22:41900000-42000000",
    )
    parser.add_argument(
        "--transcript-filter",
        default="all",
        help=(
            "Isoform filter(s), comma-separated with OR behavior: all, basic, "
            "mane, canonical, appris_principal, protein_coding, coding, "
            "noncoding, or type:TRANSCRIPT_TYPE (default: all)"
        ),
    )
    parser.add_argument(
        "--reference-gene-type",
        default="protein_coding",
        help=(
            "Reference-mode gene type(s), comma-separated with OR behavior: "
            "all, protein_coding, noncoding, pseudogene, miRNA, lncRNA, or "
            "type:GENE_TYPE (default: protein_coding)"
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("transcript_isoforms.pdf"),
        help="Editable PDF output (default: transcript_isoforms.pdf)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Isoform output directory. Multiple genes are always written here "
            "as separate GENE.pdf files"
        ),
    )
    parser.add_argument("--unit", choices=("kb", "Mb"), default="kb")
    parser.add_argument(
        "--splice-file",
        type=Path,
        help="Optional two-column splice annotation file; omit for no highlighting",
    )
    parser.add_argument(
        "--splice-tolerance",
        type=int,
        default=2,
        metavar="BP",
        help="Exon-coordinate matching tolerance (default: 2 bp)",
    )
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.output.suffix.lower() != ".pdf":
        parser.error("--output must end in .pdf")
    if args.splice_tolerance < 0:
        parser.error("--splice-tolerance must be zero or greater")
    if args.plot_type == "isoforms":
        if bool(args.genes) == bool(args.gene_file):
            parser.error(
                "--plot-type isoforms requires exactly one of --genes or --gene-file"
            )
        if args.region:
            parser.error("--region is only valid with --plot-type reference")
        if args.genes:
            genes = [item.strip() for item in args.genes.split(",") if item.strip()]
            if not genes:
                parser.error("--genes must contain at least one gene name or ID")
            genes = list(dict.fromkeys(genes))
    else:
        if not args.region:
            parser.error("--region is required when --plot-type reference")
        if args.genes or args.gene_file:
            parser.error(
                "--genes and --gene-file are only valid with --plot-type isoforms"
            )
        if args.output_dir:
            parser.error("--output-dir is only valid with --plot-type isoforms")
        if args.transcript_filter.casefold() != "all":
            parser.error(
                "--transcript-filter is only valid with --plot-type isoforms"
            )

    try:
        transcript_filters = normalize_transcript_filters(args.transcript_filter)
        reference_gene_types = normalize_reference_gene_types(
            args.reference_gene_type
        )
        gene_splice_annotations = None
        if args.plot_type == "isoforms" and args.gene_file:
            genes, gene_splice_annotations = read_gene_table(args.gene_file)
        splice_annotations = (
            read_splice_annotations(args.splice_file) if args.splice_file else None
        )
        if args.plot_type == "isoforms":
            gene_models = read_genes(args.gtf_file, genes)
            gene_models = filter_gene_transcripts(
                gene_models, transcript_filters
            )
            if len(gene_models) > 1 or args.output_dir:
                output_dir = args.output_dir or args.output.parent
                output_paths = write_gene_pdfs(
                    gene_models,
                    output_dir,
                    args.unit,
                    splice_annotations,
                    args.splice_tolerance,
                    gene_splice_annotations,
                )
            else:
                write_pdf(
                    gene_models,
                    args.output,
                    args.unit,
                    splice_annotations,
                    args.splice_tolerance,
                    gene_splice_annotations,
                )
                output_paths = {next(iter(gene_models)): args.output}
        else:
            reference = read_reference_region(
                args.gtf_file,
                parse_region(args.region),
                reference_gene_types,
            )
            write_reference_pdf(
                reference,
                args.output,
                args.unit,
                splice_annotations,
                args.splice_tolerance,
            )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    if args.plot_type == "isoforms":
        transcript_count = sum(
            len(gene["transcripts"]) for gene in gene_models.values()
        )
        print(
            f"Saved {transcript_count} transcript model(s) for {len(gene_models)} "
            f"gene(s) to {len(output_paths)} PDF file(s)"
        )
        for gene_name, output_path in output_paths.items():
            print(f"  {gene_name}: {output_path}")
    else:
        print(
            f"Saved {len(reference['genes'])} representative reference transcript(s) "
            f"to {args.output}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
