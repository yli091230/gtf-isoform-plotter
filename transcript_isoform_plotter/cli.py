"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .gtf import parse_region, read_genes, read_reference_region
from .plotting import write_pdf, write_reference_pdf
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
        "--output",
        "-o",
        type=Path,
        default=Path("transcript_isoforms.pdf"),
        help="Editable PDF output (default: transcript_isoforms.pdf)",
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

    try:
        gene_splice_annotations = None
        if args.plot_type == "isoforms" and args.gene_file:
            genes, gene_splice_annotations = read_gene_table(args.gene_file)
        splice_annotations = (
            read_splice_annotations(args.splice_file) if args.splice_file else None
        )
        if args.plot_type == "isoforms":
            gene_models = read_genes(args.gtf_file, genes)
            write_pdf(
                gene_models,
                args.output,
                args.unit,
                splice_annotations,
                args.splice_tolerance,
                gene_splice_annotations,
            )
        else:
            reference = read_reference_region(args.gtf_file, parse_region(args.region))
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
            f"gene(s) to {args.output}"
        )
    else:
        print(
            f"Saved {len(reference['genes'])} representative reference transcript(s) "
            f"to {args.output}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
