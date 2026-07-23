"""GTF parsing utilities."""

from __future__ import annotations

import gzip
import re
from collections import OrderedDict
from pathlib import Path
from typing import TextIO


def parse_attributes(text: str) -> dict[str, str]:
    """Parse GTF column 9 into a dictionary."""
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def open_gtf(path: Path) -> TextIO:
    """Open an uncompressed or gzip-compressed GTF as text."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def _new_transcript(attrs: dict, strand: str, tags=None) -> dict:
    transcript_id = attrs["transcript_id"]
    return {
        # Use the stable transcript accession in plots (for example ENST...),
        # rather than provider-specific display names such as SEPTIN3-207.
        "label": transcript_id,
        "transcript_id": transcript_id,
        "strand": strand,
        "exons": [],
        "cds": [],
        "span": None,
        "tags": set(tags or []),
        "fallback_span": False,
    }


def _finalize_transcript(transcript_id: str, transcript: dict) -> None:
    if transcript["exons"]:
        transcript["exons"].sort()
        transcript["cds"].sort()
    elif transcript["span"]:
        transcript["exons"] = [transcript["span"]]
        transcript["fallback_span"] = True
    else:
        raise ValueError(f"Transcript {transcript_id} has no coordinates")


def read_genes(path: Path, gene_queries: list[str]) -> OrderedDict[str, dict]:
    """Read requested genes and their transcript models in one GTF pass.

    A query can match ``gene_name`` or ``gene_id`` (case-insensitively). GTFs
    containing only transcript features are supported by treating each
    transcript span as one block.
    """
    query_lookup = {query.casefold(): query for query in gene_queries}
    genes: OrderedDict[str, dict] = OrderedDict(
        (query, {"chrom": None, "transcripts": OrderedDict()})
        for query in gene_queries
    )

    with open_gtf(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}, line {line_number}: expected 9 GTF columns")
            feature = fields[2]
            if feature not in {"transcript", "exon", "CDS"}:
                continue

            attrs = parse_attributes(fields[8])
            identifiers = [attrs.get("gene_name", ""), attrs.get("gene_id", "")]
            matched_query = next(
                (query_lookup[value.casefold()] for value in identifiers
                 if value.casefold() in query_lookup),
                None,
            )
            if matched_query is None:
                continue

            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            gene = genes[matched_query]
            if gene["chrom"] not in (None, chrom):
                raise ValueError(
                    f"Gene {matched_query!r} occurs on multiple chromosomes: "
                    f"{gene['chrom']} and {chrom}"
                )
            gene["chrom"] = chrom
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            transcript = gene["transcripts"].setdefault(
                transcript_id,
                _new_transcript(attrs, strand, re.findall(r'tag\s+"([^"]+)"', fields[8])),
            )
            transcript["tags"].update(re.findall(r'tag\s+"([^"]+)"', fields[8]))
            if feature == "exon":
                transcript["exons"].append((start, end))
            elif feature == "CDS":
                transcript["cds"].append((start, end))
            else:
                transcript["span"] = (start, end)

    missing = [query for query, gene in genes.items() if not gene["transcripts"]]
    if missing:
        available = ", ".join(missing)
        raise ValueError(f"Gene name/ID not found in {path}: {available}")

    for gene in genes.values():
        for transcript_id, transcript in gene["transcripts"].items():
            _finalize_transcript(transcript_id, transcript)
    return genes


def parse_region(text: str) -> tuple[str, int, int]:
    """Parse ``chrom:start-end``; commas in coordinates are allowed."""
    match = re.fullmatch(r"([^:]+):(\d[\d,]*)-(\d[\d,]*)", text.strip())
    if not match:
        raise ValueError(
            "--region must use chrom:start-end, for example chr22:41900000-42000000"
        )
    chrom = match.group(1)
    start = int(match.group(2).replace(",", ""))
    end = int(match.group(3).replace(",", ""))
    if start >= end:
        raise ValueError("--region start must be smaller than its end")
    return chrom, start, end


def _representative_rank(transcript: dict) -> tuple:
    """Rank curated/canonical annotations before length as a final fallback."""
    accession = transcript["transcript_id"].split(".", 1)[0]
    tags = transcript["tags"]
    # Prefer the full transcript bounds for length ranking. In region mode,
    # exon records entirely outside the requested interval are intentionally
    # not loaded because they cannot appear on the plot.
    bounds = transcript["span"] or (
        transcript["exons"][0][0],
        transcript["exons"][-1][1],
    )
    span = bounds[1] - bounds[0] + 1
    return (
        accession.startswith(("NM_", "NR_")),
        "MANE_Select" in tags,
        "Ensembl_canonical" in tags,
        "basic" in tags,
        span,
    )


def read_reference_region(path: Path, region: tuple[str, int, int]) -> dict:
    """Select one best curated/canonical transcript per gene in a region."""
    chrom, region_start, region_end = region
    genes: OrderedDict[str, dict] = OrderedDict()

    with open_gtf(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}, line {line_number}: expected 9 GTF columns")
            if fields[0] != chrom or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            start, end = int(fields[3]), int(fields[4])
            if end < region_start or start > region_end:
                continue
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            gene_label = attrs.get("gene_name") or attrs.get("gene_id")
            if not gene_label:
                continue
            gene = genes.setdefault(gene_label, {"transcripts": OrderedDict()})
            tags = re.findall(r'tag\s+"([^"]+)"', fields[8])
            transcript = gene["transcripts"].setdefault(
                transcript_id, _new_transcript(attrs, fields[6], tags)
            )
            transcript["tags"].update(tags)
            if fields[2] == "exon":
                transcript["exons"].append((start, end))
            elif fields[2] == "CDS":
                transcript["cds"].append((start, end))
            else:
                transcript["span"] = (start, end)

    if not genes:
        raise ValueError(
            f"No transcript annotations found in {chrom}:{region_start}-{region_end}"
        )

    selected = OrderedDict()
    for gene_label, gene in genes.items():
        valid = []
        for transcript_id, transcript in gene["transcripts"].items():
            try:
                _finalize_transcript(transcript_id, transcript)
            except ValueError:
                continue
            valid.append(transcript)
        if valid:
            selected[gene_label] = max(valid, key=_representative_rank)
    if not selected:
        raise ValueError(
            f"No usable transcript annotations found in {chrom}:{region_start}-{region_end}"
        )
    return {
        "chrom": chrom,
        "start": region_start,
        "end": region_end,
        "genes": selected,
    }
