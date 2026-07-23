"""GTF parsing utilities."""

from __future__ import annotations

import gzip
import re
from collections import OrderedDict
from pathlib import Path
from typing import TextIO

SUPPORTED_TRANSCRIPT_FILTERS = {
    "all",
    "basic",
    "mane",
    "canonical",
    "appris_principal",
    "protein_coding",
    "coding",
    "noncoding",
}

SUPPORTED_REFERENCE_GENE_TYPES = {
    "all",
    "protein_coding",
    "noncoding",
    "pseudogene",
    "mirna",
    "lncrna",
}


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
        "transcript_type": attrs.get("transcript_type", ""),
        "gene_type": (
            attrs.get("gene_type")
            or attrs.get("gene_biotype")
            or attrs.get("transcript_type", "")
        ),
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


def normalize_transcript_filters(text: str) -> list[str]:
    """Parse and validate comma-separated transcript filters."""
    aliases = {
        "mane_select": "mane",
        "ensembl_canonical": "canonical",
        "appris": "appris_principal",
    }
    filters = []
    for value in text.split(","):
        value = aliases.get(value.strip().casefold(), value.strip().casefold())
        if not value:
            continue
        if value not in SUPPORTED_TRANSCRIPT_FILTERS and not value.startswith("type:"):
            choices = ", ".join(sorted(SUPPORTED_TRANSCRIPT_FILTERS))
            raise ValueError(
                f"Unknown transcript filter {value!r}. Use one of {choices}, "
                "or type:TRANSCRIPT_TYPE"
            )
        if value == "type:":
            raise ValueError("type: transcript filter requires a transcript type value")
        if value not in filters:
            filters.append(value)
    if not filters:
        raise ValueError("--transcript-filter must contain at least one filter")
    return filters


def _transcript_matches_filter(transcript: dict, filter_name: str) -> bool:
    tags = transcript["tags"]
    if filter_name == "all":
        return True
    if filter_name == "basic":
        return "basic" in tags
    if filter_name == "mane":
        return "MANE_Select" in tags
    if filter_name == "canonical":
        return "Ensembl_canonical" in tags
    if filter_name == "appris_principal":
        return any(tag.startswith("appris_principal") for tag in tags)
    if filter_name == "protein_coding":
        return transcript["transcript_type"] == "protein_coding"
    if filter_name == "coding":
        return bool(transcript["cds"])
    if filter_name == "noncoding":
        return not transcript["cds"]
    if filter_name.startswith("type:"):
        requested_type = filter_name.split(":", 1)[1]
        return transcript["transcript_type"].casefold() == requested_type.casefold()
    return False


def filter_gene_transcripts(genes: OrderedDict, filters: list[str]) -> OrderedDict:
    """Filter isoforms with OR semantics while preserving GTF ordering."""
    if "all" in filters:
        return genes
    filtered = OrderedDict()
    empty_genes = []
    for gene_name, gene in genes.items():
        transcripts = OrderedDict(
            (
                transcript_id,
                transcript,
            )
            for transcript_id, transcript in gene["transcripts"].items()
            if any(
                _transcript_matches_filter(transcript, filter_name)
                for filter_name in filters
            )
        )
        if not transcripts:
            empty_genes.append(gene_name)
            continue
        filtered[gene_name] = {**gene, "transcripts": transcripts}
    if empty_genes:
        raise ValueError(
            "No transcripts passed --transcript-filter for: "
            + ", ".join(empty_genes)
        )
    return filtered


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


def normalize_reference_gene_types(text: str) -> list[str]:
    """Parse comma-separated reference gene-biotype filters."""
    aliases = {
        "protein-coding": "protein_coding",
        "protein coding": "protein_coding",
        "mirna": "mirna",
        "lincrna": "lncrna",
        "other": "noncoding",
    }
    filters = []
    for raw_value in text.split(","):
        value = raw_value.strip()
        normalized = aliases.get(value.casefold(), value.casefold())
        if not normalized:
            continue
        if (
            normalized not in SUPPORTED_REFERENCE_GENE_TYPES
            and not normalized.startswith("type:")
        ):
            choices = ", ".join(sorted(SUPPORTED_REFERENCE_GENE_TYPES))
            raise ValueError(
                f"Unknown reference gene type {value!r}. Use one of {choices}, "
                "or type:GENE_TYPE"
            )
        if normalized == "type:":
            raise ValueError("type: reference gene filter requires a gene type value")
        if normalized not in filters:
            filters.append(normalized)
    if not filters:
        raise ValueError("--reference-gene-type must contain at least one filter")
    return filters


def _matches_reference_gene_type(transcript: dict, filter_name: str) -> bool:
    gene_type = transcript.get("gene_type", "")
    normalized_type = gene_type.casefold()
    if filter_name == "all":
        return True
    if not normalized_type:
        return False
    if filter_name == "protein_coding":
        return normalized_type == "protein_coding"
    if filter_name == "noncoding":
        return normalized_type != "protein_coding"
    if filter_name == "pseudogene":
        return "pseudogene" in normalized_type
    if filter_name == "mirna":
        return normalized_type in {"mirna", "mirna_gene"}
    if filter_name == "lncrna":
        return normalized_type in {"lncrna", "lincrna"}
    if filter_name.startswith("type:"):
        return normalized_type == filter_name.split(":", 1)[1].casefold()
    return False


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


def read_reference_region(
    path: Path,
    region: tuple[str, int, int],
    gene_type_filters=None,
) -> dict:
    """Select one best curated/canonical transcript per gene in a region."""
    gene_type_filters = gene_type_filters or ["all"]
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
            if any(
                _matches_reference_gene_type(transcript, filter_name)
                for filter_name in gene_type_filters
            ):
                valid.append(transcript)
        if valid:
            selected[gene_label] = max(valid, key=_representative_rank)
    if not selected:
        raise ValueError(
            "No reference genes passed --reference-gene-type "
            f"{','.join(gene_type_filters)} in "
            f"{chrom}:{region_start}-{region_end}"
        )
    return {
        "chrom": chrom,
        "start": region_start,
        "end": region_end,
        "genes": selected,
    }
