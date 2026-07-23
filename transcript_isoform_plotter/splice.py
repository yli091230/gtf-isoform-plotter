"""Splice annotation parsing and matching."""

from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path


def _coordinate_pair(value: str, path: Path, line_number: int, column_number: int):
    """Parse one comma-separated genomic coordinate pair."""
    values = re.findall(r"-?\d+", value)
    if len(values) != 2:
        raise ValueError(
            f"{path}, line {line_number}, column {column_number}: "
            "expected two comma-separated coordinates"
        )
    return tuple(sorted((int(values[0]), int(values[1]))))


def read_splice_annotations(path: Path) -> list[dict[str, tuple[int, int]]]:
    """Read tab-separated affected-exon and splice-boundary coordinate pairs."""
    annotations = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 2:
                if not re.search(r"\d", stripped):
                    continue
                raise ValueError(
                    f"{path}, line {line_number}: expected two tab-separated columns"
                )
            if not re.search(r"\d", stripped):
                continue
            pairs = [
                _coordinate_pair(column, path, line_number, column_number)
                for column_number, column in enumerate(columns, start=1)
            ]
            if pairs:
                annotations.append({"exon": pairs[0], "boundary": pairs[1]})
    if not annotations:
        raise ValueError(f"No splice annotations were found in {path}")
    return annotations


def read_gene_table(path: Path):
    """Read genes and optional per-gene splice annotations from a TSV file.

    Column 1 is the gene name/ID. Column 2 is an optional affected-exon pair.
    Column 3 is an optional splice-boundary pair. Repeated genes are retained
    once in plotting order while their annotations are accumulated.
    """
    genes = []
    annotations_by_gene = OrderedDict()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) > 3:
                raise ValueError(
                    f"{path}, line {line_number}: expected one to three tab-separated columns"
                )
            columns += [""] * (3 - len(columns))
            gene = columns[0].strip()
            if not gene:
                raise ValueError(f"{path}, line {line_number}: gene column is empty")
            if gene.casefold() in {"gene", "gene_name", "gene_id"}:
                continue
            if gene not in annotations_by_gene:
                genes.append(gene)
                annotations_by_gene[gene] = []

            exon_text = columns[1].strip()
            boundary_text = columns[2].strip()
            if boundary_text and not exon_text:
                raise ValueError(
                    f"{path}, line {line_number}: column 3 requires an exon pair in column 2"
                )
            if exon_text:
                annotation = {
                    "exon": _coordinate_pair(exon_text, path, line_number, 2),
                    "boundary": (
                        _coordinate_pair(boundary_text, path, line_number, 3)
                        if boundary_text
                        else None
                    ),
                }
                annotations_by_gene[gene].append(annotation)

    if not genes:
        raise ValueError(f"No genes were found in {path}")
    return genes, annotations_by_gene


def matching_annotations(exon, annotations, tolerance, min_overlap=0.80):
    """Return annotations associated with an exon.

    A match is accepted when both ends agree within ``tolerance`` bp, or when
    at least ``min_overlap`` of the supplied affected-exon interval overlaps
    the GTF exon. The overlap rule supports partial exon intervals produced by
    splice-event tools rather than requiring complete GTF exon bounds.
    """
    start, end = exon
    matches = []
    for item in annotations:
        annotation_start, annotation_end = item["exon"]
        endpoint_match = (
            abs(start - annotation_start) <= tolerance
            and abs(end - annotation_end) <= tolerance
        )
        overlap = max(
            0,
            min(end, annotation_end) - max(start, annotation_start) + 1,
        )
        annotation_length = annotation_end - annotation_start + 1
        overlap_match = overlap / annotation_length >= min_overlap
        if endpoint_match or overlap_match:
            matches.append(item)
    return matches
