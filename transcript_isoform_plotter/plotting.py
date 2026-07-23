"""Editable vector plotting for transcript isoforms."""

from __future__ import annotations

import math
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator

from .splice import matching_annotations

TRACK_BLUE = "#0000cc"
AFFECTED_RED = "#e31a1c"
BOUNDARY_BLUE = "#9bd7f0"
FALLBACK_ORANGE = "#e68613"
EXON_HEIGHT = 0.23
UTR_HEIGHT = 0.10


def _affected_bounds(genomic_exon, annotation, tolerance):
    """Return affected bounds and whether both exon endpoints matched."""
    exon_start, exon_end = genomic_exon
    region_start, region_end = annotation["exon"]
    endpoint_match = (
        abs(exon_start - region_start) <= tolerance
        and abs(exon_end - region_end) <= tolerance
    )
    if endpoint_match:
        return exon_start, exon_end, True
    return max(exon_start, region_start), min(exon_end, region_end), False


def _snap_coordinate(coordinate, candidates, tolerance):
    """Snap to the nearest candidate within tolerance, otherwise preserve it."""
    nearby = [
        candidate
        for candidate in candidates
        if abs(candidate - coordinate) <= tolerance
    ]
    return (
        min(nearby, key=lambda value: (abs(value - coordinate), value))
        if nearby
        else coordinate
    )


def add_direction_chevrons(
    ax,
    left,
    right,
    y,
    strand,
    color=TRACK_BLUE,
    spacing=None,
    origin=0,
):
    """Draw fine chevrons on a shared genomic grid.

    ``spacing`` and ``origin`` are set once per plot so every transcript uses
    identical arrow spacing and alignment, regardless of transcript length.
    """
    span = right - left
    if span <= 0:
        return
    spacing = spacing or span / 90
    if spacing <= 0:
        return
    half_width = spacing * 0.14
    half_height = 0.056
    segments = []
    first_index = math.floor((left - origin) / spacing) + 1
    x = origin + first_index * spacing
    while x < right:
        side = x - half_width if strand == "+" else x + half_width
        segments.extend(
            [[(side, y - half_height), (x, y)],
             [(side, y + half_height), (x, y)]]
        )
        x += spacing
    ax.add_collection(
        LineCollection(segments, colors=color, linewidths=0.28, capstyle="butt", zorder=2)
    )


def _tick_step(span: float, unit: str) -> float:
    candidates = (
        [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]
        if unit == "kb"
        else [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
    )
    return min(candidates, key=lambda step: abs(span / step - 7))


def _nice_scale_length(span: float) -> float:
    """Return a 1/2/5 scale-bar length near eight percent of the view."""
    target = max(span * 0.08, 1e-12)
    power = 10 ** math.floor(math.log10(target))
    choices = [power, 2 * power, 5 * power, 10 * power]
    return min(choices, key=lambda value: abs(value - target))


def create_gene_figure(
    gene_name: str,
    gene: dict,
    unit: str,
    splice_annotations=None,
    splice_tolerance: int = 2,
):
    """Create one gene/one page transcript figure."""
    transcripts: OrderedDict = gene["transcripts"]
    count = len(transcripts)
    fig_height = max(1.8, 0.38 * count + 1.25)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    divisor = 1_000 if unit == "kb" else 1_000_000
    all_exons = [exon for tx in transcripts.values() for exon in tx["exons"]]
    known_exons = [
        exon
        for transcript in transcripts.values()
        if not transcript["fallback_span"]
        for exon in transcript["exons"]
    ]
    genomic_min = min(start for start, _ in all_exons) / divisor
    genomic_max = max(end for _, end in all_exons) / divisor
    genomic_span = genomic_max - genomic_min
    chevron_spacing = genomic_span / 90

    # A splice boundary is a gene-level event. Draw it once as a continuous
    # vertical band spanning the full transcript panel rather than repeating a
    # separate rectangle on every transcript row.
    for annotation in splice_annotations or []:
        if annotation.get("boundary") is None:
            continue
        boundary_start_bp, boundary_end_bp = annotation["boundary"]
        boundary_start_bp = _snap_coordinate(
            boundary_start_bp,
            [start for start, _ in known_exons],
            splice_tolerance,
        )
        boundary_end_bp = _snap_coordinate(
            boundary_end_bp,
            [end for _, end in known_exons],
            splice_tolerance,
        )
        boundary_start = boundary_start_bp / divisor
        boundary_end = (boundary_end_bp + 1) / divisor
        ax.axvspan(
            boundary_start,
            boundary_end,
            ymin=0,
            ymax=1,
            facecolor=BOUNDARY_BLUE,
            edgecolor="none",
            alpha=0.55,
            zorder=0.5,
        )

    for row, transcript in enumerate(transcripts.values()):
        y = count - 1 - row
        track_color = FALLBACK_ORANGE if transcript["fallback_span"] else TRACK_BLUE
        genomic_exons = transcript["exons"]
        exons = [(start / divisor, end / divisor) for start, end in genomic_exons]
        ax.hlines(y, exons[0][0], exons[-1][1], color=track_color, lw=0.75, zorder=1)
        add_direction_chevrons(
            ax,
            exons[0][0],
            exons[-1][1],
            y,
            transcript["strand"],
            track_color,
            chevron_spacing,
            genomic_min,
        )

        matched_annotations = []
        for genomic_exon, (start, end) in zip(genomic_exons, exons):
            matches = (
                matching_annotations(genomic_exon, splice_annotations, splice_tolerance)
                if splice_annotations and not transcript["fallback_span"]
                else []
            )
            ax.add_patch(
                Rectangle(
                    (
                        start,
                        y
                        - (
                            EXON_HEIGHT
                            if transcript["fallback_span"]
                            else UTR_HEIGHT
                        )
                        / 2,
                    ),
                    end - start + 1 / divisor,
                    EXON_HEIGHT if transcript["fallback_span"] else UTR_HEIGHT,
                    facecolor=track_color,
                    edgecolor="white" if transcript["fallback_span"] else "none",
                    linewidth=0.25,
                    hatch="//" if transcript["fallback_span"] else None,
                    zorder=3,
                )
            )
            for match in matches:
                matched_item = (match, genomic_exon)
                if matched_item not in matched_annotations:
                    matched_annotations.append(matched_item)
                affected_start, affected_end, _ = _affected_bounds(
                    genomic_exon, match, splice_tolerance
                )
                if affected_start <= affected_end:
                    affected_width = (affected_end - affected_start + 1) / divisor
                    ax.add_patch(
                        Rectangle(
                            (
                                affected_start / divisor,
                                y - UTR_HEIGHT / 2,
                            ),
                            affected_width,
                            UTR_HEIGHT,
                            facecolor=AFFECTED_RED,
                            edgecolor="none",
                            zorder=3.2,
                        )
                    )

        # Draw coding segments over the narrow exon base. Exonic sequence
        # outside these CDS intervals remains visible as narrow UTR.
        for cds_start_bp, cds_end_bp in transcript.get("cds", []):
            cds_start = cds_start_bp / divisor
            cds_end = cds_end_bp / divisor
            ax.add_patch(
                Rectangle(
                    (cds_start, y - EXON_HEIGHT / 2),
                    cds_end - cds_start + 1 / divisor,
                    EXON_HEIGHT,
                    facecolor=track_color,
                    edgecolor="none",
                    zorder=3.5,
                )
            )
            for match, genomic_exon in matched_annotations:
                if cds_start_bp > genomic_exon[1] or cds_end_bp < genomic_exon[0]:
                    continue
                (
                    affected_exon_start,
                    affected_exon_end,
                    _,
                ) = _affected_bounds(
                    genomic_exon, match, splice_tolerance
                )
                affected_start = max(cds_start_bp, affected_exon_start)
                affected_end = min(cds_end_bp, affected_exon_end)
                if affected_start <= affected_end:
                    affected_width = (affected_end - affected_start + 1) / divisor
                    ax.add_patch(
                        Rectangle(
                            (
                                affected_start / divisor,
                                y - EXON_HEIGHT / 2,
                            ),
                            affected_width,
                            EXON_HEIGHT,
                            facecolor=AFFECTED_RED,
                            edgecolor="none",
                            zorder=3.6,
                        )
                    )

    labels = [tx["label"] for tx in transcripts.values()]
    ax.set_yticks(range(count - 1, -1, -1), labels=labels)
    ax.set_ylim(-0.62, count - 0.22)
    padding = genomic_span * 0.05
    ax.set_xlim(genomic_min - padding, genomic_max + padding)
    tick_step = _tick_step(genomic_span + 2 * padding, unit)
    ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    decimals = (1 if tick_step < 1 else 0) if unit == "kb" else 3
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.{decimals}f}"))
    ax.set_xlabel(f"{gene['chrom']} coordinate ({unit})", fontsize=10, labelpad=3)
    ax.set_title(gene_name, fontsize=11, pad=12)

    bar_length = 1 if unit == "kb" else 0.001
    if genomic_span >= bar_length:
        bar_right = genomic_max
        bar_left = bar_right - bar_length
        bar_y = count - 0.13
        cap = 0.055
        ax.plot([bar_left, bar_right], [bar_y, bar_y], color="0.45", lw=0.8, clip_on=False)
        ax.plot([bar_left, bar_left], [bar_y - cap, bar_y + cap], color="0.45", lw=0.8, clip_on=False)
        ax.plot([bar_right, bar_right], [bar_y - cap, bar_y + cap], color="0.45", lw=0.8, clip_on=False)
        ax.text((bar_left + bar_right) / 2, bar_y + 0.075, "1 kb", ha="center", fontsize=8)

    legend_handles = []
    if splice_annotations:
        legend_handles.append(
            Patch(facecolor=AFFECTED_RED, label="Specified affected region")
        )
        if any(item.get("boundary") is not None for item in splice_annotations):
            legend_handles.append(
                Patch(facecolor=BOUNDARY_BLUE, alpha=0.55, label="Splicing boundary"),
            )
    if any(tx["fallback_span"] for tx in transcripts.values()):
        legend_handles.append(
            Patch(
                facecolor=FALLBACK_ORANGE,
                edgecolor="white",
                hatch="//",
                label="Transcript span only (exons unavailable)",
            )
        )
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="lower left",
            bbox_to_anchor=(0, 1.015),
            ncol=min(3, len(legend_handles)),
            frameon=False,
            fontsize=8,
            handlelength=1.2,
            columnspacing=1.5,
            borderaxespad=0,
        )

    ax.tick_params(axis="y", length=0, labelsize=9, colors=TRACK_BLUE, pad=5)
    ax.tick_params(axis="x", labelsize=8, length=3, width=0.7, pad=2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.7)
    fig.tight_layout()
    return fig


def create_reference_figure(
    reference: dict,
    unit: str,
    splice_annotations=None,
    splice_tolerance: int = 2,
):
    """Create a reference track with one representative transcript per gene."""
    genes = reference["genes"]
    fig, ax = plt.subplots(figsize=(12, 3.4))
    divisor = 1_000 if unit == "kb" else 1_000_000
    region_min = reference["start"] / divisor
    region_max = reference["end"] / divisor
    span = region_max - region_min
    minimum_feature_width = max(span * 0.0015, 1 / divisor)
    chevron_spacing = span / 90
    label_lane_ends = {"+": [float("-inf")] * 4, "-": [float("-inf")] * 4}

    for gene_label, transcript in genes.items():
        # Reference genes share two structural rows, separated by strand.
        y = 1.0 if transcript["strand"] == "+" else 0.0
        track_color = FALLBACK_ORANGE if transcript["fallback_span"] else TRACK_BLUE
        genomic_exons = transcript["exons"]
        exons = [(start / divisor, end / divisor) for start, end in genomic_exons]
        ax.hlines(y, exons[0][0], exons[-1][1], color=track_color, lw=0.7, zorder=1)
        add_direction_chevrons(
            ax,
            exons[0][0],
            exons[-1][1],
            y,
            transcript["strand"],
            track_color,
            chevron_spacing,
            region_min,
        )
        for genomic_exon, (start, end) in zip(genomic_exons, exons):
            matches = (
                matching_annotations(genomic_exon, splice_annotations, splice_tolerance)
                if splice_annotations
                else []
            )
            for match in matches:
                if match.get("boundary") is None:
                    continue
                boundary_start, boundary_end = match["boundary"]
                boundary_start /= divisor
                boundary_end /= divisor
                width = max(boundary_end - boundary_start, span * 0.0012)
                center = (boundary_start + boundary_end) / 2
                ax.add_patch(
                    Rectangle(
                        (center - width / 2, y - EXON_HEIGHT * 0.78),
                        width,
                        EXON_HEIGHT * 1.56,
                        facecolor=BOUNDARY_BLUE,
                        edgecolor="none",
                        alpha=0.55,
                        zorder=0.5,
                    )
                )
            ax.add_patch(
                Rectangle(
                    (start, y - EXON_HEIGHT / 2),
                    max(end - start + 1 / divisor, minimum_feature_width),
                    EXON_HEIGHT,
                    facecolor=AFFECTED_RED if matches else track_color,
                    edgecolor="white" if transcript["fallback_span"] else "none",
                    linewidth=0.25,
                    hatch="//" if transcript["fallback_span"] else None,
                    zorder=3,
                )
            )

        visible_left = max(exons[0][0], region_min)
        visible_right = min(exons[-1][1], region_max)
        if visible_left <= visible_right:
            # Approximate label width in genomic units, then greedily place
            # nearby names into separate label lanes to reduce collisions.
            label_width = span * (0.006 * len(gene_label) + 0.008)
            center = (visible_left + visible_right) / 2
            center = min(
                max(center, region_min + label_width / 2),
                region_max - label_width / 2,
            )
            label_left = center - label_width / 2
            label_right = center + label_width / 2
            lanes = label_lane_ends[transcript["strand"]]
            lane = next(
                (
                    index
                    for index, previous_right in enumerate(lanes)
                    if label_left > previous_right + span * 0.004
                ),
                min(range(len(lanes)), key=lanes.__getitem__),
            )
            lanes[lane] = label_right
            label_y = (
                y + 0.20 + lane * 0.12
                if transcript["strand"] == "+"
                else y - 0.20 - lane * 0.12
            )
            ax.text(
                center,
                label_y,
                gene_label,
                color=TRACK_BLUE,
                fontsize=7.0,
                fontstyle="italic",
                ha="center",
                va="bottom" if transcript["strand"] == "+" else "top",
                clip_on=True,
            )

    ax.set_yticks([1, 0], labels=["Forward (+)", "Reverse (−)"])
    ax.set_ylim(-0.72, 1.82)
    ax.set_xlim(region_min, region_max)
    tick_step = _tick_step(span, unit)
    ax.xaxis.set_major_locator(MultipleLocator(tick_step))
    decimals = (1 if tick_step < 1 else 0) if unit == "kb" else 3
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:.{decimals}f}")
    )
    ax.set_xlabel(f"{reference['chrom']} coordinate ({unit})", fontsize=10, labelpad=3)
    ax.set_title("Representative curated reference transcripts", fontsize=11, pad=12)

    scale_length = _nice_scale_length(span)
    scale_right = region_max - span * 0.015
    scale_left = scale_right - scale_length
    scale_y = 1.68
    cap = 0.045
    ax.plot(
        [scale_left, scale_right],
        [scale_y, scale_y],
        color="0.35",
        lw=0.8,
        clip_on=False,
    )
    ax.plot(
        [scale_left, scale_left],
        [scale_y - cap, scale_y + cap],
        color="0.35",
        lw=0.8,
        clip_on=False,
    )
    ax.plot(
        [scale_right, scale_right],
        [scale_y - cap, scale_y + cap],
        color="0.35",
        lw=0.8,
        clip_on=False,
    )
    ax.text(
        (scale_left + scale_right) / 2,
        scale_y + 0.06,
        f"{scale_length:g} {unit}",
        ha="center",
        va="bottom",
        fontsize=8,
    )

    legend_handles = []
    if splice_annotations:
        legend_handles.append(Patch(facecolor=AFFECTED_RED, label="Affected exon"))
        if any(item.get("boundary") is not None for item in splice_annotations):
            legend_handles.append(
                Patch(facecolor=BOUNDARY_BLUE, alpha=0.55, label="Splicing boundary"),
            )
    if any(tx["fallback_span"] for tx in genes.values()):
        legend_handles.append(
            Patch(
                facecolor=FALLBACK_ORANGE,
                edgecolor="white",
                hatch="//",
                label="Transcript span only (exons unavailable)",
            )
        )
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="lower left",
            bbox_to_anchor=(0, 1.015),
            ncol=min(3, len(legend_handles)),
            frameon=False,
            fontsize=8,
            handlelength=1.2,
            columnspacing=1.5,
            borderaxespad=0,
        )

    ax.tick_params(axis="y", length=0, labelsize=8, colors="0.25", pad=7)
    ax.tick_params(axis="x", labelsize=8, length=3, width=0.7, pad=2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.7)
    fig.tight_layout()
    return fig


def write_pdf(
    genes: OrderedDict,
    output: Path,
    unit: str,
    splice_annotations=None,
    splice_tolerance: int = 2,
    gene_splice_annotations=None,
) -> None:
    """Write one editable PDF page per requested gene."""
    settings = {"pdf.fonttype": 42, "ps.fonttype": 42}
    with mpl.rc_context(settings), PdfPages(output) as pdf:
        metadata = {"Title": "Transcript isoform structures", "Creator": "gtf-isoform-plotter"}
        pdf.infodict().update(metadata)
        for gene_name, gene in genes.items():
            annotations = list(splice_annotations or [])
            if gene_splice_annotations:
                annotations.extend(gene_splice_annotations.get(gene_name, []))
            fig = create_gene_figure(
                gene_name, gene, unit, annotations or None, splice_tolerance
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def safe_gene_filename(gene_name: str) -> str:
    """Convert a gene name/ID into a portable filename stem."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._")
    return stem or "gene"


def write_gene_pdfs(
    genes: OrderedDict,
    output_dir: Path,
    unit: str,
    splice_annotations=None,
    splice_tolerance: int = 2,
    gene_splice_annotations=None,
) -> OrderedDict:
    """Write one editable, gene-named PDF for each isoform model."""
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = OrderedDict()
    used_names = {}
    for gene_name, gene in genes.items():
        filename = f"{safe_gene_filename(gene_name)}.pdf"
        folded = filename.casefold()
        if folded in used_names:
            raise ValueError(
                f"Gene names {used_names[folded]!r} and {gene_name!r} produce "
                f"the same output filename {filename!r}"
            )
        used_names[folded] = gene_name
        output = output_dir / filename
        write_pdf(
            OrderedDict([(gene_name, gene)]),
            output,
            unit,
            splice_annotations,
            splice_tolerance,
            gene_splice_annotations,
        )
        outputs[gene_name] = output
    return outputs


def write_reference_pdf(
    reference: dict,
    output: Path,
    unit: str,
    splice_annotations=None,
    splice_tolerance: int = 2,
) -> None:
    """Write an editable one-page representative reference-gene track."""
    settings = {"pdf.fonttype": 42, "ps.fonttype": 42}
    with mpl.rc_context(settings), PdfPages(output) as pdf:
        pdf.infodict().update(
            {
                "Title": "Representative curated reference transcripts",
                "Creator": "gtf-isoform-plotter",
            }
        )
        fig = create_reference_figure(
            reference, unit, splice_annotations, splice_tolerance
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
