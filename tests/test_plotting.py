from collections import OrderedDict

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from transcript_isoform_plotter.plotting import (
    AFFECTED_RED,
    EXON_HEIGHT,
    UTR_HEIGHT,
    add_direction_chevrons,
    create_gene_figure,
    create_reference_figure,
)


def test_chevrons_share_one_coordinate_grid_across_transcripts():
    fig, ax = plt.subplots()
    add_direction_chevrons(ax, 0.2, 5.2, 1, "+", spacing=1, origin=0)
    add_direction_chevrons(ax, 2.2, 8.2, 0, "+", spacing=1, origin=0)
    first_tips = [segment[1][0] for segment in ax.collections[0].get_segments()[::2]]
    second_tips = [segment[1][0] for segment in ax.collections[1].get_segments()[::2]]
    assert first_tips == [1, 2, 3, 4, 5]
    assert second_tips == [3, 4, 5, 6, 7, 8]
    plt.close(fig)


def test_isoform_plot_draws_utr_narrow_and_cds_tall():
    gene = {
        "chrom": "chr22",
        "transcripts": OrderedDict(
            [
                (
                    "TX1",
                    {
                        "label": "TX1",
                        "transcript_id": "TX1",
                        "strand": "+",
                        "exons": [(100, 200), (400, 500)],
                        "cds": [(150, 200), (400, 450)],
                        "span": (100, 500),
                        "tags": set(),
                        "fallback_span": False,
                    },
                )
            ]
        ),
    }
    fig = create_gene_figure("GENE1", gene, "kb")
    heights = [patch.get_height() for patch in fig.axes[0].patches]
    assert heights.count(UTR_HEIGHT) == 2
    assert heights.count(EXON_HEIGHT) == 2
    plt.close(fig)


def test_splice_boundary_is_one_full_panel_band():
    transcripts = OrderedDict()
    for transcript_id, exons in [
        ("TX1", [(100, 200), (400, 500)]),
        ("TX2", [(100, 200), (600, 700)]),
    ]:
        transcripts[transcript_id] = {
            "label": transcript_id,
            "transcript_id": transcript_id,
            "strand": "+",
            "exons": exons,
            "cds": [],
            "span": (exons[0][0], exons[-1][1]),
            "tags": set(),
            "fallback_span": False,
        }
    annotation = [
        {"exon": (400, 450), "boundary": (190, 410)}
    ]
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": transcripts},
        "kb",
        splice_annotations=annotation,
    )
    boundary_bands = [
        patch for patch in fig.axes[0].patches if patch.get_alpha() == 0.55
    ]
    assert len(boundary_bands) == 1
    assert boundary_bands[0].get_y() == 0
    assert boundary_bands[0].get_height() == 1
    assert abs(boundary_bands[0].get_x() - 0.190) < 1e-9
    assert abs(boundary_bands[0].get_width() - 0.221) < 1e-9
    plt.close(fig)


def test_narrow_boundary_keeps_literal_width():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 100100)],
        "cds": [],
        "span": (100, 100100),
        "tags": set(),
        "fallback_span": False,
    }
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
        splice_annotations=[{"exon": (100, 200), "boundary": (500, 501)}],
    )
    boundary_band = next(
        patch for patch in fig.axes[0].patches if patch.get_alpha() == 0.55
    )
    assert abs(boundary_band.get_x() - 0.500) < 1e-9
    assert abs(boundary_band.get_width() - 0.002) < 1e-9
    plt.close(fig)


def test_boundary_endpoints_snap_to_exon_boundaries_within_tolerance():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 100100), (200000, 200071)],
        "cds": [],
        "span": (100, 200071),
        "tags": set(),
        "fallback_span": False,
    }
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
        splice_annotations=[
            {"exon": (200000, 200071), "boundary": (199998, 200073)}
        ],
    )
    patches = fig.axes[0].patches
    boundary_band = next(patch for patch in patches if patch.get_alpha() == 0.55)
    matching_exon = next(
        patch
        for patch in patches
        if abs(patch.get_x() - 200.0) < 1e-9 and patch.get_alpha() is None
    )
    boundary_right = boundary_band.get_x() + boundary_band.get_width()
    exon_right = matching_exon.get_x() + matching_exon.get_width()
    assert abs(boundary_band.get_x() - 200.0) < 1e-9
    assert abs(boundary_right - exon_right) < 1e-9
    plt.close(fig)


def test_isoform_x_axis_has_five_percent_flanks():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 1100)],
        "cds": [],
        "span": (100, 1100),
        "tags": set(),
        "fallback_span": False,
    }
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
    )
    left, right = fig.axes[0].get_xlim()
    assert abs(left - 0.050) < 1e-9
    assert abs(right - 1.150) < 1e-9
    plt.close(fig)


def test_only_intersection_with_known_exon_is_colored_red():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 200)],
        "cds": [(130, 190)],
        "span": (100, 200),
        "tags": set(),
        "fallback_span": False,
    }
    annotation = [{"exon": (150, 175), "boundary": None}]
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
        splice_annotations=annotation,
    )
    red_patches = [
        patch
        for patch in fig.axes[0].patches
        if patch.get_facecolor() == to_rgba(AFFECTED_RED)
    ]
    # One narrow exon overlay and one taller CDS overlay, both restricted to
    # the supplied 150-175 interval rather than the complete 100-200 exon.
    assert len(red_patches) == 2
    assert all(abs(patch.get_x() - 0.150) < 1e-9 for patch in red_patches)
    assert all(abs(patch.get_width() - 0.026) < 1e-9 for patch in red_patches)
    plt.close(fig)


def test_endpoint_match_within_tolerance_colors_full_exon():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 200)],
        "cds": [(130, 190)],
        "span": (100, 200),
        "tags": set(),
        "fallback_span": False,
    }
    annotation = [{"exon": (98, 202), "boundary": None}]
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
        splice_annotations=annotation,
        splice_tolerance=2,
    )
    red_patches = [
        patch
        for patch in fig.axes[0].patches
        if patch.get_facecolor() == to_rgba(AFFECTED_RED)
    ]
    assert len(red_patches) == 2
    utr_patch = next(patch for patch in red_patches if patch.get_height() == UTR_HEIGHT)
    cds_patch = next(patch for patch in red_patches if patch.get_height() == EXON_HEIGHT)
    assert abs(utr_patch.get_x() - 0.100) < 1e-9
    assert abs(utr_patch.get_width() - 0.101) < 1e-9
    assert abs(cds_patch.get_x() - 0.130) < 1e-9
    assert abs(cds_patch.get_width() - 0.061) < 1e-9
    plt.close(fig)


def test_short_exon_and_full_match_keep_literal_proportional_width():
    transcript = {
        "label": "TX1",
        "transcript_id": "TX1",
        "strand": "+",
        "exons": [(100, 140), (100000, 100100)],
        "cds": [],
        "span": (100, 100100),
        "tags": set(),
        "fallback_span": False,
    }
    annotation = [{"exon": (99, 140), "boundary": None}]
    fig = create_gene_figure(
        "GENE1",
        {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
        "kb",
        splice_annotations=annotation,
    )
    patches = fig.axes[0].patches
    blue_exon = next(
        patch
        for patch in patches
        if patch.get_x() == 0.1 and patch.get_facecolor() != to_rgba(AFFECTED_RED)
    )
    red_exon = next(
        patch
        for patch in patches
        if patch.get_facecolor() == to_rgba(AFFECTED_RED)
    )
    assert red_exon.get_x() == blue_exon.get_x()
    assert abs(red_exon.get_width() - blue_exon.get_width()) < 1e-12
    assert abs(blue_exon.get_width() - 0.041) < 1e-9
    plt.close(fig)


def test_no_red_for_nonoverlap_or_transcript_span_fallback():
    annotations = [{"exon": (300, 350), "boundary": None}]
    for fallback in (False, True):
        transcript = {
            "label": "TX1",
            "transcript_id": "TX1",
            "strand": "+",
            "exons": [(100, 200)],
            "cds": [],
            "span": (100, 200),
            "tags": set(),
            "fallback_span": fallback,
        }
        # For the fallback case, use an interval inside the transcript span;
        # it must still remain uncolored because no exon is known.
        selected_annotations = (
            [{"exon": (120, 150), "boundary": None}]
            if fallback
            else annotations
        )
        fig = create_gene_figure(
            "GENE1",
            {"chrom": "chr1", "transcripts": OrderedDict([("TX1", transcript)])},
            "kb",
            splice_annotations=selected_annotations,
        )
        assert not any(
            patch.get_facecolor() == to_rgba(AFFECTED_RED)
            for patch in fig.axes[0].patches
        )
        plt.close(fig)


def test_reference_plot_uses_two_strand_rows_gene_labels_and_unit_scale():
    reference = {
        "chrom": "chr22",
        "start": 100,
        "end": 2100,
        "genes": OrderedDict(
            [
                (
                    "FORWARD_GENE",
                    {
                        "transcript_id": "NM_1",
                        "label": "NM_1",
                        "strand": "+",
                        "exons": [(150, 300), (500, 700)],
                        "span": (150, 700),
                        "tags": set(),
                        "fallback_span": False,
                    },
                ),
                (
                    "REVERSE_GENE",
                    {
                        "transcript_id": "NM_2",
                        "label": "NM_2",
                        "strand": "-",
                        "exons": [(1200, 1400), (1700, 2000)],
                        "span": (1200, 2000),
                        "tags": set(),
                        "fallback_span": False,
                    },
                ),
            ]
        ),
    }
    fig = create_reference_figure(reference, "kb")
    ax = fig.axes[0]
    text = {item.get_text(): item for item in ax.texts}
    assert text["FORWARD_GENE"].get_position()[1] > 1
    assert text["REVERSE_GENE"].get_position()[1] < 0
    assert any(item.get_text().endswith(" kb") for item in ax.texts)
    assert [tick.get_text() for tick in ax.get_yticklabels()] == [
        "Forward (+)",
        "Reverse (−)",
    ]
    plt.close(fig)

    fig_mb = create_reference_figure(reference, "Mb")
    assert any(item.get_text().endswith(" Mb") for item in fig_mb.axes[0].texts)
    plt.close(fig_mb)
