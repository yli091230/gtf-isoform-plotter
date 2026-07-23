from collections import OrderedDict

import matplotlib.pyplot as plt

from transcript_isoform_plotter.plotting import (
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
