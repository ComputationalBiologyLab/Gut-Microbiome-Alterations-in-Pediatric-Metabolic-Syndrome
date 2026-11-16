#!/usr/bin/env python3
"""PRISMA flow diagram for the three-cohort pediatric MetS meta-analysis."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle

ROOT = Path(os.environ.get("MS_ROOT", ".")).resolve()  # project directory
FIGS = ROOT / "figures"
META = ROOT / "meta_analysis"

COUNTS = {
    "pubmed": 111,
    "bioproject": 40,
    "citation": 3,
    "identified": 154,
    "duplicates": 28,
    "screened": 126,
    "screen_excluded": 108,
    "fulltext": 18,
    "unavailable": 2,
    "cohort_articles": 16,
    "elig_excluded": 13,
    "included": 3,
}

BANNER = "#EFB33C"
SIDEBAR = "#A9C5E1"

FIG_W = 9.9
BAR_X, BAR_W = 0.10, 0.40
MAIN_X, MAIN_W = 0.78, 3.72
RIGHT_X, RIGHT_W = 5.35, 4.35
PAD = 0.22
LH_MAIN = 0.28
LH_RIGHT = 0.32
ROW_GAP = 0.58
TOP_MARGIN = 0.95
BOTTOM_MARGIN = 0.30

ROWS = [
    [
        "Records identified from PubMed and",
        f"NCBI BioProject databases (n = {COUNTS['identified']})",
        f"Duplicates removed (n = {COUNTS['duplicates']})",
    ],
    [
        "Records screened by title and",
        f"abstract (n = {COUNTS['screened']})",
    ],
    [
        "Full-text articles assessed for",
        f"eligibility (n = {COUNTS['fulltext']})",
    ],
    [
        "Articles addressing MetS in children",
        f"and adolescents (n = {COUNTS['cohort_articles']})",
    ],
    [
        "Studies included in Meta-",
        f"analysis (n = {COUNTS['included']})",
    ],
]

# right box: (anchor row index, arrow row index, title lines, reasons)
RIGHT_BOXES = [
    (
        0,
        1,
        [f"Records excluded (n = {COUNTS['screen_excluded']}), with reasons:"],
        [
            "Wrong population or sample type: 52",
            "No MetS case-control comparison: 31",
            "Review/commentary/no primary data: 14",
            "Non-16S or non-fecal microbiome: 11",
        ],
    ),
    (
        2,
        2,
        [
            "Studies excluded due to data",
            f"unavailability (n = {COUNTS['unavailable']})",
        ],
        [],
    ),
    (
        3,
        3,
        [f"Studies excluded (n = {COUNTS['elig_excluded']}), with reasons:"],
        [
            "Adult-only cohort: 4",
            "No public SRA/BioProject data: 2",
            "Shotgun-only / non-comparable 16S: 2",
            "Overlapping or duplicate cohort: 2",
            "Animal model: 1",
            "MetS arm n < 10: 1",
            "Risk-factor gradient design: 1",
        ],
    ),
]

STAGES = [("Identification", 0, 0), ("Screening", 1, 3), ("Included", 4, 4)]


def main_height(lines: list[str]) -> float:
    return max(0.88, PAD * 2 + len(lines) * LH_MAIN)


def right_height(title_lines: list[str], reasons: list[str]) -> float:
    slots = len(title_lines) + (1 + len(reasons) if reasons else 0)
    return PAD * 2 + slots * LH_RIGHT


def draw_rect(ax, x, y, w, h, face="white", edge="black", lw=1.4):
    ax.add_patch(
        Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, linewidth=lw, zorder=2)
    )


def v_arrow(ax, x, y_from, y_to):
    ax.add_patch(
        FancyArrow(
            x,
            y_from,
            0,
            y_to - y_from,
            width=0.006,
            head_width=0.11,
            head_length=0.13,
            length_includes_head=True,
            color="black",
            zorder=3,
        )
    )


def h_arrow(ax, x_from, y, x_to):
    ax.add_patch(
        FancyArrow(
            x_from,
            y,
            x_to - x_from,
            0,
            width=0.006,
            head_width=0.11,
            head_length=0.13,
            length_includes_head=True,
            color="black",
            zorder=3,
        )
    )


def write_counts_tsv() -> None:
    rows = [
        ("stage", "item", "count"),
        ("identification", "pubmed_records", COUNTS["pubmed"]),
        ("identification", "bioproject_records", COUNTS["bioproject"]),
        ("identification", "citation_records", COUNTS["citation"]),
        ("identification", "records_identified_total", COUNTS["identified"]),
        ("identification", "duplicates_removed", COUNTS["duplicates"]),
        ("screening", "records_screened", COUNTS["screened"]),
        ("screening", "records_excluded", COUNTS["screen_excluded"]),
        ("eligibility", "fulltext_assessed", COUNTS["fulltext"]),
        ("eligibility", "excluded_data_unavailable", COUNTS["unavailable"]),
        ("eligibility", "articles_pediatric_mets", COUNTS["cohort_articles"]),
        ("eligibility", "studies_excluded", COUNTS["elig_excluded"]),
        ("included", "studies_in_meta_analysis", COUNTS["included"]),
    ]
    with (META / "PRISMA_flow_counts.tsv").open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write("\t".join(str(c) for c in r) + "\n")


def make_figure() -> Path:
    main_hs = [main_height(lines) for lines in ROWS]
    right_hs = {i: right_height(t, r) for i, (_a, _ar, t, r) in enumerate(RIGHT_BOXES)}

    # Provisional top-down layout with y measured downward from 0
    tops = []
    cursor = 0.0
    for h in main_hs:
        tops.append(cursor)
        cursor += h + ROW_GAP
    content_depth = cursor - ROW_GAP

    # Keep the incoming arrow clear of a right box's bottom edge
    for i, (anchor, arrow_row, _t, _r) in enumerate(RIGHT_BOXES):
        arrow_depth = tops[arrow_row] + main_hs[arrow_row] / 2 - tops[anchor]
        right_hs[i] = max(right_hs[i], arrow_depth + 0.40)

    for i, (anchor, _arrow_row, _t, _r) in enumerate(RIGHT_BOXES):
        content_depth = max(content_depth, tops[anchor] + right_hs[i])

    fig_h = content_depth + TOP_MARGIN + BOTTOM_MARGIN + 0.55
    fig = plt.figure(figsize=(FIG_W, fig_h))
    # Axes fills the canvas so one data unit equals one inch
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    def Y(depth: float) -> float:
        return fig_h - TOP_MARGIN - 0.55 - depth

    # Header banner
    banner_h = 0.42
    banner_x = MAIN_X - 0.16
    banner_w = (RIGHT_X + RIGHT_W) - banner_x
    draw_rect(ax, banner_x, fig_h - 0.30 - banner_h, banner_w, banner_h, face=BANNER, lw=1.0)
    ax.text(
        banner_x + banner_w / 2,
        fig_h - 0.30 - banner_h / 2,
        "Identification of studies via databases and registers",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        zorder=4,
    )

    # Main boxes
    for i, lines in enumerate(ROWS):
        h = main_hs[i]
        y_top = Y(tops[i])
        y_bot = y_top - h
        draw_rect(ax, MAIN_X, y_bot, MAIN_W, h)
        first = y_top - PAD - LH_MAIN * 0.72
        for j, line in enumerate(lines):
            ax.text(
                MAIN_X + MAIN_W / 2,
                first - j * LH_MAIN,
                line,
                ha="center",
                va="center",
                fontsize=9.2,
                fontweight="bold",
                zorder=4,
            )
        if i < len(ROWS) - 1:
            v_arrow(ax, MAIN_X + MAIN_W / 2, y_bot, Y(tops[i + 1]))

    # Right boxes
    for i, (anchor, arrow_row, title_lines, reasons) in enumerate(RIGHT_BOXES):
        h = right_hs[i]
        y_top = Y(tops[anchor])
        y_bot = y_top - h
        draw_rect(ax, RIGHT_X, y_bot, RIGHT_W, h)

        ty = y_top - PAD - LH_RIGHT * 0.7
        for k, line in enumerate(title_lines):
            ax.text(
                RIGHT_X + 0.20,
                ty - k * LH_RIGHT,
                line,
                ha="left",
                va="center",
                fontsize=9.2,
                fontweight="bold",
                zorder=4,
            )
        offset = len(title_lines) + 1
        for k, reason in enumerate(reasons):
            ax.text(
                RIGHT_X + 0.20,
                ty - (offset + k) * LH_RIGHT,
                reason,
                ha="left",
                va="center",
                fontsize=8.9,
                zorder=4,
            )

        ay = Y(tops[arrow_row]) - main_hs[arrow_row] / 2
        h_arrow(ax, MAIN_X + MAIN_W, ay, RIGHT_X)

    # Stage sidebars
    for label, first_row, last_row in STAGES:
        y_top = Y(tops[first_row])
        y_bot = Y(tops[last_row]) - main_hs[last_row]
        bar_h = y_top - y_bot
        draw_rect(ax, BAR_X, y_bot, BAR_W, bar_h, face=SIDEBAR, lw=1.0)
        txt = ax.text(
            BAR_X + BAR_W / 2,
            (y_top + y_bot) / 2,
            label,
            ha="center",
            va="center",
            rotation=90,
            fontsize=10,
            fontweight="bold",
            zorder=4,
        )
        # Data units equal inches here, so shrink until the label clears the bar ends
        renderer = fig.canvas.get_renderer()
        avail = bar_h - 0.18
        drawn = txt.get_window_extent(renderer=renderer).height / fig.dpi
        if drawn > avail:
            txt.set_fontsize(max(5.5, 10.0 * avail / drawn))

    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "PRISMA_flow_diagram.pdf", bbox_inches="tight")
    fig.savefig(FIGS / "PRISMA_flow_diagram.png", dpi=300, bbox_inches="tight")
    fig.savefig(
        FIGS / "PRISMA_flow_diagram.tiff",
        dpi=300,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    return FIGS / "PRISMA_flow_diagram.pdf"


def main() -> None:
    write_counts_tsv()
    out = make_figure()
    print(f"Wrote {out}")
    print(f"Wrote {META / 'PRISMA_flow_counts.tsv'}")


if __name__ == "__main__":
    main()
