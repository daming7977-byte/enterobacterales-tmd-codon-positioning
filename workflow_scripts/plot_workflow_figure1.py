from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


OUT = Path(
    "work/final_results/FIGURE_1_workflow.pdf"
)


fig, ax = plt.subplots(
    figsize=(15, 7)
)

ax.set_xlim(0, 15)
ax.set_ylim(0, 7)
ax.axis("off")


def add_box(
    x,
    y,
    w,
    h,
    title,
    text,
):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.03",
        linewidth=1.5,
        fill=False
    )

    ax.add_patch(box)

    ax.text(
        x + w / 2,
        y + h - 0.28,
        title,
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold"
    )

    ax.text(
        x + w / 2,
        y + h / 2 - 0.1,
        text,
        ha="center",
        va="center",
        fontsize=10,
        wrap=True
    )


def add_arrow(x1, y1, x2, y2):

    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="->",
        mutation_scale=15,
        linewidth=1.5
    )

    ax.add_patch(arrow)


# ------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------

add_box(
    0.4, 4.2, 2.2, 1.7,
    "1. Species panel",
    "10 Enterobacterales\nreference genomes"
)

add_box(
    3.0, 4.2, 2.2, 1.7,
    "2. Strict orthology",
    "1,836 one-to-one groups\npresent in ≥8 species"
)

add_box(
    5.6, 4.2, 2.2, 1.7,
    "3. Membrane QC",
    "Reviewed topology +\nDeepTMHMM completion\n225 final families"
)

add_box(
    8.2, 4.2, 2.2, 1.7,
    "4. Codon metric",
    "Species-specific\nsynonymous codon\nadaptation weights"
)

add_box(
    10.8, 4.2, 2.2, 1.7,
    "5. Low-adaptation\nsegments",
    "Bottom-decile codons\n≥3 consecutive codons"
)

add_box(
    13.4, 4.2, 1.2, 1.7,
    "6. TMD\nanchors",
    "Start\nEnd\nCenter"
)


add_arrow(2.6, 5.05, 3.0, 5.05)
add_arrow(5.2, 5.05, 5.6, 5.05)
add_arrow(7.8, 5.05, 8.2, 5.05)
add_arrow(10.4, 5.05, 10.8, 5.05)
add_arrow(13.0, 5.05, 13.4, 5.05)


# ------------------------------------------------------------
# Analysis layer
# ------------------------------------------------------------

add_box(
    4.0, 1.5, 3.0, 1.6,
    "Cross-species conservation",
    "Variance of low-adaptation\nsegment position relative to\nhomologous TMD boundaries"
)

add_box(
    7.6, 1.5, 3.0, 1.6,
    "Primary null model",
    "1,000 synonymous-\ncomposition-preserving\npermutations"
)

add_box(
    11.2, 1.5, 3.0, 1.6,
    "Robustness / specificity",
    "Exclude TMD1\nFixed-presence null\nNon-TMD pseudo-anchors\nMatched soluble families"
)


add_arrow(
    14.0, 4.2,
    6.6, 3.1
)

add_arrow(
    7.0, 2.3,
    7.6, 2.3
)

add_arrow(
    10.6, 2.3,
    11.2, 2.3
)


# ------------------------------------------------------------
# Frozen study numbers
# ------------------------------------------------------------

ax.text(
    7.5,
    6.55,
    "Study design for evolutionary analysis of TMD-relative synonymous codon positioning",
    ha="center",
    va="center",
    fontsize=16,
    fontweight="bold"
)

ax.text(
    7.5,
    0.55,
    "Analysis frozen 2026-08-31  |  Primary qualifying rule: species_total ≥8 and species_with_feature ≥3",
    ha="center",
    va="center",
    fontsize=10
)


plt.tight_layout()

fig.savefig(
    OUT,
    bbox_inches="tight"
)

print("Output:", OUT)