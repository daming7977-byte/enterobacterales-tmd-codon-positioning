from pathlib import Path
import csv
import matplotlib.pyplot as plt

INFILE = Path(
    "work/codon/representative_family_plot_table.tsv"
)

OUT = Path(
    "work/codon/FIGURE_2_representative_families.pdf"
)

GROUP_ORDER = [
    ("orthogroup_540", "TMD1"),
    ("orthogroup_36", "TMD2"),
    ("orthogroup_105", "TMD8"),
]

rows = []

with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

fig, axes = plt.subplots(
    1,
    3,
    figsize=(13, 5),
    sharey=False
)

for ax, (group, tmd_label) in zip(
    axes,
    GROUP_ORDER
):

    subset = [
        r for r in rows
        if r["orthogroup"] == group
    ]

    subset.sort(
        key=lambda r: r["genome"]
    )

    # Plot relative to TMD start = 0
    for y, row in enumerate(subset):

        tmd_start = int(row["tmd_start"])
        tmd_end = int(row["tmd_end"])

        seg_start = int(row["segment_start"])
        seg_end = int(row["segment_end"])

        rel_tmd_start = 0
        rel_tmd_end = tmd_end - tmd_start

        rel_seg_start = seg_start - tmd_start
        rel_seg_end = seg_end - tmd_start

        # baseline
        ax.plot(
            [-50, 30],
            [y, y],
            linewidth=0.8
        )

        # TMD
        ax.plot(
            [rel_tmd_start, rel_tmd_end],
            [y, y],
            linewidth=8,
            solid_capstyle="butt"
        )

        # low-adaptation segment
        ax.plot(
            [rel_seg_start, rel_seg_end],
            [y, y],
            linewidth=8,
            solid_capstyle="butt"
        )

    ax.axvline(
        0,
        linestyle="--",
        linewidth=1
    )

    ax.set_yticks(
        range(len(subset))
    )

    ax.set_yticklabels(
        [r["genome"] for r in subset],
        fontsize=8
    )

    ax.set_xlabel(
        "Position relative to TMD start (aa)"
    )

    ax.set_title(
        f"{group}\n{tmd_label}"
    )

    ax.set_xlim(-50, 30)

    ax.invert_yaxis()

axes[0].set_ylabel(
    "Species / genome"
)

fig.suptitle(
    "Representative conserved low-adaptation segments relative to TMD boundaries"
)

plt.tight_layout()

fig.savefig(
    OUT,
    bbox_inches="tight"
)

print("Output:", OUT)