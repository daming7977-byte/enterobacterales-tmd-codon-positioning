from pathlib import Path
import csv
import matplotlib.pyplot as plt


INFILE = Path(
    "work/final_results/membrane_vs_soluble_family_paired.tsv"
)

OUT = Path(
    "work/final_results/FIGURE_4_matched_soluble.pdf"
)


rows = {
    "start": [],
    "end": [],
}

with INFILE.open() as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        label = row["anchor_type"]

        rows[label].append(
            (
                row["membrane_orthogroup"],
                float(
                    row[
                        "membrane_family_median_variance"
                    ]
                ),
                float(
                    row[
                        "soluble_family_median_variance"
                    ]
                ),
            )
        )


fig, axes = plt.subplots(
    1,
    2,
    figsize=(9, 5)
)


for ax, label in zip(
    axes,
    ["start", "end"]
):

    data = rows[label]

    # Draw each matched family pair
    for i, (group, membrane, soluble) in enumerate(data):

        ax.plot(
            [0, 1],
            [membrane, soluble],
            marker="o",
            linewidth=1,
            alpha=0.7
        )

    ax.set_yscale("symlog", linthresh=1)

    ax.set_xticks(
        [0, 1]
    )

    ax.set_xticks(
        [0, 1]
    )

    ax.set_xticklabels(
        ["Membrane", "Matched soluble"]
    )

    ax.set_ylabel(
    "Family-level median positional variance (aa²; symlog scale)"
    )

    if label == "start":

        ax.set_title(
            "A  TMD start"
        )

        ax.text(
            0.5,
            0.97,
            "17 matched families\n"
            "12/17 soluble > membrane\n"
            "sign-flip P = 0.0102\n"
            "sign test P = 0.0717",
            transform=ax.transAxes,
            ha="center",
            va="top"
        )

    else:

        ax.set_title(
            "B  TMD end"
        )

        ax.text(
            0.5,
            0.97,
            "17 matched families\n"
            "10/17 soluble > membrane\n"
            "sign-flip P = 0.0738\n"
            "sign test P = 0.3145",
            transform=ax.transAxes,
            ha="center",
            va="top"
        )


fig.suptitle(
    "Matched soluble-family comparison of positional constraint"
)

plt.tight_layout()

fig.savefig(
    OUT,
    bbox_inches="tight"
)

print(
    "START matched families:",
    len(rows["start"])
)

print(
    "END matched families:",
    len(rows["end"])
)

print(
    "Output:",
    OUT
)