from pathlib import Path
from collections import defaultdict
import csv
import statistics
import matplotlib.pyplot as plt


FINAL = Path("work/final_results")

START_NULL = FINAL / "synonymous_shuffle_null_start_summary_1000perm.tsv"
END_NULL = FINAL / "synonymous_shuffle_null_end_summary_1000perm.tsv"

WITHOUT_TMD1_UNITS = Path(
    "work/codon/synonymous_shuffle_null_units_1000perm.tsv"
)

NON_TMD_NULL = FINAL / "non_tmd_anchor_null_symmetric_1000perm.tsv"

OUT = FINAL / "figure3_global_conservation_v2.pdf"


def read_summary_medians(path):

    values = []

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            for key, text in row.items():

                if (
                    text not in (None, "")
                    and "median" in key.lower()
                ):

                    try:
                        values.append(float(text))
                        break
                    except ValueError:
                        pass

    return values


def build_without_tmd1_null(path):

    by_permutation = defaultdict(list)

    with path.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            tmd_index = int(
                row["tmd_index"]
            )

            if tmd_index == 1:
                continue

            species_total = int(
                row["species_total"]
            )

            species_with = int(
                row["species_with_feature"]
            )

            if (
                species_total < 8
                or species_with < 3
            ):
                continue

            variance = float(
                row["variance"]
            )

            permutation = int(
                row["permutation"]
            )

            by_permutation[
                permutation
            ].append(variance)

    medians = []

    for permutation in sorted(
        by_permutation
    ):

        values = by_permutation[
            permutation
        ]

        if values:

            medians.append(
                statistics.median(values)
            )

    return medians


def read_non_tmd_null(path):

    values = []

    with path.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            values.append(
                float(
                    row["median_variance"]
                )
            )

    return values


# ------------------------------------------------------------
# Load frozen null distributions
# ------------------------------------------------------------

start_null = read_summary_medians(
    START_NULL
)

end_null = read_summary_medians(
    END_NULL
)

without_tmd1_null = (
    build_without_tmd1_null(
        WITHOUT_TMD1_UNITS
    )
)

non_tmd_null = read_non_tmd_null(
    NON_TMD_NULL
)


# ------------------------------------------------------------
# Frozen observed statistics
# ------------------------------------------------------------

observed_start = 22.7
observed_end = 34.55

observed_without_tmd1 = 17.208333

observed_non_tmd_comparison = 47.2917


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, axes = plt.subplots(
    2,
    2,
    figsize=(10, 8)
)


def plot_null(
    ax,
    null_values,
    observed,
    title,
):

    ax.hist(
        null_values,
        bins=30
    )

    ax.axvline(
        observed,
        linestyle="--",
        linewidth=2
    )

    ax.set_title(title)

    ax.set_xlabel(
        "Null median variance (aa²)"
    )

    ax.set_ylabel(
        "Permutation count"
    )

    ax.text(
        0.98,
        0.95,
        f"Observed = {observed:.2f}\n"
        "p = 0.000999",
        transform=ax.transAxes,
        ha="right",
        va="top"
    )


plot_null(
    axes[0, 0],
    start_null,
    observed_start,
    "A  TMD start — synonymous shuffle",
)

plot_null(
    axes[0, 1],
    end_null,
    observed_end,
    "B  TMD end — synonymous shuffle",
)

plot_null(
    axes[1, 0],
    without_tmd1_null,
    observed_without_tmd1,
    "C  Excluding TMD1",
)

plot_null(
    axes[1, 1],
    non_tmd_null,
    observed_non_tmd_comparison,
    "D  Homologous non-TMD pseudo-anchors",
)


fig.suptitle(
    "Evolutionary positional constraint of low-adaptation segments relative to TMDs"
)

plt.tight_layout()

fig.savefig(
    OUT,
    bbox_inches="tight"
)


print(
    "START null permutations:",
    len(start_null)
)

print(
    "END null permutations:",
    len(end_null)
)

print(
    "Without-TMD1 null permutations:",
    len(without_tmd1_null)
)

print(
    "Non-TMD null permutations:",
    len(non_tmd_null)
)

print(
    "Without-TMD1 null median:",
    statistics.median(
        without_tmd1_null
    )
)

print(
    "Non-TMD null median:",
    statistics.median(
        non_tmd_null
    )
)

print(
    "Output:",
    OUT
)