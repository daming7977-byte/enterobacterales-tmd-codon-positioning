from pathlib import Path
from collections import defaultdict
import csv
import statistics

START_ASSIGN = Path(
    "work/codon/exact_species_final_soluble_segment_nearest_pseudo_start.tsv"
)

END_ASSIGN = Path(
    "work/codon/exact_species_final_soluble_segment_nearest_pseudo_end.tsv"
)


START_ANCHORS = Path(
    "work/codon/soluble_pseudo_anchor_start.tsv"
)

END_ANCHORS = Path(
    "work/codon/soluble_pseudo_anchor_end.tsv"
)

OUT_START = Path(
    "work/codon/exact_species_final_soluble_pseudo_start_cross_species_conservation.tsv"
)

OUT_END = Path(
    "work/codon/exact_species_final_soluble_pseudo_end_cross_species_conservation.tsv"
)

QC_THRESHOLD = 0.05


def build_conservation(
    assignment_file,
    anchor_file,
    outfile,
    label,
):

    # --------------------------------------------------------
    # Species totals from QC-passing anchors
    # --------------------------------------------------------

    species_total_by_unit = {}

    with anchor_file.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            diff = float(
                row["absolute_relative_difference"]
            )

            if diff > QC_THRESHOLD:
                continue

            key = (
                row["membrane_orthogroup"],
                int(row["tmd_index"]),
            )

            species_total_by_unit[key] = int(
                row["species_total"]
            )

    # --------------------------------------------------------
    # Keep one closest feature per species per unit
    # --------------------------------------------------------

    features = defaultdict(dict)

    with assignment_file.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        distance_col = (
            f"distance_to_pseudo_{label}"
        )

        for row in reader:

            key = (
                row["membrane_orthogroup"],
                int(row["pseudo_tmd_index"]),
            )

            genome = row["genome"]

            distance = float(
                row[distance_col]
            )

            previous = features[key].get(
                genome
            )

            if (
                previous is None
                or abs(distance) < abs(previous)
            ):
                features[key][genome] = distance

    # --------------------------------------------------------
    # Build table
    # --------------------------------------------------------

    rows_out = []

    for key in sorted(
        species_total_by_unit
    ):

        membrane_group, tmd_index = key

        species_total = (
            species_total_by_unit[key]
        )

        vals = list(
            features.get(
                key,
                {}
            ).values()
        )

        species_with_feature = len(vals)

        presence_fraction = (
            species_with_feature
            / species_total
            if species_total
            else 0
        )

        if len(vals) >= 2:

            median_distance = (
                statistics.median(vals)
            )

            mean_distance = (
                statistics.mean(vals)
            )

            sd = statistics.stdev(vals)

            variance = (
                statistics.variance(vals)
            )

        elif len(vals) == 1:

            median_distance = vals[0]
            mean_distance = vals[0]
            sd = ""
            variance = ""

        else:

            median_distance = ""
            mean_distance = ""
            sd = ""
            variance = ""

        rows_out.append(
            (
                membrane_group,
                tmd_index,
                species_total,
                species_with_feature,
                presence_fraction,
                median_distance,
                mean_distance,
                sd,
                variance,
            )
        )

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    with outfile.open("w") as out:

        out.write(
            "membrane_orthogroup\t"
            "pseudo_tmd_index\t"
            "species_total\t"
            "species_with_feature\t"
            "presence_fraction\t"
            f"median_distance_{label}\t"
            f"mean_distance_{label}\t"
            f"sd_distance_{label}\t"
            f"variance_distance_{label}\n"
        )

        for row in rows_out:

            formatted = []

            for x in row:

                if isinstance(x, float):
                    formatted.append(
                        f"{x:.6f}"
                    )
                else:
                    formatted.append(
                        str(x)
                    )

            out.write(
                "\t".join(formatted)
                + "\n"
            )

    # --------------------------------------------------------
    # Primary qualifying units
    # --------------------------------------------------------

    qualifying = [
        r
        for r in rows_out
        if r[2] >= 8
        and r[3] >= 3
        and r[8] != ""
    ]

    variances = [
        float(r[8])
        for r in qualifying
    ]

    print()
    print(
        f"=== SOLUBLE {label.upper()} ==="
    )

    print(
        "Pseudo-anchor units:",
        len(rows_out)
    )

    print(
        "Qualifying units:",
        len(qualifying)
    )

    if variances:

        print(
            "Median variance:",
            statistics.median(
                variances
            )
        )

        print(
            "Mean variance:",
            statistics.mean(
                variances
            )
        )

        print(
            "Min variance:",
            min(variances)
        )

        print(
            "Max variance:",
            max(variances)
        )

    print(
        "Output:",
        outfile
    )


build_conservation(
    START_ASSIGN,
    START_ANCHORS,
    OUT_START,
    "start",
)

build_conservation(
    END_ASSIGN,
    END_ANCHORS,
    OUT_END,
    "end",
)