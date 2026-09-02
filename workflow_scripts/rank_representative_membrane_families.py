from pathlib import Path
from collections import defaultdict
import csv

START = Path(
    "work/codon/modal_tmd_start_cross_species_conservation.tsv"
)

END = Path(
    "work/codon/modal_tmd_end_cross_species_conservation.tsv"
)

OUT = Path(
    "work/codon/representative_membrane_family_ranking.tsv"
)

groups = defaultdict(
    lambda: {
        "start_units": [],
        "end_units": [],
        "modal_tmd_count": None,
    }
)


def load_table(path, label):

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            group = row["orthogroup"]

            modal = int(
                row["modal_tmd_count"]
            )

            groups[group][
                "modal_tmd_count"
            ] = modal

            species_total = int(
                row["species_total"]
            )

            species_with = int(
                row["species_with_feature"]
            )

            var_text = row[
                f"variance_distance_{label}"
            ]

            if (
                species_total >= 8
                and species_with >= 3
                and var_text != ""
            ):

                groups[group][
                    f"{label}_units"
                ].append(
                    (
                        int(row["tmd_index"]),
                        species_total,
                        species_with,
                        float(var_text),
                    )
                )


load_table(
    START,
    "start"
)

load_table(
    END,
    "end"
)


rows = []

for group, info in groups.items():

    start_units = info["start_units"]
    end_units = info["end_units"]

    if not start_units and not end_units:
        continue

    modal = info["modal_tmd_count"]

    all_vars = (
        [x[3] for x in start_units]
        + [x[3] for x in end_units]
    )

    best_var = min(all_vars)

    max_species_with = max(
        [x[2] for x in start_units + end_units]
    )

    n_qualifying = (
        len(start_units)
        + len(end_units)
    )

    score = (
        best_var
        + 5 * max(0, modal - 8)
        - 2 * n_qualifying
        - 2 * max_species_with
    )

    rows.append(
        (
            group,
            modal,
            len(start_units),
            len(end_units),
            n_qualifying,
            max_species_with,
            best_var,
            score,
        )
    )


rows.sort(
    key=lambda x: (
        x[7],
        x[6],
        x[0],
    )
)


with OUT.open("w") as out:

    out.write(
        "orthogroup\t"
        "modal_tmd_count\t"
        "start_qualifying_units\t"
        "end_qualifying_units\t"
        "total_qualifying_units\t"
        "max_species_with_feature\t"
        "best_variance\t"
        "ranking_score\n"
    )

    for row in rows:

        out.write(
            "\t".join(
                map(str, row)
            )
            + "\n"
        )


print(
    "Families ranked:",
    len(rows)
)

print()
print("Top 20:")

for row in rows[:20]:
    print(
        row[0],
        "TMDs=", row[1],
        "start=", row[2],
        "end=", row[3],
        "units=", row[4],
        "max_species=", row[5],
        "best_var=", row[6],
        "score=", row[7],
    )

print()
print("Output:", OUT)