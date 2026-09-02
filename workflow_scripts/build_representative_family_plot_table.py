from pathlib import Path
import csv

SEGMENTS = Path(
    "work/codon/segment_nearest_tmd.tsv"
)

OUT = Path(
    "work/codon/representative_family_plot_table.tsv"
)

# family -> TMD index selected for main figure
SELECTED = {
    "orthogroup_540": 1,
    "orthogroup_36": 2,
    "orthogroup_105": 8,
}

rows = []

with SEGMENTS.open() as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        group = row["orthogroup"]

        if group not in SELECTED:
            continue

        target_tmd = SELECTED[group]

        if int(row["nearest_tmd_index"]) != target_tmd:
            continue

        rows.append(
            (
                group,
                row["genome"],
                row["protein_id"],
                int(row["nearest_tmd_index"]),
                int(row["nearest_tmd_start"]),
                int(row["nearest_tmd_end"]),
                int(row["segment_start"]),
                int(row["segment_end"]),
                float(row["segment_center"]),
                float(row["distance_to_tmd_start"]),
                float(row["distance_to_tmd_end"]),
                row["relation"],
            )
        )

rows.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[6],
    )
)

with OUT.open("w") as out:

    out.write(
        "orthogroup\t"
        "genome\t"
        "protein_id\t"
        "tmd_index\t"
        "tmd_start\t"
        "tmd_end\t"
        "segment_start\t"
        "segment_end\t"
        "segment_center\t"
        "distance_to_tmd_start\t"
        "distance_to_tmd_end\t"
        "relation\n"
    )

    for row in rows:
        out.write(
            "\t".join(
                map(str, row)
            )
            + "\n"
        )

print("Rows written:", len(rows))

for group in SELECTED:

    group_rows = [
        r for r in rows
        if r[0] == group
    ]

    species = {
        r[1]
        for r in group_rows
    }

    print()
    print(group)
    print("  selected TMD:", SELECTED[group])
    print("  rows:", len(group_rows))
    print("  species with segment:", len(species))

print()
print("Output:", OUT)