from pathlib import Path
from collections import defaultdict
import csv
import statistics

INPUT = Path(
    "work/codon/modal_tmd_start_cross_species_conservation.tsv"
)

OUT = Path(
    "work/codon/leave_one_orthogroup_out_start.tsv"
)

# ------------------------------------------------------------
# Load qualifying start-anchor units
# ------------------------------------------------------------

units = []

with INPUT.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        species_total = int(row["species_total"])
        species_with_feature = int(row["species_with_feature"])
        var = row["variance_distance_start"]

        if (
            species_total >= 8
            and species_with_feature >= 3
            and var != ""
        ):
            units.append(
                (
                    row["orthogroup"],
                    int(row["tmd_index"]),
                    float(var),
                )
            )

# ------------------------------------------------------------
# Full-data median
# ------------------------------------------------------------

all_variances = [x[2] for x in units]
full_median = statistics.median(all_variances)

groups = sorted({
    x[0]
    for x in units
})

# ------------------------------------------------------------
# Leave one orthogroup out
# ------------------------------------------------------------

rows_out = []

for group in groups:

    retained = [
        x[2]
        for x in units
        if x[0] != group
    ]

    removed_n = sum(
        1
        for x in units
        if x[0] == group
    )

    loo_median = statistics.median(retained)

    change = loo_median - full_median

    rows_out.append(
        (
            group,
            removed_n,
            len(retained),
            loo_median,
            change,
        )
    )

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

with OUT.open("w") as out:

    out.write(
        "removed_orthogroup\t"
        "removed_units\t"
        "retained_units\t"
        "loo_median_variance\t"
        "change_from_full_median\n"
    )

    for row in rows_out:
        out.write(
            "\t".join(map(str, row)) + "\n"
        )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

loo_values = [
    row[3]
    for row in rows_out
]

changes = [
    row[4]
    for row in rows_out
]

print("Qualifying units:", len(units))
print("Orthogroups:", len(groups))
print("Full median variance:", full_median)

print()
print("Leave-one-out minimum median:", min(loo_values))
print("Leave-one-out maximum median:", max(loo_values))
print("Largest absolute median change:", max(abs(x) for x in changes))

print()
print("Most influential removals:")

for row in sorted(
    rows_out,
    key=lambda x: abs(x[4]),
    reverse=True
)[:15]:

    print(
        row[0],
        "removed_units=",
        row[1],
        "loo_median=",
        row[3],
        "change=",
        row[4],
    )

print()
print("Output:", OUT)