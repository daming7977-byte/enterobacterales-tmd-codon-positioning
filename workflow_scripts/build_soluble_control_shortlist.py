from pathlib import Path
import csv
from collections import defaultdict

INPUT = Path(
    "work/codon/soluble_control_matching_space.tsv"
)

OUT_MATCHES = Path(
    "work/codon/soluble_control_top5_matches.tsv"
)

OUT_GROUPS = Path(
    "work/codon/soluble_control_shortlist_groups.txt"
)

TOLERANCE = 0.05
TOP_N = 5

# ------------------------------------------------------------
# Load membrane and candidate groups
# ------------------------------------------------------------

membrane = []
candidates = []

with INPUT.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        item = {
            "group": row["orthogroup"],
            "coverage": int(row["species_coverage"]),
            "length": float(row["median_protein_length"]),
        }

        if row["membrane_group"] == "True":
            membrane.append(item)
        else:
            candidates.append(item)

# ------------------------------------------------------------
# Build top-N closest matches
# ------------------------------------------------------------

rows_out = []
unique_shortlist = set()
unmatched = []

for m in membrane:

    valid = []

    for c in candidates:

        if c["coverage"] != m["coverage"]:
            continue

        relative_difference = (
            abs(c["length"] - m["length"])
            / m["length"]
        )

        if relative_difference <= TOLERANCE:

            valid.append(
                (
                    relative_difference,
                    abs(c["length"] - m["length"]),
                    c["group"],
                    c["length"],
                )
            )

    valid.sort(
        key=lambda x: (
            x[0],
            x[1],
            x[2],
        )
    )

    selected = valid[:TOP_N]

    if not selected:
        unmatched.append(m["group"])
        continue

    for rank, item in enumerate(selected, start=1):

        (
            relative_difference,
            absolute_difference,
            candidate_group,
            candidate_length,
        ) = item

        rows_out.append(
            (
                m["group"],
                m["coverage"],
                m["length"],
                rank,
                candidate_group,
                candidate_length,
                absolute_difference,
                relative_difference,
            )
        )

        unique_shortlist.add(candidate_group)

# ------------------------------------------------------------
# Write match table
# ------------------------------------------------------------

with OUT_MATCHES.open("w") as out:

    out.write(
        "membrane_orthogroup\t"
        "species_coverage\t"
        "membrane_median_length\t"
        "candidate_rank\t"
        "candidate_orthogroup\t"
        "candidate_median_length\t"
        "absolute_length_difference\t"
        "relative_length_difference\n"
    )

    for row in rows_out:

        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")

# ------------------------------------------------------------
# Write unique shortlist
# ------------------------------------------------------------

with OUT_GROUPS.open("w") as out:

    for group in sorted(unique_shortlist):
        out.write(group + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

per_membrane = defaultdict(int)

for row in rows_out:
    per_membrane[row[0]] += 1

print("Membrane groups:", len(membrane))
print("Membrane groups with >=1 shortlist candidate:", len(per_membrane))
print("Total membrane-candidate pairs:", len(rows_out))
print("Unique candidate orthogroups:", len(unique_shortlist))
print("Unmatched membrane groups:", len(unmatched))

print()
print("Unmatched:")
for group in unmatched:
    print(group)

print()
print("Outputs:")
print(OUT_MATCHES)
print(OUT_GROUPS)