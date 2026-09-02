from pathlib import Path
from collections import defaultdict
import csv

MATCHING_SPACE = Path(
    "work/codon/soluble_control_matching_space.tsv"
)

MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

ALREADY_SCREENED = Path(
    "work/codon/soluble_control_shortlist_groups.txt"
)

OUT_PAIRS = Path(
    "work/codon/exact_species_full_soluble_candidates.tsv"
)

OUT_GROUPS = Path(
    "work/codon/exact_species_full_soluble_candidate_groups.txt"
)

OUT_NEW_GROUPS = Path(
    "work/codon/exact_species_new_unscreened_groups.txt"
)

TOLERANCE = 0.05

# ------------------------------------------------------------
# Load exact species sets
# ------------------------------------------------------------

species_sets = defaultdict(set)

with MEMBERS.open() as f:
    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")
        genome = full_id.split("|", 1)[0]

        species_sets[group].add(genome)

# ------------------------------------------------------------
# Load group lengths / membrane status
# ------------------------------------------------------------

membrane = []
candidates = []

with MATCHING_SPACE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        item = {
            "group": row["orthogroup"],
            "length": float(row["median_protein_length"]),
        }

        if row["membrane_group"] == "True":
            membrane.append(item)
        else:
            candidates.append(item)

# ------------------------------------------------------------
# Build FULL exact-species candidate space
# ------------------------------------------------------------

rows = []
candidate_groups = set()
per_membrane = defaultdict(int)

for m in membrane:

    mset = species_sets[m["group"]]

    valid = []

    for c in candidates:

        cset = species_sets[c["group"]]

        if cset != mset:
            continue

        relative_difference = (
            abs(c["length"] - m["length"])
            / m["length"]
        )

        if relative_difference > TOLERANCE:
            continue

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

    for rank, item in enumerate(valid, start=1):

        (
            relative_difference,
            absolute_difference,
            candidate_group,
            candidate_length,
        ) = item

        rows.append(
            (
                m["group"],
                len(mset),
                m["length"],
                rank,
                candidate_group,
                candidate_length,
                absolute_difference,
                relative_difference,
            )
        )

        candidate_groups.add(candidate_group)
        per_membrane[m["group"]] += 1

# ------------------------------------------------------------
# Already DeepTMHMM-screened groups
# ------------------------------------------------------------

with ALREADY_SCREENED.open() as f:
    already_screened = {
        x.strip()
        for x in f
        if x.strip()
    }

new_groups = candidate_groups - already_screened

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

with OUT_PAIRS.open("w") as out:

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

    for row in rows:

        out.write(
            "\t".join(
                [
                    str(row[0]),
                    str(row[1]),
                    f"{row[2]:.6f}",
                    str(row[3]),
                    str(row[4]),
                    f"{row[5]:.6f}",
                    f"{row[6]:.6f}",
                    f"{row[7]:.6f}",
                ]
            )
            + "\n"
        )

with OUT_GROUPS.open("w") as out:
    for g in sorted(candidate_groups):
        out.write(g + "\n")

with OUT_NEW_GROUPS.open("w") as out:
    for g in sorted(new_groups):
        out.write(g + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Membrane groups:", len(membrane))
print(
    "Membrane groups with >=1 exact-species ±5% candidate:",
    len(per_membrane)
)

print(
    "Coverage fraction:",
    len(per_membrane) / len(membrane)
)

print(
    "Total membrane-candidate pairs:",
    len(rows)
)

print(
    "Unique exact-species candidate orthogroups:",
    len(candidate_groups)
)

print(
    "Already DeepTMHMM-screened candidate groups:",
    len(candidate_groups & already_screened)
)

print(
    "NEW candidate groups needing DeepTMHMM:",
    len(new_groups)
)

unmatched = sorted(
    {m["group"] for m in membrane}
    - set(per_membrane)
)

print("Membrane groups with no exact-species candidate:", len(unmatched))

if unmatched:
    for g in unmatched:
        print(" ", g)

print()
print("Outputs:")
print(OUT_PAIRS)
print(OUT_GROUPS)
print(OUT_NEW_GROUPS)