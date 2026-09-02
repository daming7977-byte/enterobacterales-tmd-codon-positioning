from pathlib import Path
from collections import defaultdict
import csv

MATCHES = Path(
    "work/codon/soluble_control_top5_matches.tsv"
)

ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

OUT = Path(
    "work/codon/strict_soluble_one_to_one_matches.tsv"
)

# ------------------------------------------------------------
# Load all-GLOB candidate groups
# ------------------------------------------------------------

with ALL_GLOB.open() as f:
    all_glob = {
        line.strip()
        for line in f
        if line.strip()
    }

# ------------------------------------------------------------
# Build candidate lists
# ------------------------------------------------------------

candidates = defaultdict(list)

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        soluble = row["candidate_orthogroup"]

        if soluble not in all_glob:
            continue

        membrane = row["membrane_orthogroup"]

        candidates[membrane].append({
            "membrane_orthogroup": membrane,
            "species_coverage": int(row["species_coverage"]),
            "membrane_median_length": float(
                row["membrane_median_length"]
            ),
            "candidate_rank": int(row["candidate_rank"]),
            "candidate_orthogroup": soluble,
            "candidate_median_length": float(
                row["candidate_median_length"]
            ),
            "absolute_length_difference": float(
                row["absolute_length_difference"]
            ),
            "relative_length_difference": float(
                row["relative_length_difference"]
            ),
        })

for membrane in candidates:
    candidates[membrane].sort(
        key=lambda r: (
            r["relative_length_difference"],
            r["candidate_rank"],
            r["candidate_orthogroup"],
        )
    )

# ------------------------------------------------------------
# Deterministic augmenting-path matching
#
# Process the most constrained membrane groups first.
# Candidate order favors smaller relative length difference.
# ------------------------------------------------------------

soluble_to_membrane = {}
chosen_row = {}

order = sorted(
    candidates,
    key=lambda m: (
        len(candidates[m]),
        candidates[m][0]["relative_length_difference"],
        m,
    )
)

def try_match(membrane, seen_soluble, seen_membrane):

    if membrane in seen_membrane:
        return False

    seen_membrane.add(membrane)

    for row in candidates[membrane]:

        soluble = row["candidate_orthogroup"]

        if soluble in seen_soluble:
            continue

        seen_soluble.add(soluble)

        if soluble not in soluble_to_membrane:

            soluble_to_membrane[soluble] = membrane
            chosen_row[membrane] = row
            return True

        old_membrane = soluble_to_membrane[soluble]

        if try_match(
            old_membrane,
            seen_soluble,
            seen_membrane,
        ):
            soluble_to_membrane[soluble] = membrane
            chosen_row[membrane] = row
            return True

    return False


for membrane in order:
    try_match(membrane, set(), set())

# Reconstruct exact chosen row for each final pair
final_rows = []

for soluble, membrane in soluble_to_membrane.items():

    matches = [
        row
        for row in candidates[membrane]
        if row["candidate_orthogroup"] == soluble
    ]

    assert len(matches) == 1

    final_rows.append(matches[0])

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

final_rows.sort(
    key=lambda r: r["membrane_orthogroup"]
)

membrane_ids = {
    r["membrane_orthogroup"]
    for r in final_rows
}

soluble_ids = {
    r["candidate_orthogroup"]
    for r in final_rows
}

all_membranes = set(candidates)

unmatched = sorted(
    all_membranes - membrane_ids
)

print("Final one-to-one matches:", len(final_rows))
print("Unique membrane groups:", len(membrane_ids))
print("Unique soluble groups:", len(soluble_ids))
print("Unmatched among 220 candidate-bearing groups:", len(unmatched))

if unmatched:
    print("Unmatched IDs:")
    for x in unmatched:
        print(" ", x)

# Rank distribution
rank_counts = defaultdict(int)

for row in final_rows:
    rank_counts[row["candidate_rank"]] += 1

print("Candidate rank distribution:")

for rank in sorted(rank_counts):
    print(
        f"  rank {rank}: {rank_counts[rank]}"
    )

# Length-difference summary
diffs = sorted(
    r["relative_length_difference"]
    for r in final_rows
)

def median(values):
    n = len(values)
    if n % 2:
        return values[n // 2]
    return (
        values[n // 2 - 1]
        + values[n // 2]
    ) / 2

print(
    "Median relative length difference:",
    median(diffs)
)

print(
    "Maximum relative length difference:",
    max(diffs)
)

print(
    "Mean relative length difference:",
    sum(diffs) / len(diffs)
)

# ------------------------------------------------------------
# Write final table
# ------------------------------------------------------------

fields = [
    "membrane_orthogroup",
    "species_coverage",
    "membrane_median_length",
    "candidate_rank",
    "candidate_orthogroup",
    "candidate_median_length",
    "absolute_length_difference",
    "relative_length_difference",
]

with OUT.open("w") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for row in final_rows:

        writer.writerow(
            {
                field: row[field]
                for field in fields
            }
        )

print("Output:", OUT)