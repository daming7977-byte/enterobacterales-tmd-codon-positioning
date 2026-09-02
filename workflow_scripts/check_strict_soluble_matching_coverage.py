from pathlib import Path
from collections import defaultdict
import csv

MATCHES = Path(
    "work/codon/soluble_control_top5_matches.tsv"
)

ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

OUT_BEST = Path(
    "work/codon/strict_soluble_best_matches.tsv"
)

OUT_SUMMARY = Path(
    "work/codon/strict_soluble_matching_summary.tsv"
)

# ------------------------------------------------------------
# Load strict all-GLOB groups
# ------------------------------------------------------------

all_glob = set()

with ALL_GLOB.open() as f:
    for line in f:
        g = line.strip()
        if g:
            all_glob.add(g)

print("All-GLOB candidate orthogroups:", len(all_glob))

# ------------------------------------------------------------
# Read top-5 length-matched candidates
# ------------------------------------------------------------

by_membrane = defaultdict(list)

with MATCHES.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        membrane = row["membrane_orthogroup"]
        candidate = row["candidate_orthogroup"]

        row["candidate_rank"] = int(row["candidate_rank"])
        row["species_coverage"] = int(row["species_coverage"])
        row["membrane_median_length"] = float(
            row["membrane_median_length"]
        )
        row["candidate_median_length"] = float(
            row["candidate_median_length"]
        )
        row["absolute_length_difference"] = float(
            row["absolute_length_difference"]
        )
        row["relative_length_difference"] = float(
            row["relative_length_difference"]
        )

        if candidate in all_glob:
            by_membrane[membrane].append(row)

# ------------------------------------------------------------
# Count coverage and choose deterministic best match
# ------------------------------------------------------------

best_rows = []
candidate_counts = {}

for membrane, rows in by_membrane.items():

    rows = sorted(
        rows,
        key=lambda r: (
            r["candidate_rank"],
            r["relative_length_difference"],
            r["candidate_orthogroup"],
        )
    )

    candidate_counts[membrane] = len(rows)

    best_rows.append(rows[0])

# Number of membrane groups represented in original top5 file
all_membrane = set()

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        all_membrane.add(row["membrane_orthogroup"])

print("Membrane groups represented in top5 table:", len(all_membrane))
print(
    "Membrane groups with >=1 strict all-GLOB top5 match:",
    len(best_rows)
)

if all_membrane:
    print(
        "Coverage among originally matchable groups:",
        len(best_rows) / len(all_membrane)
    )

# ------------------------------------------------------------
# Candidate-count distribution
# ------------------------------------------------------------

dist = defaultdict(int)

for membrane in all_membrane:
    n = candidate_counts.get(membrane, 0)
    dist[n] += 1

print("Strict candidate counts per membrane group:")

for n in sorted(dist):
    print(f"  {n} candidates: {dist[n]} groups")

# ------------------------------------------------------------
# Write deterministic best-match table
# ------------------------------------------------------------

fieldnames = [
    "membrane_orthogroup",
    "species_coverage",
    "membrane_median_length",
    "candidate_rank",
    "candidate_orthogroup",
    "candidate_median_length",
    "absolute_length_difference",
    "relative_length_difference",
]

with OUT_BEST.open("w") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for row in sorted(
        best_rows,
        key=lambda r: r["membrane_orthogroup"]
    ):
        writer.writerow(
            {k: row[k] for k in fieldnames}
        )

# ------------------------------------------------------------
# Write summary for every membrane group
# ------------------------------------------------------------

with OUT_SUMMARY.open("w") as out:

    out.write(
        "membrane_orthogroup\t"
        "strict_all_GLOB_candidates_in_top5\t"
        "has_strict_match\n"
    )

    for membrane in sorted(all_membrane):

        n = candidate_counts.get(membrane, 0)

        out.write(
            f"{membrane}\t{n}\t{n > 0}\n"
        )

print("Best-match table:", OUT_BEST)
print("Summary table:", OUT_SUMMARY)