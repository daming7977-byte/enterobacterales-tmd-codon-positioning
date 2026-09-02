from pathlib import Path
from collections import defaultdict, Counter
import csv

MATCHES = Path(
    "work/codon/soluble_control_top5_matches.tsv"
)

ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

# ------------------------------------------------------------
# Load all-GLOB groups
# ------------------------------------------------------------

with ALL_GLOB.open() as f:
    all_glob = {
        line.strip()
        for line in f
        if line.strip()
    }

# ------------------------------------------------------------
# Build membrane -> strict candidate list
# ------------------------------------------------------------

candidates = defaultdict(list)

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        soluble = row["candidate_orthogroup"]

        if soluble not in all_glob:
            continue

        membrane = row["membrane_orthogroup"]

        candidates[membrane].append(
            (
                int(row["candidate_rank"]),
                float(row["relative_length_difference"]),
                soluble,
            )
        )

for membrane in candidates:
    candidates[membrane].sort()

membranes = sorted(candidates)

print("Membrane groups with strict candidates:", len(membranes))

# ------------------------------------------------------------
# Naive best-match reuse
# ------------------------------------------------------------

naive = {
    m: candidates[m][0][2]
    for m in membranes
}

reuse = Counter(naive.values())

print("Unique soluble groups under naive best match:", len(reuse))
print("Maximum reuse of one soluble group:", max(reuse.values()))

print("Reuse distribution:")

reuse_dist = Counter(reuse.values())

for n in sorted(reuse_dist):
    print(
        f"  used by {n} membrane group(s): "
        f"{reuse_dist[n]} soluble groups"
    )

# ------------------------------------------------------------
# Maximum-cardinality one-to-one bipartite matching
# Deterministic order:
# candidate rank -> relative length difference -> ID
# ------------------------------------------------------------

soluble_to_membrane = {}

def try_match(membrane, seen):

    for rank, diff, soluble in candidates[membrane]:

        if soluble in seen:
            continue

        seen.add(soluble)

        if (
            soluble not in soluble_to_membrane
            or try_match(
                soluble_to_membrane[soluble],
                seen
            )
        ):
            soluble_to_membrane[soluble] = membrane
            return True

    return False


matched = 0

# Start with groups having fewer choices first
order = sorted(
    membranes,
    key=lambda m: (
        len(candidates[m]),
        candidates[m][0][0],
        candidates[m][0][1],
        m,
    )
)

for membrane in order:

    if try_match(membrane, set()):
        matched += 1

print()
print("Maximum one-to-one matches within strict top5:", matched)
print("Unmatched membrane groups:", len(membranes) - matched)

matched_membranes = set(soluble_to_membrane.values())

if len(membranes) - matched:
    print("Unmatched IDs:")
    for m in sorted(set(membranes) - matched_membranes):
        print(" ", m)