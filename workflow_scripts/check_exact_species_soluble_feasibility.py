from pathlib import Path
from collections import defaultdict
import csv

MATCHES = Path(
    "work/codon/soluble_control_top5_matches.tsv"
)

ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

# ------------------------------------------------------------
# Species sets
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
# Strict all-GLOB groups
# ------------------------------------------------------------

with ALL_GLOB.open() as f:
    all_glob = {
        x.strip()
        for x in f
        if x.strip()
    }

# ------------------------------------------------------------
# Exact-species candidates inside existing top5
# ------------------------------------------------------------

candidates = defaultdict(list)
all_membranes = set()

with MATCHES.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        membrane = row["membrane_orthogroup"]
        soluble = row["candidate_orthogroup"]

        all_membranes.add(membrane)

        if soluble not in all_glob:
            continue

        if species_sets[membrane] != species_sets[soluble]:
            continue

        candidates[membrane].append(
            (
                int(row["candidate_rank"]),
                float(row["relative_length_difference"]),
                soluble,
            )
        )

for m in candidates:
    candidates[m].sort()

print("Original membrane groups in top5 table:", len(all_membranes))
print(
    "Membrane groups with >=1 exact-species all-GLOB candidate:",
    len(candidates)
)

print(
    "Fraction:",
    len(candidates) / len(all_membranes)
)

# Candidate number distribution
dist = defaultdict(int)

for m in all_membranes:
    dist[len(candidates.get(m, []))] += 1

print("Exact-species strict candidate counts:")

for n in sorted(dist):
    print(
        f"  {n} candidates: {dist[n]} groups"
    )

# ------------------------------------------------------------
# Maximum-cardinality one-to-one matching
# ------------------------------------------------------------

soluble_to_membrane = {}

def try_match(membrane, seen):

    for rank, diff, soluble in candidates.get(membrane, []):

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


order = sorted(
    candidates,
    key=lambda m: (
        len(candidates[m]),
        candidates[m][0][1],
        m,
    )
)

matched = 0

for m in order:

    if try_match(m, set()):
        matched += 1

print()
print(
    "Maximum one-to-one exact-species matches within current top5:",
    matched
)

print(
    "Fraction of original 222:",
    matched / len(all_membranes)
)