from pathlib import Path
from collections import defaultdict
import csv

MATCHES = Path(
    "work/codon/strict_soluble_one_to_one_min_cost_matches.tsv"
)

MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

# ------------------------------------------------------------
# Load species sets for every orthogroup
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
# Compare final matched pairs
# ------------------------------------------------------------

total = 0
exact = 0

mismatches = []

with MATCHES.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        membrane = row["membrane_orthogroup"]
        soluble = row["candidate_orthogroup"]

        mset = species_sets[membrane]
        sset = species_sets[soluble]

        total += 1

        if mset == sset:
            exact += 1

        else:
            mismatches.append(
                (
                    membrane,
                    soluble,
                    sorted(mset - sset),
                    sorted(sset - mset),
                )
            )

print("Final matched pairs:", total)
print("Exact species-set matches:", exact)
print("Species-set mismatches:", len(mismatches))

if total:
    print(
        "Exact species-set fraction:",
        exact / total
    )

if mismatches:

    print()
    print("First 20 mismatches:")

    for membrane, soluble, missing_in_soluble, extra_in_soluble in mismatches[:20]:

        print(
            membrane,
            soluble,
            "missing_in_soluble=",
            ",".join(missing_in_soluble),
            "extra_in_soluble=",
            ",".join(extra_in_soluble),
        )