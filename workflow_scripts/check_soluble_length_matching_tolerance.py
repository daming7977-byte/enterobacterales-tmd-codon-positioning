from pathlib import Path
import csv

INPUT = Path("work/codon/soluble_control_matching_space.tsv")

TOLERANCES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

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

print("Membrane groups:", len(membrane))
print("Candidate groups:", len(candidates))
print()
print("tolerance\tmatchable\tfraction\tmedian_candidates_per_membrane")

for tol in TOLERANCES:

    counts = []

    for m in membrane:

        valid = [
            c for c in candidates
            if c["coverage"] == m["coverage"]
            and abs(c["length"] - m["length"]) / m["length"] <= tol
        ]

        counts.append(len(valid))

    matchable = sum(n > 0 for n in counts)

    nonzero = sorted(n for n in counts if n > 0)

    if nonzero:
        median_n = nonzero[len(nonzero)//2]
    else:
        median_n = 0

    print(
        f"{tol:.2f}\t"
        f"{matchable}\t"
        f"{matchable/len(membrane):.4f}\t"
        f"{median_n}"
    )