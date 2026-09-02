from pathlib import Path
import csv

MATCHES = Path(
    "work/codon/exact_species_final_one_to_one_min_cost_matches.tsv"
)

ORTHOLOGY = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

OUT = Path(
    "work/codon/exact_species_final_soluble_members.tsv"
)

# 1. Read final soluble orthogroups
soluble_groups = set()

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        soluble_groups.add(
            row["candidate_orthogroup"]
        )

print("Final soluble orthogroups:", len(soluble_groups))

# 2. Extract members from strict orthology table
rows = []

with ORTHOLOGY.open() as f:
    for line in f:
        line = line.rstrip("\n")

        if not line:
            continue

        orthogroup, full_id = line.split("\t")

        if orthogroup not in soluble_groups:
            continue

        if "|" not in full_id:
            raise ValueError(
                f"Unexpected full_id: {full_id}"
            )

        genome, protein_id = full_id.split("|", 1)

        rows.append(
            (
                orthogroup,
                genome,
                protein_id,
                full_id,
            )
        )

rows.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[2],
    )
)

groups_found = {
    row[0]
    for row in rows
}

missing_groups = (
    soluble_groups
    - groups_found
)

with OUT.open("w") as out:

    out.write(
        "orthogroup\t"
        "genome\t"
        "protein_id\t"
        "full_id\n"
    )

    for row in rows:
        out.write(
            "\t".join(row)
            + "\n"
        )

print("Groups extracted:", len(groups_found))
print("Protein members:", len(rows))
print("Missing groups:", len(missing_groups))

if missing_groups:
    print("Missing:")
    for g in sorted(missing_groups):
        print(g)

print("Output:", OUT)