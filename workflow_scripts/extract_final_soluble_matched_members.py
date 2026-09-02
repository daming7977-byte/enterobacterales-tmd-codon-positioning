from pathlib import Path
import csv

MATCHES = Path(
    "work/codon/strict_soluble_one_to_one_min_cost_matches.tsv"
)

ALL_MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

PROTEIN_FASTA = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

OUT_MEMBERS = Path(
    "work/codon/final_soluble_matched_members.tsv"
)

OUT_FASTA = Path(
    "work/codon/final_soluble_matched_all.faa"
)

# ------------------------------------------------------------
# Load final 217 soluble orthogroups
# ------------------------------------------------------------

selected_groups = set()

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        selected_groups.add(
            row["candidate_orthogroup"]
        )

print(
    "Selected final soluble orthogroups:",
    len(selected_groups)
)

# ------------------------------------------------------------
# Extract members
# ------------------------------------------------------------

members = []

with ALL_MEMBERS.open() as f:

    for line in f:

        line = line.rstrip("\n")

        if not line:
            continue

        orthogroup, full_id = line.split("\t")

        if orthogroup in selected_groups:
            members.append(
                (orthogroup, full_id)
            )

print(
    "Final soluble members:",
    len(members)
)

# ------------------------------------------------------------
# Write member table
# ------------------------------------------------------------

with OUT_MEMBERS.open("w") as out:

    out.write(
        "orthogroup\tfull_id\n"
    )

    for orthogroup, full_id in sorted(members):

        out.write(
            f"{orthogroup}\t{full_id}\n"
        )

# ------------------------------------------------------------
# Extract protein sequences
# ------------------------------------------------------------

wanted = {
    full_id
    for _, full_id in members
}

found = set()

with PROTEIN_FASTA.open() as inp, \
     OUT_FASTA.open("w") as out:

    keep = False

    for line in inp:

        if line.startswith(">"):

            full_id = (
                line[1:]
                .strip()
                .split()[0]
            )

            keep = full_id in wanted

            if keep:
                found.add(full_id)
                out.write(line)

        elif keep:
            out.write(line)

print(
    "Protein sequences written:",
    len(found)
)

print(
    "Missing protein sequences:",
    len(wanted - found)
)

assert len(selected_groups) == 217
assert len(found) == len(wanted)

print("Members table:", OUT_MEMBERS)
print("FASTA:", OUT_FASTA)