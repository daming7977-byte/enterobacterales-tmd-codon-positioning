from pathlib import Path
import csv

BEST_MATCHES = Path(
    "work/codon/strict_soluble_best_matches.tsv"
)

ALL_MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

OUT_MEMBERS = Path(
    "work/codon/strict_soluble_matched_members.tsv"
)

OUT_FASTA = Path(
    "work/codon/strict_soluble_matched_all.faa"
)

PROTEIN_FASTA = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

# ------------------------------------------------------------
# Load selected soluble orthogroups
# ------------------------------------------------------------

selected_groups = set()

with BEST_MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        selected_groups.add(row["candidate_orthogroup"])

print("Selected soluble orthogroups:", len(selected_groups))

# ------------------------------------------------------------
# Extract group members
# eligible_one_to_one_strict.tsv has no header:
# orthogroup    full_id
# ------------------------------------------------------------

members = []

with ALL_MEMBERS.open() as f:
    for line in f:
        line = line.rstrip("\n")

        if not line:
            continue

        group, full_id = line.split("\t")

        if group in selected_groups:
            members.append((group, full_id))

print("Matched soluble members:", len(members))

# ------------------------------------------------------------
# Write member table
# ------------------------------------------------------------

with OUT_MEMBERS.open("w") as out:
    out.write("orthogroup\tfull_id\n")

    for group, full_id in sorted(members):
        out.write(f"{group}\t{full_id}\n")

# ------------------------------------------------------------
# Load requested protein IDs
# ------------------------------------------------------------

wanted = {full_id for _, full_id in members}

# ------------------------------------------------------------
# Extract protein FASTA
# ------------------------------------------------------------

found = set()

with PROTEIN_FASTA.open() as inp, OUT_FASTA.open("w") as out:

    keep = False

    for line in inp:

        if line.startswith(">"):
            full_id = line[1:].strip().split()[0]

            keep = full_id in wanted

            if keep:
                found.add(full_id)
                out.write(line)

        elif keep:
            out.write(line)

print("Protein sequences written:", len(found))
print("Missing protein sequences:", len(wanted - found))

if wanted - found:
    print("Missing IDs:")
    for x in sorted(wanted - found):
        print(x)

print("Members table:", OUT_MEMBERS)
print("FASTA:", OUT_FASTA)