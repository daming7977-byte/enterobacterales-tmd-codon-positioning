from pathlib import Path
from collections import defaultdict

GROUPS = Path(
    "work/codon/exact_species_new_unscreened_groups.txt"
)

MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

PROTEINS = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

OUT_MEMBERS = Path(
    "work/codon/exact_species_new_candidate_members.tsv"
)

OUT_FASTA = Path(
    "work/codon/exact_species_new_candidates_all.faa"
)

# ------------------------------------------------------------
# Load target groups
# ------------------------------------------------------------

with GROUPS.open() as f:
    target_groups = {
        x.strip()
        for x in f
        if x.strip()
    }

print("Target orthogroups:", len(target_groups))

# ------------------------------------------------------------
# Extract target members
# ------------------------------------------------------------

members = []

with MEMBERS.open() as f:
    for line in f:

        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")

        if group in target_groups:
            members.append(
                (group, full_id)
            )

print("Target protein members:", len(members))

# ------------------------------------------------------------
# Write member table
# ------------------------------------------------------------

with OUT_MEMBERS.open("w") as out:

    out.write(
        "orthogroup\tfull_id\n"
    )

    for group, full_id in sorted(members):
        out.write(
            f"{group}\t{full_id}\n"
        )

# ------------------------------------------------------------
# Extract protein FASTA
# ------------------------------------------------------------

wanted = {
    full_id
    for _, full_id in members
}

found = set()

with PROTEINS.open() as inp, \
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

print("Sequences written:", len(found))
print("Missing:", len(wanted - found))

assert len(target_groups) == 571
assert len(found) == len(wanted)

print("Members:", OUT_MEMBERS)
print("FASTA:", OUT_FASTA)