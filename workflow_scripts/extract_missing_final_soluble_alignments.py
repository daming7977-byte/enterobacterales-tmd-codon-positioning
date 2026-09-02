from pathlib import Path
from collections import defaultdict
import csv

MATCHES = Path(
    "work/codon/exact_species_final_one_to_one_min_cost_matches.tsv"
)

MEMBERS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

PROTEINS = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

ALIGN_OUT = Path(
    "work/codon/soluble_alignments/output"
)

INPUT_OUT = Path(
    "work/codon/soluble_alignments/missing_input"
)

INPUT_OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Final soluble groups
# ------------------------------------------------------------

final_groups = set()

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        final_groups.add(
            row["candidate_orthogroup"]
        )

missing_groups = {
    g
    for g in final_groups
    if not (
        ALIGN_OUT / f"{g}.faa"
    ).exists()
}

print("Final soluble groups:", len(final_groups))
print("Missing alignment groups:", len(missing_groups))

# ------------------------------------------------------------
# Members of missing groups
# ------------------------------------------------------------

group_members = defaultdict(list)

with MEMBERS.open() as f:
    for line in f:

        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")

        if group in missing_groups:
            group_members[group].append(full_id)

wanted = {
    full_id
    for ids in group_members.values()
    for full_id in ids
}

# ------------------------------------------------------------
# Load protein sequences
# ------------------------------------------------------------

seqs = {}
current = None
parts = []

with PROTEINS.open() as f:

    for line in f:

        line = line.rstrip("\n")

        if line.startswith(">"):

            if current is not None and current in wanted:
                seqs[current] = "".join(parts)

            current = line[1:].split()[0]
            parts = []

        else:
            parts.append(line.strip())

    if current is not None and current in wanted:
        seqs[current] = "".join(parts)

# ------------------------------------------------------------
# Write one FASTA per orthogroup
# ------------------------------------------------------------

written = 0

for group in sorted(missing_groups):

    out = INPUT_OUT / f"{group}.faa"

    with out.open("w") as f:

        for full_id in sorted(group_members[group]):

            seq = seqs[full_id]

            f.write(f">{full_id}\n")

            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

            written += 1

print("FASTA files written:", len(missing_groups))
print("Sequences written:", written)
print("Output:", INPUT_OUT)