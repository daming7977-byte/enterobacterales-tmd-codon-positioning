from pathlib import Path
from collections import defaultdict
import statistics

GROUPS = Path(
    "work/orthology/groups/eligible_one_to_one_strict.tsv"
)

FASTA = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

MEMBRANE_GROUPS = Path(
    "work/topology/topology_pass225_groups.txt"
)

OUT = Path(
    "work/codon/soluble_control_matching_space.tsv"
)

# ------------------------------------------------------------
# Load protein lengths
# ------------------------------------------------------------

lengths = {}

current = None
seq = []

with FASTA.open() as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if current is not None:
                lengths[current] = len("".join(seq))

            current = line[1:].split()[0]
            seq = []

        else:
            seq.append(line)

    if current is not None:
        lengths[current] = len("".join(seq))

print("Protein lengths loaded:", len(lengths))

# ------------------------------------------------------------
# Load strict orthogroups
# ------------------------------------------------------------

members = defaultdict(list)

with GROUPS.open() as f:

    for line in f:

        line = line.rstrip("\n")

        if not line:
            continue

        group, full_id = line.split("\t")

        members[group].append(full_id)

print("Strict orthogroups:", len(members))

# ------------------------------------------------------------
# Membrane groups
# ------------------------------------------------------------

membrane_groups = set()

with MEMBRANE_GROUPS.open() as f:

    for line in f:

        group = line.strip()

        if group:
            membrane_groups.add(group)

print("Membrane groups:", len(membrane_groups))

# ------------------------------------------------------------
# Summarize all groups
# ------------------------------------------------------------

rows = []

missing_lengths = 0

for group, ids in sorted(members.items()):

    group_lengths = []

    genomes = set()

    for full_id in ids:

        genome = full_id.split("|", 1)[0]
        genomes.add(genome)

        if full_id in lengths:
            group_lengths.append(lengths[full_id])
        else:
            missing_lengths += 1

    if not group_lengths:
        continue

    median_length = statistics.median(group_lengths)

    rows.append(
        (
            group,
            len(genomes),
            len(ids),
            median_length,
            group in membrane_groups,
        )
    )

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

with OUT.open("w") as out:

    out.write(
        "orthogroup\t"
        "species_coverage\t"
        "member_count\t"
        "median_protein_length\t"
        "membrane_group\n"
    )

    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

mem = [
    row for row in rows
    if row[4]
]

nonmem = [
    row for row in rows
    if not row[4]
]

print()
print("Membrane groups represented:", len(mem))
print("Non-membrane candidate groups:", len(nonmem))
print("Missing protein lengths:", missing_lengths)

if mem:
    print()
    print(
        "Membrane median group length:",
        statistics.median(row[3] for row in mem)
    )

    print(
        "Membrane length range:",
        min(row[3] for row in mem),
        "-",
        max(row[3] for row in mem)
    )

if nonmem:
    print()
    print(
        "Candidate median group length:",
        statistics.median(row[3] for row in nonmem)
    )

    print(
        "Candidate length range:",
        min(row[3] for row in nonmem),
        "-",
        max(row[3] for row in nonmem)
    )

print()
print("Species coverage counts:")

for cov in [8, 9, 10]:

    m = sum(
        row[1] == cov
        for row in mem
    )

    n = sum(
        row[1] == cov
        for row in nonmem
    )

    print(
        cov,
        "membrane=",
        m,
        "candidate=",
        n
    )

print()
print("Output:", OUT)