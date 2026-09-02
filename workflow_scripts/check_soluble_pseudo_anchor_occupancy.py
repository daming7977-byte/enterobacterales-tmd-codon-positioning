from pathlib import Path
from collections import defaultdict, Counter
import csv
import math

MATCHES = Path(
    "work/codon/exact_species_final_one_to_one_min_cost_matches.tsv"
)

TMD_TABLE = Path(
    "work/codon/tmd_segment_integration.tsv"
)

ALIGN_DIR = Path(
    "work/codon/soluble_alignments/output"
)

# ------------------------------------------------------------
# Final membrane -> soluble pairs
# ------------------------------------------------------------

pairs = {}

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        pairs[row["membrane_orthogroup"]] = (
            row["candidate_orthogroup"]
        )

print("Final matched pairs:", len(pairs))

# ------------------------------------------------------------
# Determine modal TMD count for each membrane orthogroup
# ------------------------------------------------------------

protein_tmds = defaultdict(set)

with TMD_TABLE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        group = row["orthogroup"]

        if group not in pairs:
            continue

        full_id = row["full_id"]
        tmd_index = int(row["tmd_index"])

        protein_tmds[(group, full_id)].add(tmd_index)

group_counts = defaultdict(list)

for (group, full_id), indices in protein_tmds.items():
    group_counts[group].append(len(indices))

modal_tmd_count = {}

for group, counts in group_counts.items():

    freq = Counter(counts)

    modal = sorted(
        freq.items(),
        key=lambda x: (-x[1], x[0])
    )[0][0]

    modal_tmd_count[group] = modal

print(
    "Membrane groups with modal TMD count:",
    len(modal_tmd_count)
)

# ------------------------------------------------------------
# FASTA alignment reader
# ------------------------------------------------------------

def read_fasta(path):

    seqs = {}
    current = None
    parts = []

    with path.open() as f:

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if current is not None:
                    seqs[current] = "".join(parts)

                current = line[1:].split()[0]
                parts = []

            else:
                parts.append(line.strip())

        if current is not None:
            seqs[current] = "".join(parts)

    return seqs

# ------------------------------------------------------------
# Test occupancy thresholds
# ------------------------------------------------------------

thresholds = [1.00, 0.90, 0.80]

results = {
    t: {
        "feasible": 0,
        "not_feasible": 0,
        "candidate_counts": [],
    }
    for t in thresholds
}

problems = []

for membrane_group, soluble_group in sorted(pairs.items()):

    path = ALIGN_DIR / f"{soluble_group}.faa"

    seqs = read_fasta(path)

    nseq = len(seqs)

    lengths = {
        len(seq)
        for seq in seqs.values()
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"Alignment length mismatch: {soluble_group}"
        )

    aln_len = next(iter(lengths))

    required_anchors = modal_tmd_count[membrane_group]

    occupancy = []

    for col in range(aln_len):

        occupied = sum(
            seq[col] != "-"
            for seq in seqs.values()
        )

        occupancy.append(occupied)

    for threshold in thresholds:

        min_occ = math.ceil(
            threshold * nseq
        )

        candidate_columns = sum(
            x >= min_occ
            for x in occupancy
        )

        results[threshold][
            "candidate_counts"
        ].append(candidate_columns)

        if candidate_columns >= required_anchors:
            results[threshold]["feasible"] += 1

        else:
            results[threshold]["not_feasible"] += 1

            problems.append(
                (
                    threshold,
                    membrane_group,
                    soluble_group,
                    nseq,
                    required_anchors,
                    candidate_columns,
                )
            )

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

print()

for threshold in thresholds:

    r = results[threshold]

    vals = sorted(
        r["candidate_counts"]
    )

    n = len(vals)

    median = (
        vals[n // 2]
        if n % 2
        else (
            vals[n // 2 - 1]
            + vals[n // 2]
        ) / 2
    )

    print(
        f"Occupancy >= {threshold:.0%}"
    )

    print(
        "  feasible:",
        r["feasible"],
        "/",
        len(pairs)
    )

    print(
        "  not feasible:",
        r["not_feasible"]
    )

    print(
        "  median candidate columns:",
        median
    )

    print(
        "  minimum candidate columns:",
        min(vals)
    )

print()

if problems:

    print("First 20 failures:")

    for row in problems[:20]:
        print(*row)