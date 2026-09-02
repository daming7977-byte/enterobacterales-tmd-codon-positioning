from pathlib import Path
from collections import defaultdict, Counter

tmd_table = Path("work/codon/tmd_segment_integration.tsv")
qc_file = Path("work/codon/cds_protein_qc.tsv")

DISTANCES = [5, 10, 15, 20, 25, 30]

# ------------------------------------------------------------
# Load protein TMD architecture
# ------------------------------------------------------------

protein_tmd_rows = defaultdict(list)

with tmd_table.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        full_id = fields[3]

        tmd_index = int(fields[4])
        start = int(fields[5])
        end = int(fields[6])
        center = float(fields[7])

        protein_tmd_rows[(group, full_id)].append(
            (genome, protein_id, tmd_index, start, end, center)
        )

# ------------------------------------------------------------
# Keep modal-topology proteins
# ------------------------------------------------------------

group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():
    group_counts[group].append(len(rows))

modal_count = {}

for group, counts in group_counts.items():
    modal_count[group] = Counter(counts).most_common(1)[0][0]

eligible = {}

for (group, full_id), rows in protein_tmd_rows.items():

    if len(rows) == modal_count[group]:
        eligible[(group, full_id)] = sorted(
            rows,
            key=lambda x: x[2]
        )

# ------------------------------------------------------------
# Protein lengths
# ------------------------------------------------------------

protein_lengths = {}

with qc_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")
        full_id = fields[3]
        protein_lengths[full_id] = int(fields[4])

# ------------------------------------------------------------
# Test distance thresholds
# ------------------------------------------------------------

def has_valid_position(length, intervals, minimum_distance):

    for pos in range(1, length + 1):

        valid = True

        for start, end in intervals:

            if start <= pos <= end:
                valid = False
                break

            distance = min(
                abs(pos - start),
                abs(pos - end)
            )

            if distance < minimum_distance:
                valid = False
                break

        if valid:
            return True

    return False


total_units = sum(len(rows) for rows in eligible.values())

print("Eligible proteins:", len(eligible))
print("Total protein-TMD units:", total_units)
print()
print("distance\tmatchable\tfraction")

for distance in DISTANCES:

    matchable = 0

    for (group, full_id), rows in eligible.items():

        length = protein_lengths[full_id]

        intervals = [
            (row[3], row[4])
            for row in rows
        ]

        valid_exists = has_valid_position(
            length,
            intervals,
            distance
        )

        if valid_exists:
            matchable += len(rows)

    print(
        f"{distance}\t"
        f"{matchable}\t"
        f"{matchable / total_units:.4f}"
    )