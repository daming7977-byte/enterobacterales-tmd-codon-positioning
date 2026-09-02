from pathlib import Path
from collections import defaultdict, Counter

tmd_table = Path("work/codon/tmd_segment_integration.tsv")
out_file = Path("work/codon/non_tmd_control_space_qc.tsv")

MIN_DISTANCE = 30

# ------------------------------------------------------------
# Recover per-protein TMD architecture
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
        tmd_start = int(fields[5])
        tmd_end = int(fields[6])
        tmd_center = float(fields[7])

        protein_tmd_rows[(group, full_id)].append(
            (
                genome,
                protein_id,
                tmd_index,
                tmd_start,
                tmd_end,
                tmd_center,
            )
        )

# ------------------------------------------------------------
# Keep only modal-topology proteins
# ------------------------------------------------------------

group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():
    group_counts[group].append(len(rows))

modal_count = {}

for group, counts in group_counts.items():
    modal_count[group] = Counter(counts).most_common(1)[0][0]

eligible = {}

for (group, full_id), rows in protein_tmd_rows.items():

    if len(rows) != modal_count[group]:
        continue

    rows = sorted(rows, key=lambda x: x[2])

    eligible[(group, full_id)] = rows

# ------------------------------------------------------------
# Approximate protein length
#
# Use the maximum coordinate seen among TMDs as a lower bound.
# For actual control construction we will use true protein length.
# Here we instead recover it from CDS QC.
# ------------------------------------------------------------

protein_lengths = {}

qc_file = Path("work/codon/cds_protein_qc.tsv")

with qc_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        genome = fields[1]
        protein_id = fields[2]
        full_id = fields[3]
        protein_length = int(fields[4])

        protein_lengths[full_id] = protein_length

# ------------------------------------------------------------
# Valid control positions
#
# A position is valid if:
# - inside protein bounds
# - outside every TMD
# - >= MIN_DISTANCE aa from every TMD boundary
# ------------------------------------------------------------

def valid_positions(length, tmd_intervals):

    valid = []

    for pos in range(1, length + 1):

        ok = True

        for start, end in tmd_intervals:

            # Inside TMD
            if start <= pos <= end:
                ok = False
                break

            # Too close to either TMD boundary
            distance = min(
                abs(pos - start),
                abs(pos - end),
            )

            if distance < MIN_DISTANCE:
                ok = False
                break

        if ok:
            valid.append(pos)

    return valid

# ------------------------------------------------------------
# For each protein-TMD:
# ask whether at least one valid non-TMD anchor exists.
#
# We also record the valid position closest to the
# relative position of the real TMD center.
# ------------------------------------------------------------

rows_out = []

total_tmds = 0
matchable_tmds = 0
proteins_matchable_all = 0

for (group, full_id), rows in sorted(eligible.items()):

    length = protein_lengths.get(full_id)

    if length is None:
        continue

    genome = rows[0][0]
    protein_id = rows[0][1]

    intervals = [
        (r[3], r[4])
        for r in rows
    ]

    valid = valid_positions(
        length,
        intervals,
    )

    all_matchable = True

    for (
        genome,
        protein_id,
        tmd_index,
        tmd_start,
        tmd_end,
        tmd_center,
    ) in rows:

        total_tmds += 1

        if valid:

            real_relative = tmd_center / length

            best = min(
                valid,
                key=lambda p: abs(
                    (p / length) - real_relative
                ),
            )

            relative_difference = abs(
                (best / length) - real_relative
            )

            matchable = True
            matchable_tmds += 1

        else:

            best = ""
            relative_difference = ""
            matchable = False
            all_matchable = False

        rows_out.append(
            (
                group,
                genome,
                protein_id,
                full_id,
                tmd_index,
                length,
                tmd_start,
                tmd_end,
                tmd_center,
                len(valid),
                matchable,
                best,
                relative_difference,
            )
        )

    if all_matchable:
        proteins_matchable_all += 1

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

with out_file.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\tfull_id\t"
        "tmd_index\tprotein_length\t"
        "tmd_start\ttmd_end\ttmd_center\t"
        "valid_control_positions\t"
        "matchable\t"
        "best_control_position\t"
        "relative_position_difference\n"
    )

    for row in rows_out:

        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Eligible modal-topology proteins:", len(eligible))
print("Protein-TMD units:", total_tmds)
print("Matchable protein-TMD units:", matchable_tmds)
print(
    "Matchable fraction:",
    matchable_tmds / total_tmds
    if total_tmds else 0
)
print(
    "Proteins with all TMDs matchable:",
    proteins_matchable_all
)
print("Output:", out_file)