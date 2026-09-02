from pathlib import Path
from collections import defaultdict, Counter

tmd_table = Path("work/codon/tmd_segment_integration.tsv")
qc_file = Path("work/codon/cds_protein_qc.tsv")

out_file = Path("work/codon/non_tmd_control_anchors.tsv")

MIN_DISTANCE = 15

# ------------------------------------------------------------
# Load TMD architecture
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
# Modal topology
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
        protein_length = int(fields[4])

        protein_lengths[full_id] = protein_length

# ------------------------------------------------------------
# Valid non-TMD positions
# ------------------------------------------------------------

def get_valid_positions(length, intervals):

    valid = []

    for pos in range(1, length + 1):

        ok = True

        for start, end in intervals:

            if start <= pos <= end:
                ok = False
                break

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
# Match each real TMD to nearest valid relative-position anchor
#
# Codon features are NOT consulted here.
# ------------------------------------------------------------

rows_out = []

for (group, full_id), rows in sorted(eligible.items()):

    length = protein_lengths[full_id]

    intervals = [
        (r[3], r[4])
        for r in rows
    ]

    valid = get_valid_positions(
        length,
        intervals
    )

    for (
        genome,
        protein_id,
        tmd_index,
        tmd_start,
        tmd_end,
        tmd_center,
    ) in rows:

        if not valid:
            continue

        real_relative = tmd_center / length

        anchor = min(
            valid,
            key=lambda p: (
                abs((p / length) - real_relative),
                p,
            )
        )

        anchor_relative = anchor / length

        relative_difference = (
            anchor_relative - real_relative
        )

        nearest_tmd_boundary_distance = min(
            min(
                abs(anchor - start),
                abs(anchor - end)
            )
            for start, end in intervals
        )

        rows_out.append(
            (
                group,
                genome,
                protein_id,
                full_id,
                tmd_index,
                length,
                tmd_center,
                real_relative,
                anchor,
                anchor_relative,
                relative_difference,
                nearest_tmd_boundary_distance,
            )
        )

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

with out_file.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\tfull_id\t"
        "matched_tmd_index\tprotein_length\t"
        "real_tmd_center\treal_tmd_relative_position\t"
        "control_anchor\tcontrol_relative_position\t"
        "relative_position_difference\t"
        "nearest_tmd_boundary_distance\n"
    )

    for row in rows_out:

        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")

print("Eligible modal-topology proteins:", len(eligible))
print("Control anchors generated:", len(rows_out))
print("Output:", out_file)