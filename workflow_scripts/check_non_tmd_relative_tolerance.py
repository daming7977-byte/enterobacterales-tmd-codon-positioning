from pathlib import Path
from collections import defaultdict, Counter

ALIGN_DIR = Path("work/topology/alignment/output")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")

MIN_DISTANCE = 15
MIN_MEMBER_FRACTION = 0.80

TOLERANCES = [0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25]

# ------------------------------------------------------------
# Load TMD architecture
# ------------------------------------------------------------

protein_tmd_rows = defaultdict(list)

with TMD_TABLE.open() as f:
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

        protein_tmd_rows[(group, full_id)].append(
            (genome, protein_id, tmd_index, start, end)
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
    if len(rows) == modal_count[group]:
        eligible[(group, full_id)] = sorted(
            rows,
            key=lambda x: x[2]
        )

# ------------------------------------------------------------
# FASTA
# ------------------------------------------------------------

def read_fasta(path):

    seqs = {}
    current = None
    chunks = []

    with path.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)

                current = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line)

        if current is not None:
            seqs[current] = "".join(chunks)

    return seqs


def coordinate_maps(aligned_seq):

    forward = {}
    reverse = {}

    residue_pos = 0

    for aln_pos, aa in enumerate(aligned_seq, start=1):
        if aa != "-":
            residue_pos += 1
            forward[residue_pos] = aln_pos
            reverse[aln_pos] = residue_pos

    return forward, reverse


def valid_non_tmd(pos, intervals):

    for start, end in intervals:

        if start <= pos <= end:
            return False

        distance = min(
            abs(pos - start),
            abs(pos - end)
        )

        if distance < MIN_DISTANCE:
            return False

    return True


# ------------------------------------------------------------
# Store each TMD unit and its valid candidate columns
# ------------------------------------------------------------

units = []

for group in sorted(modal_count):

    aln_path = ALIGN_DIR / f"{group}.faa"

    if not aln_path.exists():
        continue

    seqs = read_fasta(aln_path)

    members = []

    for full_id, aligned_seq in seqs.items():

        key = (group, full_id)

        if key not in eligible:
            continue

        rows = eligible[key]

        intervals = [
            (r[3], r[4])
            for r in rows
        ]

        forward, reverse = coordinate_maps(aligned_seq)

        members.append(
            (
                aligned_seq,
                rows,
                intervals,
                forward,
                reverse,
            )
        )

    if not members:
        continue

    aln_length = len(members[0][0])

    candidate_columns = []

    for aln_pos in range(1, aln_length + 1):

        valid_members = 0

        for (
            aligned_seq,
            rows,
            intervals,
            forward,
            reverse,
        ) in members:

            if aln_pos not in reverse:
                continue

            residue_pos = reverse[aln_pos]

            if valid_non_tmd(residue_pos, intervals):
                valid_members += 1

        fraction = valid_members / len(members)

        if fraction >= MIN_MEMBER_FRACTION:
            candidate_columns.append(aln_pos)

    modal_n_tmd = modal_count[group]

    for tmd_index in range(1, modal_n_tmd + 1):

        centers = []

        for (
            aligned_seq,
            rows,
            intervals,
            forward,
            reverse,
        ) in members:

            row = rows[tmd_index - 1]

            start = row[3]
            end = row[4]

            center = int(round((start + end) / 2))

            if center in forward:
                centers.append(forward[center])

        if not centers:
            continue

        centers.sort()

        n = len(centers)

        if n % 2:
            median_center = centers[n // 2]
        else:
            median_center = (
                centers[n // 2 - 1]
                + centers[n // 2]
            ) / 2

        tmd_relative = median_center / aln_length

        candidate_relative = [
            c / aln_length
            for c in candidate_columns
        ]

        units.append(
            (
                group,
                tmd_index,
                tmd_relative,
                candidate_relative,
            )
        )

# ------------------------------------------------------------
# Test tolerances
# ------------------------------------------------------------

print("Homologous TMD units:", len(units))
print()
print("tolerance\tmatchable\tfraction")

for tolerance in TOLERANCES:

    matchable = 0

    for (
        group,
        tmd_index,
        tmd_relative,
        candidate_relative,
    ) in units:

        ok = any(
            abs(c - tmd_relative) <= tolerance
            for c in candidate_relative
        )

        if ok:
            matchable += 1

    print(
        f"{tolerance:.3f}\t"
        f"{matchable}\t"
        f"{matchable / len(units):.4f}"
    )
    