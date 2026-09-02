from pathlib import Path
from collections import defaultdict, Counter

ALIGN_DIR = Path("work/topology/alignment/output")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")

OUT_FILE = Path("work/codon/non_tmd_candidate_pools.tsv")

MIN_DISTANCE = 15
MIN_MEMBER_FRACTION = 0.80

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
            (
                genome,
                protein_id,
                tmd_index,
                start,
                end,
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
# FASTA reader
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

# ------------------------------------------------------------
# Coordinate mapping
# ------------------------------------------------------------

def build_reverse_map(aligned_seq):

    reverse = {}
    residue_pos = 0

    for aln_pos, aa in enumerate(aligned_seq, start=1):

        if aa != "-":
            residue_pos += 1
            reverse[aln_pos] = residue_pos

    return reverse

# ------------------------------------------------------------
# Valid non-TMD residue?
# ------------------------------------------------------------

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
# Build homologous non-TMD candidate pools
# ------------------------------------------------------------

rows_out = []

groups_with_candidates = 0
total_candidate_columns = 0

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

        reverse = build_reverse_map(aligned_seq)

        members.append(
            (
                full_id,
                aligned_seq,
                intervals,
                reverse,
            )
        )

    if not members:
        continue

    aln_length = len(members[0][1])

    candidate_count_this_group = 0

    for aln_pos in range(1, aln_length + 1):

        residue_members = 0
        valid_members = 0

        for (
            full_id,
            aligned_seq,
            intervals,
            reverse,
        ) in members:

            if aln_pos not in reverse:
                continue

            residue_members += 1

            residue_pos = reverse[aln_pos]

            if valid_non_tmd(
                residue_pos,
                intervals
            ):
                valid_members += 1

        occupancy_fraction = (
            residue_members / len(members)
        )

        valid_fraction = (
            valid_members / len(members)
        )

        if (
            occupancy_fraction >= MIN_MEMBER_FRACTION
            and
            valid_fraction >= MIN_MEMBER_FRACTION
        ):

            candidate_count_this_group += 1
            total_candidate_columns += 1

            rows_out.append(
                (
                    group,
                    aln_pos,
                    aln_length,
                    len(members),
                    residue_members,
                    valid_members,
                    occupancy_fraction,
                    valid_fraction,
                )
            )

    if candidate_count_this_group > 0:
        groups_with_candidates += 1

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

with OUT_FILE.open("w") as out:

    out.write(
        "orthogroup\t"
        "alignment_position\t"
        "alignment_length\t"
        "modal_members\t"
        "residue_members\t"
        "valid_non_tmd_members\t"
        "occupancy_fraction\t"
        "valid_fraction\n"
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

candidate_counts = Counter()

for row in rows_out:
    candidate_counts[row[0]] += 1

print("Orthogroups with candidate pools:", groups_with_candidates)
print("Total candidate alignment columns:", total_candidate_columns)

if candidate_counts:

    vals = sorted(candidate_counts.values())

    print("Minimum candidates per matched orthogroup:", min(vals))
    print("Median candidates per matched orthogroup:", vals[len(vals)//2])
    print("Maximum candidates per matched orthogroup:", max(vals))

print("Output:", OUT_FILE)