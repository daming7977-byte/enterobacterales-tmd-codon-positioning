from pathlib import Path
from collections import defaultdict, Counter
import math

ALIGN_DIR = Path("work/topology/alignment/output")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")

MIN_DISTANCE = 15
MIN_MEMBER_FRACTION = 0.80
N_BINS = 10

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
# Residue <-> alignment maps
# ------------------------------------------------------------

def residue_to_alignment_map(aligned_seq):

    mapping = {}
    residue_pos = 0

    for aln_pos, aa in enumerate(aligned_seq, start=1):

        if aa != "-":
            residue_pos += 1
            mapping[residue_pos] = aln_pos

    return mapping


def alignment_to_residue_map(aligned_seq):

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
# Relative-position bin
# ------------------------------------------------------------

def get_bin(position, alignment_length):

    rel = (position - 1) / alignment_length

    b = int(rel * N_BINS)

    if b >= N_BINS:
        b = N_BINS - 1

    return b

# ------------------------------------------------------------
# Main feasibility analysis
# ------------------------------------------------------------

results = []

total_tmd_units = 0
same_bin_matchable = 0

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

        forward = residue_to_alignment_map(aligned_seq)
        reverse = alignment_to_residue_map(aligned_seq)

        members.append(
            (
                full_id,
                aligned_seq,
                rows,
                intervals,
                forward,
                reverse,
            )
        )

    if not members:
        continue

    aln_length = len(members[0][1])
    modal_n_tmd = modal_count[group]

    # --------------------------------------------------------
    # Find valid homologous non-TMD alignment columns
    # --------------------------------------------------------

    candidate_columns = []

    for aln_pos in range(1, aln_length + 1):

        valid_members = 0

        for (
            full_id,
            aligned_seq,
            rows,
            intervals,
            forward,
            reverse,
        ) in members:

            if aln_pos not in reverse:
                continue

            residue_pos = reverse[aln_pos]

            if valid_non_tmd(
                residue_pos,
                intervals
            ):
                valid_members += 1

        fraction = valid_members / len(members)

        if fraction >= MIN_MEMBER_FRACTION:

            candidate_columns.append(
                (
                    aln_pos,
                    get_bin(aln_pos, aln_length),
                    valid_members,
                    fraction,
                )
            )

    # --------------------------------------------------------
    # TMD units
    # --------------------------------------------------------

    for tmd_index in range(1, modal_n_tmd + 1):

        centers = []

        for (
            full_id,
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

        total_tmd_units += 1

        centers.sort()

        n = len(centers)

        if n % 2:
            median_center = centers[n // 2]
        else:
            median_center = (
                centers[n // 2 - 1]
                + centers[n // 2]
            ) / 2

        tmd_bin = get_bin(
            median_center,
            aln_length
        )

        same_bin_candidates = [
            x for x in candidate_columns
            if x[1] == tmd_bin
        ]

        matchable = bool(same_bin_candidates)

        if matchable:
            same_bin_matchable += 1

        results.append(
            (
                group,
                tmd_index,
                len(members),
                aln_length,
                median_center,
                tmd_bin,
                len(candidate_columns),
                len(same_bin_candidates),
                matchable,
            )
        )

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

out_file = Path(
    "work/codon/non_tmd_decile_feasibility.tsv"
)

with out_file.open("w") as out:

    out.write(
        "orthogroup\t"
        "tmd_index\t"
        "modal_members\t"
        "alignment_length\t"
        "median_tmd_center_alignment_position\t"
        "tmd_decile\t"
        "all_valid_non_tmd_columns\t"
        "same_decile_candidate_columns\t"
        "matchable\n"
    )

    for row in results:
        out.write(
            "\t".join(map(str, row)) + "\n"
        )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Homologous TMD units evaluated:", total_tmd_units)
print("Same-decile matchable units:", same_bin_matchable)

if total_tmd_units:
    print(
        "Same-decile matchable fraction:",
        same_bin_matchable / total_tmd_units
    )

print("Output:", out_file)