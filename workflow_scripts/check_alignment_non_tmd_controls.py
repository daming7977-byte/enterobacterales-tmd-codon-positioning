from pathlib import Path
from collections import defaultdict, Counter

ALIGN_DIR = Path("work/topology/alignment/output")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")
CDS_QC = Path("work/codon/cds_protein_qc.tsv")

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
# Sequence coordinate -> alignment coordinate
# ------------------------------------------------------------

def residue_to_alignment_map(aligned_seq):

    mapping = {}
    residue_pos = 0

    for aln_pos, aa in enumerate(aligned_seq, start=1):

        if aa != "-":
            residue_pos += 1
            mapping[residue_pos] = aln_pos

    return mapping

# ------------------------------------------------------------
# Is residue position >= MIN_DISTANCE away from every TMD?
# ------------------------------------------------------------

def residue_is_valid_non_tmd(pos, intervals):

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
# Per orthogroup feasibility
# ------------------------------------------------------------

results = []

groups_seen = set()
groups_with_candidate = 0
total_tmd_units = 0
matchable_tmd_units = 0

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

        mapping = residue_to_alignment_map(aligned_seq)

        # Reverse map: alignment position -> residue position
        reverse = {
            aln_pos: residue_pos
            for residue_pos, aln_pos in mapping.items()
        }

        members.append(
            (
                full_id,
                aligned_seq,
                rows,
                intervals,
                mapping,
                reverse,
            )
        )

    if not members:
        continue

    groups_seen.add(group)

    modal_n_tmd = modal_count[group]

    aln_length = len(members[0][1])

    # --------------------------------------------------------
    # Candidate homologous non-TMD alignment columns
    # --------------------------------------------------------

    candidate_columns = []

    for aln_pos in range(1, aln_length + 1):

        valid_members = 0
        residue_members = 0

        for (
            full_id,
            aligned_seq,
            rows,
            intervals,
            mapping,
            reverse,
        ) in members:

            aa = aligned_seq[aln_pos - 1]

            if aa == "-":
                continue

            residue_members += 1

            residue_pos = reverse[aln_pos]

            if residue_is_valid_non_tmd(
                residue_pos,
                intervals
            ):
                valid_members += 1

        member_fraction = (
            valid_members / len(members)
            if members else 0
        )

        if member_fraction >= MIN_MEMBER_FRACTION:

            candidate_columns.append(
                (
                    aln_pos,
                    valid_members,
                    member_fraction,
                )
            )

    if candidate_columns:
        groups_with_candidate += 1

    # --------------------------------------------------------
    # For every homologous TMD index:
    # find nearest candidate alignment column
    # to median TMD-center alignment position
    # --------------------------------------------------------

    for tmd_index in range(1, modal_n_tmd + 1):

        tmd_center_aln_positions = []

        for (
            full_id,
            aligned_seq,
            rows,
            intervals,
            mapping,
            reverse,
        ) in members:

            row = rows[tmd_index - 1]

            start = row[3]
            end = row[4]

            center = int(round((start + end) / 2))

            if center in mapping:
                tmd_center_aln_positions.append(
                    mapping[center]
                )

        if not tmd_center_aln_positions:
            continue

        total_tmd_units += 1

        tmd_center_aln_positions.sort()

        n = len(tmd_center_aln_positions)

        if n % 2:
            median_center = tmd_center_aln_positions[n // 2]
        else:
            median_center = (
                tmd_center_aln_positions[n // 2 - 1]
                + tmd_center_aln_positions[n // 2]
            ) / 2

        if candidate_columns:

            best = min(
                candidate_columns,
                key=lambda x: (
                    abs(x[0] - median_center),
                    x[0],
                )
            )

            best_pos = best[0]
            best_valid_members = best[1]
            best_fraction = best[2]

            matchable = True
            matchable_tmd_units += 1

            alignment_distance = abs(
                best_pos - median_center
            )

        else:

            best_pos = ""
            best_valid_members = ""
            best_fraction = ""
            alignment_distance = ""
            matchable = False

        results.append(
            (
                group,
                tmd_index,
                len(members),
                modal_n_tmd,
                median_center,
                len(candidate_columns),
                matchable,
                best_pos,
                best_valid_members,
                best_fraction,
                alignment_distance,
            )
        )

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

out_file = Path(
    "work/codon/alignment_non_tmd_control_feasibility.tsv"
)

with out_file.open("w") as out:

    out.write(
        "orthogroup\t"
        "tmd_index\t"
        "modal_members\t"
        "modal_tmd_count\t"
        "median_tmd_center_alignment_position\t"
        "candidate_non_tmd_columns\t"
        "matchable\t"
        "best_control_alignment_position\t"
        "valid_members_at_control\t"
        "valid_member_fraction\t"
        "alignment_distance\n"
    )

    for row in results:

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

print("Orthogroups evaluated:", len(groups_seen))
print(
    "Orthogroups with >=1 valid non-TMD alignment column:",
    groups_with_candidate
)
print("Homologous TMD units evaluated:", total_tmd_units)
print("Matchable homologous TMD units:", matchable_tmd_units)

if total_tmd_units:
    print(
        "Matchable fraction:",
        matchable_tmd_units / total_tmd_units
    )

print("Output:", out_file)