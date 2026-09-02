from pathlib import Path
from collections import defaultdict, Counter
import csv
import statistics

MATCHES = Path(
    "work/codon/exact_species_final_one_to_one_min_cost_matches.tsv"
)

TMD_TABLE = Path(
    "work/codon/tmd_segment_integration.tsv"
)

PROTEIN_FASTA = Path(
    "work/mmseqs_clean/all_proteins_chromosomal.faa"
)

ALIGN_DIR = Path(
    "work/codon/soluble_alignments/output"
)

OUT_START = Path(
    "work/codon/soluble_pseudo_anchor_start.tsv"
)

OUT_END = Path(
    "work/codon/soluble_pseudo_anchor_end.tsv"
)

OUT_POSITIONS = Path(
    "work/codon/soluble_pseudo_anchor_residue_positions.tsv"
)


# ------------------------------------------------------------
# FASTA reader
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
# Relative coordinate:
# first residue = 0
# last residue  = 1
# ------------------------------------------------------------

def relative_position(pos, length):

    if length <= 1:
        return 0.0

    return (pos - 1) / (length - 1)


# ------------------------------------------------------------
# Final exact-species membrane-soluble pairs
# ------------------------------------------------------------

pairs = {}

with MATCHES.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        pairs[
            row["membrane_orthogroup"]
        ] = row["candidate_orthogroup"]

print("Final matched pairs:", len(pairs))


# ------------------------------------------------------------
# Protein lengths
# ------------------------------------------------------------

all_proteins = read_fasta(PROTEIN_FASTA)

protein_lengths = {
    full_id: len(seq)
    for full_id, seq in all_proteins.items()
}


# ------------------------------------------------------------
# Read membrane TMD architecture
# ------------------------------------------------------------

protein_tmd_rows = defaultdict(list)

with TMD_TABLE.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        group = row["orthogroup"]

        if group not in pairs:
            continue

        full_id = row["full_id"]

        protein_tmd_rows[
            (group, full_id)
        ].append(
            (
                int(row["tmd_index"]),
                int(row["tmd_start"]),
                int(row["tmd_end"]),
            )
        )


# ------------------------------------------------------------
# Modal TMD count
# ------------------------------------------------------------

group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():

    indices = {
        x[0]
        for x in rows
    }

    group_counts[group].append(
        len(indices)
    )

modal_count = {}

for group, counts in group_counts.items():

    freq = Counter(counts)

    modal_count[group] = sorted(
        freq.items(),
        key=lambda x: (-x[1], x[0])
    )[0][0]

print(
    "Groups with modal TMD count:",
    len(modal_count)
)


# ------------------------------------------------------------
# Membrane relative-position targets
# only modal-topology proteins
# ------------------------------------------------------------

start_values = defaultdict(list)
end_values = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():

    if group not in modal_count:
        continue

    indices = {
        x[0]
        for x in rows
    }

    if len(indices) != modal_count[group]:
        continue

    length = protein_lengths[full_id]

    for tmd_index, start, end in rows:

        start_values[
            (group, tmd_index)
        ].append(
            relative_position(
                start,
                length
            )
        )

        end_values[
            (group, tmd_index)
        ].append(
            relative_position(
                end,
                length
            )
        )


start_targets = {}
end_targets = {}

for key, values in start_values.items():
    start_targets[key] = statistics.median(values)

for key, values in end_values.items():
    end_targets[key] = statistics.median(values)


# ------------------------------------------------------------
# Ordered minimum-cost assignment
#
# targets:
#   ordered relative TMD positions
#
# candidates:
#   ordered 100%-occupancy soluble columns
#
# Choose one unique candidate per target,
# preserving order, minimizing total abs difference.
# ------------------------------------------------------------

def ordered_match(targets, candidates):

    m = len(targets)
    n = len(candidates)

    if n < m:
        raise RuntimeError(
            f"Not enough candidate columns: "
            f"{n} for {m} anchors"
        )

    INF = float("inf")

    dp = [
        [INF] * (n + 1)
        for _ in range(m + 1)
    ]

    take = [
        [False] * (n + 1)
        for _ in range(m + 1)
    ]

    for j in range(n + 1):
        dp[0][j] = 0.0

    for i in range(1, m + 1):

        for j in range(1, n + 1):

            # Skip candidate j
            best = dp[i][j - 1]
            choose = False

            # Match target i to candidate j
            prev = dp[i - 1][j - 1]

            if prev != INF:

                cost = (
                    prev
                    + abs(
                        targets[i - 1][1]
                        - candidates[j - 1][1]
                    )
                )

                if cost < best:

                    best = cost
                    choose = True

            dp[i][j] = best
            take[i][j] = choose

    if dp[m][n] == INF:
        raise RuntimeError(
            "Ordered matching failed"
        )

    selected = []

    i = m
    j = n

    while i > 0:

        if j <= 0:
            raise RuntimeError(
                "Backtracking failed"
            )

        if take[i][j]:

            selected.append(
                (
                    targets[i - 1],
                    candidates[j - 1],
                )
            )

            i -= 1
            j -= 1

        else:
            j -= 1

    selected.reverse()

    return selected, dp[m][n]


# ------------------------------------------------------------
# Alignment column information
# ------------------------------------------------------------

def build_alignment_columns(path):

    seqs = read_fasta(path)

    lengths = {
        len(seq)
        for seq in seqs.values()
    }

    if len(lengths) != 1:
        raise RuntimeError(
            f"Alignment length mismatch: {path}"
        )

    aln_len = next(iter(lengths))

    ungapped_lengths = {
        full_id: len(seq.replace("-", ""))
        for full_id, seq in seqs.items()
    }

    residue_counter = {
        full_id: 0
        for full_id in seqs
    }

    columns = []

    per_sequence_position = {}

    for col in range(aln_len):

        residue_positions = {}

        fully_occupied = True

        for full_id, seq in seqs.items():

            char = seq[col]

            if char == "-":

                fully_occupied = False

            else:

                residue_counter[full_id] += 1

                residue_positions[
                    full_id
                ] = residue_counter[full_id]

        if not fully_occupied:
            continue

        relative_positions = []

        for full_id, residue_pos in residue_positions.items():

            relative_positions.append(
                relative_position(
                    residue_pos,
                    ungapped_lengths[full_id]
                )
            )

        median_relative = statistics.median(
            relative_positions
        )

        # alignment column is written 1-based
        alignment_column = col + 1

        columns.append(
            (
                alignment_column,
                median_relative,
            )
        )

        per_sequence_position[
            alignment_column
        ] = residue_positions

    return (
        seqs,
        columns,
        per_sequence_position
    )


# ------------------------------------------------------------
# Build start/end pseudo anchors
# ------------------------------------------------------------

start_rows = []
end_rows = []
position_rows = []

total_start_cost = 0.0
total_end_cost = 0.0

for membrane_group, soluble_group in sorted(
    pairs.items()
):

    path = (
        ALIGN_DIR
        / f"{soluble_group}.faa"
    )

    (
        soluble_seqs,
        candidate_columns,
        per_sequence_position,
    ) = build_alignment_columns(path)

    n_tmd = modal_count[membrane_group]

    start_target_list = [
        (
            i,
            start_targets[
                (membrane_group, i)
            ]
        )
        for i in range(1, n_tmd + 1)
    ]

    end_target_list = [
        (
            i,
            end_targets[
                (membrane_group, i)
            ]
        )
        for i in range(1, n_tmd + 1)
    ]

    start_match, start_cost = ordered_match(
        start_target_list,
        candidate_columns
    )

    end_match, end_cost = ordered_match(
        end_target_list,
        candidate_columns
    )

    total_start_cost += start_cost
    total_end_cost += end_cost

    for (
        (tmd_index, target_rel),
        (column, pseudo_rel),
    ) in start_match:

        start_rows.append(
            (
                membrane_group,
                soluble_group,
                tmd_index,
                n_tmd,
                target_rel,
                column,
                pseudo_rel,
                abs(target_rel - pseudo_rel),
                len(soluble_seqs),
            )
        )

        for full_id, residue_pos in (
            per_sequence_position[column].items()
        ):

            genome = full_id.split("|", 1)[0]

            position_rows.append(
                (
                    "start",
                    membrane_group,
                    soluble_group,
                    tmd_index,
                    genome,
                    full_id,
                    column,
                    residue_pos,
                )
            )

    for (
        (tmd_index, target_rel),
        (column, pseudo_rel),
    ) in end_match:

        end_rows.append(
            (
                membrane_group,
                soluble_group,
                tmd_index,
                n_tmd,
                target_rel,
                column,
                pseudo_rel,
                abs(target_rel - pseudo_rel),
                len(soluble_seqs),
            )
        )

        for full_id, residue_pos in (
            per_sequence_position[column].items()
        ):

            genome = full_id.split("|", 1)[0]

            position_rows.append(
                (
                    "end",
                    membrane_group,
                    soluble_group,
                    tmd_index,
                    genome,
                    full_id,
                    column,
                    residue_pos,
                )
            )


# ------------------------------------------------------------
# Write anchor tables
# ------------------------------------------------------------

def write_anchor_table(path, rows):

    with path.open("w") as out:

        out.write(
            "membrane_orthogroup\t"
            "soluble_orthogroup\t"
            "tmd_index\t"
            "modal_tmd_count\t"
            "membrane_target_relative_position\t"
            "soluble_alignment_column\t"
            "soluble_median_relative_position\t"
            "absolute_relative_difference\t"
            "species_total\n"
        )

        for row in rows:

            out.write(
                f"{row[0]}\t"
                f"{row[1]}\t"
                f"{row[2]}\t"
                f"{row[3]}\t"
                f"{row[4]:.8f}\t"
                f"{row[5]}\t"
                f"{row[6]:.8f}\t"
                f"{row[7]:.8f}\t"
                f"{row[8]}\n"
            )


write_anchor_table(
    OUT_START,
    start_rows
)

write_anchor_table(
    OUT_END,
    end_rows
)


# ------------------------------------------------------------
# Write per-sequence pseudo-anchor residue positions
# ------------------------------------------------------------

with OUT_POSITIONS.open("w") as out:

    out.write(
        "anchor_type\t"
        "membrane_orthogroup\t"
        "soluble_orthogroup\t"
        "tmd_index\t"
        "genome\t"
        "full_id\t"
        "alignment_column\t"
        "pseudo_anchor_residue_position\n"
    )

    for row in position_rows:

        out.write(
            "\t".join(
                map(str, row)
            )
            + "\n"
        )


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

print()
print("Start pseudo-anchors:", len(start_rows))
print("End pseudo-anchors:", len(end_rows))

print(
    "Groups represented in start:",
    len({x[0] for x in start_rows})
)

print(
    "Groups represented in end:",
    len({x[0] for x in end_rows})
)

start_diffs = [
    x[7]
    for x in start_rows
]

end_diffs = [
    x[7]
    for x in end_rows
]

print(
    "Median START relative-position difference:",
    statistics.median(start_diffs)
)

print(
    "Maximum START relative-position difference:",
    max(start_diffs)
)

print(
    "Median END relative-position difference:",
    statistics.median(end_diffs)
)

print(
    "Maximum END relative-position difference:",
    max(end_diffs)
)

print()
print("Start output:", OUT_START)
print("End output:", OUT_END)
print("Per-sequence positions:", OUT_POSITIONS)