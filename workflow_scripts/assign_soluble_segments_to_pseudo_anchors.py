from pathlib import Path
from collections import defaultdict
import csv

SEGMENTS = Path(
    "work/codon/final_soluble_low_adaptation_segments.tsv"
)

ANCHOR_POSITIONS = Path(
    "work/codon/soluble_pseudo_anchor_residue_positions.tsv"
)

START_ANCHORS = Path(
    "work/codon/soluble_pseudo_anchor_start.tsv"
)

END_ANCHORS = Path(
    "work/codon/soluble_pseudo_anchor_end.tsv"
)

OUT_START = Path(
    "work/codon/soluble_segment_nearest_pseudo_start.tsv"
)

OUT_END = Path(
    "work/codon/soluble_segment_nearest_pseudo_end.tsv"
)

QC_THRESHOLD = 0.05


# ------------------------------------------------------------
# Load which pseudo-anchors pass QC
# ------------------------------------------------------------

def load_qc_anchor_keys(path):

    keep = set()

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            diff = float(
                row["absolute_relative_difference"]
            )

            if diff <= QC_THRESHOLD:

                keep.add(
                    (
                        row["membrane_orthogroup"],
                        row["soluble_orthogroup"],
                        int(row["tmd_index"]),
                    )
                )

    return keep


start_keep = load_qc_anchor_keys(
    START_ANCHORS
)

end_keep = load_qc_anchor_keys(
    END_ANCHORS
)

print(
    "START anchors passing QC:",
    len(start_keep)
)

print(
    "END anchors passing QC:",
    len(end_keep)
)


# ------------------------------------------------------------
# Load per-sequence pseudo-anchor residue positions
# ------------------------------------------------------------

start_positions = defaultdict(list)
end_positions = defaultdict(list)

with ANCHOR_POSITIONS.open() as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        anchor_type = row["anchor_type"]

        membrane_group = row[
            "membrane_orthogroup"
        ]

        soluble_group = row[
            "soluble_orthogroup"
        ]

        tmd_index = int(
            row["tmd_index"]
        )

        genome = row["genome"]

        full_id = row["full_id"]

        anchor_position = int(
            row["pseudo_anchor_residue_position"]
        )

        key = (
            membrane_group,
            soluble_group,
            tmd_index,
        )

        record = (
            tmd_index,
            anchor_position,
            membrane_group,
            soluble_group,
        )

        seq_key = (
            soluble_group,
            genome,
            full_id.split("|", 1)[1],
        )

        if (
            anchor_type == "start"
            and key in start_keep
        ):
            start_positions[
                seq_key
            ].append(record)

        elif (
            anchor_type == "end"
            and key in end_keep
        ):
            end_positions[
                seq_key
            ].append(record)


for d in (start_positions, end_positions):

    for key in d:
        d[key].sort(
            key=lambda x: x[1]
        )


# ------------------------------------------------------------
# Load soluble segments
# ------------------------------------------------------------

segments = []

with SEGMENTS.open() as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        group = row["orthogroup"]
        genome = row["genome"]
        protein_id = row["protein_id"]

        start = int(
            row["start_codon"]
        )

        end = int(
            row["end_codon"]
        )

        center = (
            start + end
        ) / 2

        length = int(
            row["length_codons"]
        )

        mean_weight = float(
            row["mean_relative_weight"]
        )

        segments.append(
            (
                group,
                genome,
                protein_id,
                start,
                end,
                center,
                length,
                mean_weight,
            )
        )

print(
    "Soluble segments loaded:",
    len(segments)
)


# ------------------------------------------------------------
# Assign each segment to exactly ONE nearest pseudo-anchor
# ------------------------------------------------------------

def assign_segments(
    anchor_dict,
    outfile,
    label,
):

    rows_out = []

    assigned = 0
    no_anchor = 0

    for (
        group,
        genome,
        protein_id,
        seg_start,
        seg_end,
        seg_center,
        seg_length,
        mean_weight,
    ) in segments:

        seq_key = (
            group,
            genome,
            protein_id,
        )

        anchors = anchor_dict.get(
            seq_key,
            []
        )

        if not anchors:

            no_anchor += 1
            continue

        best = min(
            anchors,
            key=lambda x: (
                abs(
                    seg_center
                    - x[1]
                ),
                x[0],
            )
        )

        (
            tmd_index,
            anchor_position,
            membrane_group,
            soluble_group,
        ) = best

        distance = (
            seg_center
            - anchor_position
        )

        rows_out.append(
            (
                membrane_group,
                soluble_group,
                genome,
                protein_id,
                tmd_index,
                anchor_position,
                seg_start,
                seg_end,
                seg_center,
                seg_length,
                mean_weight,
                distance,
            )
        )

        assigned += 1

    rows_out.sort(
        key=lambda x: (
            x[0],
            x[1],
            x[2],
            x[3],
            x[4],
            x[6],
        )
    )

    with outfile.open("w") as out:

        out.write(
            "membrane_orthogroup\t"
            "soluble_orthogroup\t"
            "genome\t"
            "protein_id\t"
            "pseudo_tmd_index\t"
            "pseudo_anchor_position\t"
            "segment_start\t"
            "segment_end\t"
            "segment_center\t"
            "segment_length\t"
            "segment_mean_weight\t"
            f"distance_to_pseudo_{label}\n"
        )

        for row in rows_out:

            out.write(
                f"{row[0]}\t"
                f"{row[1]}\t"
                f"{row[2]}\t"
                f"{row[3]}\t"
                f"{row[4]}\t"
                f"{row[5]}\t"
                f"{row[6]}\t"
                f"{row[7]}\t"
                f"{row[8]:.6f}\t"
                f"{row[9]}\t"
                f"{row[10]:.6f}\t"
                f"{row[11]:.6f}\n"
            )

    print()
    print(
        f"=== {label.upper()} ==="
    )

    print(
        "Segments assigned:",
        assigned
    )

    print(
        "Segments with no eligible pseudo-anchor:",
        no_anchor
    )

    print(
        "Output:",
        outfile
    )


assign_segments(
    start_positions,
    OUT_START,
    "start",
)

assign_segments(
    end_positions,
    OUT_END,
    "end",
)