from pathlib import Path
from collections import defaultdict

members_file = Path("work/codon/topology_pass225_members.tsv")
segments_file = Path("work/codon/low_adaptation_segments.tsv")
tmd_file = Path("work/topology/deeptmhmm/anchor244/TMRs.gff3")

out_file = Path("work/codon/segment_nearest_tmd.tsv")

# ------------------------------------------------------------
# Final protein set
# ------------------------------------------------------------

protein_group = {}

with members_file.open() as f:
    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")
        protein_group[full_id] = group


# ------------------------------------------------------------
# TMDs
# ------------------------------------------------------------

tmds = defaultdict(list)

with tmd_file.open() as f:
    for line in f:

        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 4:
            fields = line.rstrip("\n").split()

        if len(fields) < 4:
            continue

        protein = fields[0]
        feature = fields[1]

        if feature != "TMhelix":
            continue

        if protein not in protein_group:
            continue

        start = int(fields[2])
        end = int(fields[3])
        center = (start + end) / 2

        tmds[protein].append((start, end, center))

for protein in tmds:
    tmds[protein].sort()


# ------------------------------------------------------------
# Assign every segment to exactly one nearest TMD
# ------------------------------------------------------------

rows = []

with segments_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]

        segment_start = int(fields[3])
        segment_end = int(fields[4])
        segment_length = int(fields[5])
        mean_weight = float(fields[6])

        segment_center = (segment_start + segment_end) / 2

        full_id = f"{genome}|{protein_id}"

        if full_id not in tmds:
            continue

        candidates = []

        for i, (start, end, center) in enumerate(tmds[full_id], start=1):

            distance_center = segment_center - center

            candidates.append(
                (
                    abs(distance_center),
                    i,
                    start,
                    end,
                    center,
                    distance_center,
                )
            )

        candidates.sort()

        _, tmd_index, tmd_start, tmd_end, tmd_center, distance_center = candidates[0]

        distance_start = segment_center - tmd_start
        distance_end = segment_center - tmd_end

        if segment_end < tmd_start:
            relation = "upstream"
        elif segment_start > tmd_end:
            relation = "downstream"
        else:
            relation = "overlap"

        rows.append(
            (
                group,
                genome,
                protein_id,
                full_id,
                segment_start,
                segment_end,
                segment_center,
                segment_length,
                mean_weight,
                tmd_index,
                tmd_start,
                tmd_end,
                tmd_center,
                distance_start,
                distance_end,
                distance_center,
                relation,
            )
        )


with out_file.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\tfull_id\t"
        "segment_start\tsegment_end\tsegment_center\t"
        "segment_length\tsegment_mean_weight\t"
        "nearest_tmd_index\tnearest_tmd_start\tnearest_tmd_end\t"
        "nearest_tmd_center\tdistance_to_tmd_start\t"
        "distance_to_tmd_end\tdistance_to_tmd_center\trelation\n"
    )

    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")


print("Segments input:", sum(1 for _ in open(segments_file)) - 1)
print("Segments assigned:", len(rows))
print("Output:", out_file)