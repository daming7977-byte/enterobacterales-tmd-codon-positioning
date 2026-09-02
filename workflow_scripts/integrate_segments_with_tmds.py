from pathlib import Path
from collections import defaultdict

members_file = Path("work/codon/topology_pass225_members.tsv")
segments_file = Path("work/codon/low_adaptation_segments.tsv")
tmd_file = Path("work/topology/deeptmhmm/anchor244/TMRs.gff3")
out_file = Path("work/codon/tmd_segment_integration.tsv")

protein_group = {}

with members_file.open() as f:
    for line in f:
        if not line.strip():
            continue
        group, full_id = line.rstrip("\n").split("\t")
        protein_group[full_id] = group

segments = defaultdict(list)

with segments_file.open() as f:
    next(f)
    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]

        start = int(fields[3])
        end = int(fields[4])
        length = int(fields[5])
        mean_weight = float(fields[6])

        full_id = f"{genome}|{protein_id}"
        center = (start + end) / 2

        segments[full_id].append({
            "start": start,
            "end": end,
            "center": center,
            "length": length,
            "mean_weight": mean_weight,
        })

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

        tmds[protein].append((start, end))

for protein in tmds:
    tmds[protein].sort()

rows = []

proteins_without_segment = 0
proteins_with_segment = 0

for full_id, group in sorted(protein_group.items()):

    if full_id not in tmds:
        continue

    genome, protein_id = full_id.split("|", 1)
    protein_segments = segments.get(full_id, [])

    if protein_segments:
        proteins_with_segment += 1
    else:
        proteins_without_segment += 1

    for tmd_index, (tmd_start, tmd_end) in enumerate(tmds[full_id], start=1):

        tmd_center = (tmd_start + tmd_end) / 2

        if not protein_segments:
            rows.append((
                group, genome, protein_id, full_id,
                tmd_index, tmd_start, tmd_end, tmd_center,
                "", "", "", "", "", "", "", "", ""
            ))
            continue

        nearest = min(
            protein_segments,
            key=lambda s: abs(s["center"] - tmd_center)
        )

        distance_to_start = nearest["center"] - tmd_start
        distance_to_end = nearest["center"] - tmd_end
        distance_to_center = nearest["center"] - tmd_center

        if nearest["end"] < tmd_start:
            relation = "upstream"
        elif nearest["start"] > tmd_end:
            relation = "downstream"
        else:
            relation = "overlap"

        rows.append((
            group,
            genome,
            protein_id,
            full_id,
            tmd_index,
            tmd_start,
            tmd_end,
            tmd_center,
            nearest["start"],
            nearest["end"],
            nearest["center"],
            nearest["length"],
            nearest["mean_weight"],
            distance_to_start,
            distance_to_end,
            distance_to_center,
            relation,
        ))

with out_file.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\tfull_id\t"
        "tmd_index\ttmd_start\ttmd_end\ttmd_center\t"
        "segment_start\tsegment_end\tsegment_center\t"
        "segment_length\tsegment_mean_weight\t"
        "distance_to_tmd_start\t"
        "distance_to_tmd_end\t"
        "distance_to_tmd_center\t"
        "relation\n"
    )

    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")

print("Analysis proteins:", len(protein_group))
print("Proteins with TMDs:", len(tmds))
print("Proteins with >=1 segment:", proteins_with_segment)
print("Proteins without segment:", proteins_without_segment)
print("Protein-TMD rows:", len(rows))
print("Output:", out_file)