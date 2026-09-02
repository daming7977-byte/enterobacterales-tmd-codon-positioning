from pathlib import Path
from collections import defaultdict, Counter
import statistics

tmd_table = Path("work/codon/tmd_segment_integration.tsv")
segment_table = Path("work/codon/segment_nearest_tmd.tsv")
out_file = Path("work/codon/modal_tmd_cross_species_conservation.tsv")

# protein -> its TMD rows
protein_tmds = defaultdict(list)

with tmd_table.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        full_id = fields[3]
        tmd_index = int(fields[4])

        protein_tmds[(group, full_id)].append(
            (genome, protein_id, tmd_index)
        )

# modal TMD count per orthogroup
group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmds.items():
    group_counts[group].append(len(rows))

modal_count = {}

for group, counts in group_counts.items():
    c = Counter(counts)
    modal_count[group] = c.most_common(1)[0][0]

# denominator: only proteins matching modal architecture
units = defaultdict(dict)
eligible_proteins = set()

for (group, full_id), rows in protein_tmds.items():

    if len(rows) != modal_count[group]:
        continue

    eligible_proteins.add(full_id)

    for genome, protein_id, tmd_index in rows:
        units[(group, tmd_index)][genome] = protein_id

# feature: nearest assigned segment, but only modal-topology proteins
features = defaultdict(dict)

with segment_table.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        full_id = fields[3]

        if full_id not in eligible_proteins:
            continue

        tmd_index = int(fields[9])
        distance_center = float(fields[15])

        key = (group, tmd_index)

        previous = features[key].get(genome)

        if previous is None or abs(distance_center) < abs(previous):
            features[key][genome] = distance_center

rows_out = []

for (group, tmd_index), species in sorted(units.items()):

    total = len(species)
    vals = list(features.get((group, tmd_index), {}).values())

    present = len(vals)
    fraction = present / total if total else 0

    median = statistics.median(vals) if vals else ""
    mean = statistics.mean(vals) if vals else ""

    if len(vals) >= 2:
        sd = statistics.stdev(vals)
        var = statistics.variance(vals)
    else:
        sd = ""
        var = ""

    rows_out.append(
        (
            group,
            tmd_index,
            modal_count[group],
            total,
            present,
            fraction,
            median,
            mean,
            sd,
            var,
        )
    )

with out_file.open("w") as out:

    out.write(
        "orthogroup\ttmd_index\tmodal_tmd_count\t"
        "species_total\tspecies_with_feature\tpresence_fraction\t"
        "median_distance_center\tmean_distance_center\t"
        "sd_distance_center\tvariance_distance_center\n"
    )

    for row in rows_out:
        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")

print("Eligible modal-topology proteins:", len(eligible_proteins))
print("Orthogroup-TMD units:", len(rows_out))
print("Units species_total >=8:", sum(r[3] >= 8 for r in rows_out))
print("Units feature >=2 species:", sum(r[4] >= 2 for r in rows_out))
print("Units feature >=3 species:", sum(r[4] >= 3 for r in rows_out))
print("Output:", out_file)