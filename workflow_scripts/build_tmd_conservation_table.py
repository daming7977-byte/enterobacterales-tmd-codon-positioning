from pathlib import Path
from collections import defaultdict
import statistics

tmd_table = Path("work/codon/tmd_segment_integration.tsv")
segment_table = Path("work/codon/segment_nearest_tmd.tsv")

out_file = Path("work/codon/tmd_cross_species_conservation.tsv")


# ------------------------------------------------------------
# Denominator:
# every protein x homologous TMD in the final analysis set
# ------------------------------------------------------------

tmd_units = defaultdict(dict)

with tmd_table.open() as f:
    header = next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        tmd_index = int(fields[4])

        key = (group, tmd_index)

        tmd_units[key][genome] = protein_id


# ------------------------------------------------------------
# Segments assigned to nearest TMD
#
# For a given protein x TMD, there may be >1 low-adaptation segment.
# For positional conservation, retain the one closest to TMD center.
# Presence remains binary for that protein x TMD.
# ------------------------------------------------------------

features = defaultdict(dict)

with segment_table.open() as f:
    header = next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]

        segment_start = int(fields[4])
        segment_end = int(fields[5])
        segment_center = float(fields[6])

        tmd_index = int(fields[9])

        distance_to_start = float(fields[13])
        distance_to_end = float(fields[14])
        distance_to_center = float(fields[15])

        key = (group, tmd_index)

        candidate = {
            "protein_id": protein_id,
            "segment_start": segment_start,
            "segment_end": segment_end,
            "segment_center": segment_center,
            "distance_start": distance_to_start,
            "distance_end": distance_to_end,
            "distance_center": distance_to_center,
        }

        previous = features[key].get(genome)

        if (
            previous is None
            or abs(distance_to_center) < abs(previous["distance_center"])
        ):
            features[key][genome] = candidate


# ------------------------------------------------------------
# Summarize cross-species conservation
# ------------------------------------------------------------

rows = []

for key in sorted(tmd_units):

    group, tmd_index = key

    total_species = len(tmd_units[key])

    present = features.get(key, {})

    species_with_feature = len(present)

    presence_fraction = (
        species_with_feature / total_species
        if total_species
        else 0
    )

    distances_center = [
        x["distance_center"]
        for x in present.values()
    ]

    distances_start = [
        x["distance_start"]
        for x in present.values()
    ]

    distances_end = [
        x["distance_end"]
        for x in present.values()
    ]

    if distances_center:

        median_center = statistics.median(distances_center)
        mean_center = statistics.mean(distances_center)

        median_start = statistics.median(distances_start)
        median_end = statistics.median(distances_end)

    else:

        median_center = ""
        mean_center = ""
        median_start = ""
        median_end = ""

    if len(distances_center) >= 2:

        sd_center = statistics.stdev(distances_center)
        variance_center = statistics.variance(distances_center)

    else:

        sd_center = ""
        variance_center = ""

    rows.append(
        (
            group,
            tmd_index,
            total_species,
            species_with_feature,
            presence_fraction,
            median_center,
            mean_center,
            sd_center,
            variance_center,
            median_start,
            median_end,
        )
    )


# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

with out_file.open("w") as out:

    out.write(
        "orthogroup\t"
        "tmd_index\t"
        "species_total\t"
        "species_with_feature\t"
        "presence_fraction\t"
        "median_distance_center\t"
        "mean_distance_center\t"
        "sd_distance_center\t"
        "variance_distance_center\t"
        "median_distance_start\t"
        "median_distance_end\n"
    )

    for row in rows:

        formatted = []

        for x in row:

            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")


print("Orthogroup-TMD units:", len(rows))

print(
    "Units with >=1 feature:",
    sum(r[3] >= 1 for r in rows)
)

print(
    "Units with feature in >=2 species:",
    sum(r[3] >= 2 for r in rows)
)

print(
    "Units with feature in >=3 species:",
    sum(r[3] >= 3 for r in rows)
)

print(
    "Units with presence fraction >=0.5:",
    sum(r[4] >= 0.5 for r in rows)
)

print(
    "Units with presence fraction >=0.8:",
    sum(r[4] >= 0.8 for r in rows)
)

print("Output:", out_file)