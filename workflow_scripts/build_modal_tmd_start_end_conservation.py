from pathlib import Path
from collections import defaultdict, Counter
import statistics

SEGMENT_ASSIGNMENT = Path("work/codon/segment_nearest_tmd.tsv")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")

OUT_START = Path("work/codon/modal_tmd_start_cross_species_conservation.tsv")
OUT_END = Path("work/codon/modal_tmd_end_cross_species_conservation.tsv")

# ------------------------------------------------------------
# Recover modal TMD count per orthogroup
# ------------------------------------------------------------

protein_tmd_counts = defaultdict(list)

with TMD_TABLE.open() as f:
    next(f)

    seen = defaultdict(set)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        full_id = fields[3]
        tmd_index = int(fields[4])

        seen[(group, full_id)].add(tmd_index)

for (group, full_id), indices in seen.items():
    protein_tmd_counts[group].append(len(indices))

modal_tmd_count = {}

for group, counts in protein_tmd_counts.items():
    modal_tmd_count[group] = Counter(counts).most_common(1)[0][0]

# ------------------------------------------------------------
# Species totals for each modal orthogroup x TMD unit
# ------------------------------------------------------------

species_present = defaultdict(set)

with TMD_TABLE.open() as f:
    next(f)

    protein_rows = defaultdict(list)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        full_id = fields[3]
        tmd_index = int(fields[4])

        protein_rows[(group, full_id)].append(
            (genome, tmd_index)
        )
eligible_proteins = set()

for (group, full_id), rows in protein_rows.items():

    n_tmd = len(rows)

    if group not in modal_tmd_count:
        continue

    if n_tmd != modal_tmd_count[group]:
        continue

    eligible_proteins.add(full_id)

    for genome, tmd_index in rows:
        species_present[(group, tmd_index)].add(genome)

# ------------------------------------------------------------
# Load feature assignments
#
# Keep one nearest feature per species for each
# orthogroup x TMD unit.
# ------------------------------------------------------------

start_features = defaultdict(dict)
end_features = defaultdict(dict)

with SEGMENT_ASSIGNMENT.open() as f:
    header = next(f).rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[idx["orthogroup"]]
        genome = fields[idx["genome"]]
        full_id = fields[idx["full_id"]]

        if full_id not in eligible_proteins:
            continue

        tmd_index = int(
            fields[idx["nearest_tmd_index"]]
        )

        distance_start = float(
            fields[idx["distance_to_tmd_start"]]
        )

        distance_end = float(
            fields[idx["distance_to_tmd_end"]]
        )

        key = (group, tmd_index)

        previous = start_features[key].get(genome)

        if (
            previous is None
            or abs(distance_start) < abs(previous)
        ):
            start_features[key][genome] = distance_start

        previous = end_features[key].get(genome)

        if (
            previous is None
            or abs(distance_end) < abs(previous)
        ):
            end_features[key][genome] = distance_end

# ------------------------------------------------------------
# Build one conservation table
# ------------------------------------------------------------

def build_table(feature_dict, out_file, label):

    rows_out = []

    for key in sorted(species_present):

        group, tmd_index = key

        species_total = len(
            species_present[key]
        )

        values_by_species = feature_dict.get(
            key,
            {}
        )

        vals = list(values_by_species.values())

        species_with_feature = len(vals)

        presence_fraction = (
            species_with_feature / species_total
            if species_total else 0
        )

        if len(vals) >= 2:
            median_distance = statistics.median(vals)
            mean_distance = statistics.mean(vals)
            sd = statistics.stdev(vals)
            var = statistics.variance(vals)
        elif len(vals) == 1:
            median_distance = vals[0]
            mean_distance = vals[0]
            sd = ""
            var = ""
        else:
            median_distance = ""
            mean_distance = ""
            sd = ""
            var = ""

        rows_out.append(
            (
                group,
                tmd_index,
                modal_tmd_count[group],
                species_total,
                species_with_feature,
                presence_fraction,
                median_distance,
                mean_distance,
                sd,
                var,
            )
        )

    with out_file.open("w") as out:

        out.write(
            "orthogroup\t"
            "tmd_index\t"
            "modal_tmd_count\t"
            "species_total\t"
            "species_with_feature\t"
            "presence_fraction\t"
            f"median_distance_{label}\t"
            f"mean_distance_{label}\t"
            f"sd_distance_{label}\t"
            f"variance_distance_{label}\n"
        )

        for row in rows_out:

            formatted = []

            for x in row:
                if isinstance(x, float):
                    formatted.append(f"{x:.6f}")
                else:
                    formatted.append(str(x))

            out.write(
                "\t".join(formatted) + "\n"
            )

    qualifying = [
        r for r in rows_out
        if r[3] >= 8
        and r[4] >= 3
        and r[9] != ""
    ]

    variances = [
        float(r[9])
        for r in qualifying
    ]

    print()
    print(f"=== {label.upper()} ANCHOR ===")
    print("Orthogroup-TMD units:", len(rows_out))
    print("Qualifying units:", len(qualifying))

    if variances:
        print(
            "Median variance:",
            statistics.median(variances)
        )

    print("Output:", out_file)


# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

build_table(
    start_features,
    OUT_START,
    "start"
)

build_table(
    end_features,
    OUT_END,
    "end"
)