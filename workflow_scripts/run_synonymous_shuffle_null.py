from pathlib import Path
from collections import defaultdict, Counter
import random
import statistics
import re

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

members_file = Path("work/codon/topology_pass225_members.tsv")
weights_file = Path("work/codon/species_codon_weights.tsv")
tmd_table = Path("work/codon/tmd_segment_integration.tsv")
observed_table = Path("work/codon/modal_tmd_cross_species_conservation.tsv")

out_summary = Path("work/codon/synonymous_shuffle_null_summary.tsv")
out_unit = Path("work/codon/synonymous_shuffle_null_units.tsv")

N_PERMUTATIONS = 1000
RANDOM_SEED = 20260825

random.seed(RANDOM_SEED)

# ------------------------------------------------------------
# Genetic code
# ------------------------------------------------------------

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",

    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",

    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",

    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# ------------------------------------------------------------
# Load target proteins
# ------------------------------------------------------------

targets = []

with members_file.open() as f:
    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")
        genome, protein_id = full_id.split("|", 1)

        targets.append((group, genome, protein_id, full_id))

# ------------------------------------------------------------
# Load codon weights
# ------------------------------------------------------------

weights = {}

with weights_file.open() as f:
    next(f)

    for line in f:
        genome, aa, codon, count, weight = line.rstrip("\n").split("\t")
        weights[(genome, codon)] = float(weight)

# ------------------------------------------------------------
# Load species thresholds
# ------------------------------------------------------------

thresholds = {}

with Path("work/codon/species_low_adaptation_thresholds.tsv").open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")
        thresholds[fields[0]] = float(fields[2])

# ------------------------------------------------------------
# Load CDS
# ------------------------------------------------------------

cds_sequences = {}

for path in Path("data").rglob("cds_from_genomic.fna"):

    m = re.search(r"(GCF_\d+\.\d+)", str(path))

    if not m:
        continue

    genome = m.group(1)

    current_id = None
    seq_parts = []

    def save_record():
        if current_id is not None:
            cds_sequences[(genome, current_id)] = "".join(seq_parts).upper()

    with path.open() as f:

        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):

                save_record()

                header = line[1:]
                seq_parts = []

                pm = re.search(r"\[protein_id=([^\]]+)\]", header)

                if pm:
                    current_id = pm.group(1)
                else:
                    current_id = header.split()[0]

            else:
                seq_parts.append(line.strip())

        save_record()

# ------------------------------------------------------------
# Load modal-topology TMDs
# Determine modal number of TMDs per orthogroup
# ------------------------------------------------------------

protein_tmd_rows = defaultdict(list)

with tmd_table.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        full_id = fields[3]

        tmd_index = int(fields[4])
        tmd_start = int(fields[5])
        tmd_end = int(fields[6])
        tmd_center = float(fields[7])

        protein_tmd_rows[(group, full_id)].append(
            (
                genome,
                protein_id,
                tmd_index,
                tmd_start,
                tmd_end,
                tmd_center,
            )
        )

group_tmd_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():
    group_tmd_counts[group].append(len(rows))

modal_tmd_count = {}

for group, counts in group_tmd_counts.items():
    c = Counter(counts)
    modal_tmd_count[group] = c.most_common(1)[0][0]

eligible_proteins = set()
tmds = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():

    if len(rows) != modal_tmd_count[group]:
        continue

    eligible_proteins.add(full_id)

    for (
        genome,
        protein_id,
        tmd_index,
        tmd_start,
        tmd_end,
        tmd_center,
    ) in rows:

        tmds[full_id].append(
            (
                group,
                genome,
                protein_id,
                tmd_index,
                tmd_start,
                tmd_end,
                tmd_center,
            )
        )

# ------------------------------------------------------------
# Denominator for orthogroup x TMD
# ------------------------------------------------------------

unit_species = defaultdict(set)

for full_id in eligible_proteins:

    for row in tmds[full_id]:

        (
            group,
            genome,
            protein_id,
            tmd_index,
            tmd_start,
            tmd_end,
            tmd_center,
        ) = row

        unit_species[(group, tmd_index)].add(genome)

# Primary units require species_total >= 8
primary_units = {
    key
    for key, species in unit_species.items()
    if len(species) >= 8
}

print("Eligible modal-topology proteins:", len(eligible_proteins))
print("Primary orthogroup-TMD units:", len(primary_units))

# ------------------------------------------------------------
# Build codon lists for eligible proteins
# ------------------------------------------------------------

protein_codons = {}

for group, genome, protein_id, full_id in targets:

    if full_id not in eligible_proteins:
        continue

    cds = cds_sequences[(genome, protein_id)]

    codons = [
        cds[i:i+3]
        for i in range(0, len(cds) - 2, 3)
    ]

    protein_codons[full_id] = {
        "group": group,
        "genome": genome,
        "protein_id": protein_id,
        "codons": codons,
    }

# ------------------------------------------------------------
# Shuffle synonymous codons within each protein
# ------------------------------------------------------------

def synonymous_shuffle(codons):

    shuffled = list(codons)

    aa_positions = defaultdict(list)

    # Exclude initiation codon at position 0
    # Exclude stop codons
    for i in range(1, len(codons)):

        codon = codons[i]
        aa = CODON_TABLE.get(codon)

        if aa is None or aa == "*":
            continue

        aa_positions[aa].append(i)

    for aa, positions in aa_positions.items():

        codon_pool = [codons[i] for i in positions]

        random.shuffle(codon_pool)

        for pos, new_codon in zip(positions, codon_pool):
            shuffled[pos] = new_codon

    return shuffled

# ------------------------------------------------------------
# Call low-adaptation segments from shuffled codons
# ------------------------------------------------------------

def call_segments(genome, codons):

    low_positions = []

    # codon_position is 1-based
    for idx, codon in enumerate(codons, start=1):

        if idx == 1:
            continue

        aa = CODON_TABLE.get(codon)

        if aa is None or aa == "*":
            continue

        w = weights.get((genome, codon))

        if w is None:
            continue

        if w <= thresholds[genome]:
            low_positions.append(idx)

    low_set = set(low_positions)

    segments = []

    current = []

    max_pos = len(codons)

    for pos in range(2, max_pos + 1):

        if pos in low_set:

            if current and pos != current[-1] + 1:
                if len(current) >= 3:
                    segments.append((current[0], current[-1]))
                current = []

            current.append(pos)

        else:

            if current:
                if len(current) >= 3:
                    segments.append((current[0], current[-1]))
                current = []

    if current and len(current) >= 3:
        segments.append((current[0], current[-1]))

    return segments

# ------------------------------------------------------------
# For each shuffled segment:
# assign to exactly one nearest TMD
# ------------------------------------------------------------

def build_unit_features():

    features = defaultdict(dict)

    for full_id, info in protein_codons.items():

        genome = info["genome"]

        shuffled = synonymous_shuffle(info["codons"])

        segments = call_segments(genome, shuffled)

        if not segments:
            continue

        protein_tmds = tmds[full_id]

        for segment_start, segment_end in segments:

            segment_center = (segment_start + segment_end) / 2

            candidates = []

            for row in protein_tmds:

                (
                    group,
                    genome2,
                    protein_id,
                    tmd_index,
                    tmd_start,
                    tmd_end,
                    tmd_center,
                ) = row

                distance_center = segment_center - tmd_center

                candidates.append(
                    (
                        abs(distance_center),
                        group,
                        tmd_index,
                        distance_center,
                    )
                )

            candidates.sort()

            _, group, tmd_index, distance_center = candidates[0]

            key = (group, tmd_index)

            # A protein may have >1 segment assigned to same TMD.
            # Keep only the closest one for positional variance.
            previous = features[key].get(genome)

            if (
                previous is None
                or abs(distance_center) < abs(previous)
            ):
                features[key][genome] = distance_center

    return features

# ------------------------------------------------------------
# Observed statistic
# Only units:
# - species_total >= 8
# - species_with_feature >= 3
#
# Primary summary = median positional variance
# ------------------------------------------------------------

observed_variances = []

with observed_table.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        tmd_index = int(fields[1])
        species_total = int(fields[3])
        species_with_feature = int(fields[4])

        if species_total < 8:
            continue

        if species_with_feature < 3:
            continue

        if fields[9] == "":
            continue

        observed_variances.append(float(fields[9]))

observed_median_variance = statistics.median(observed_variances)

print("Observed qualifying units:", len(observed_variances))
print("Observed median variance:", observed_median_variance)

# ------------------------------------------------------------
# Permutations
# ------------------------------------------------------------

summary_rows = []
unit_rows = []

for perm in range(1, N_PERMUTATIONS + 1):

    features = build_unit_features()

    variances = []

    qualifying_units = 0

    for key in primary_units:

        values = list(features.get(key, {}).values())

        if len(values) < 3:
            continue

        variance = statistics.variance(values)

        variances.append(variance)
        qualifying_units += 1

        unit_rows.append(
            (
                perm,
                key[0],
                key[1],
                len(unit_species[key]),
                len(values),
                variance,
            )
        )

    if variances:
        median_variance = statistics.median(variances)
    else:
        median_variance = float("nan")

    summary_rows.append(
        (
            perm,
            qualifying_units,
            median_variance,
        )
    )

    if perm % 10 == 0:
        print(
            "Permutation",
            perm,
            "/",
            N_PERMUTATIONS,
            "qualifying_units=",
            qualifying_units,
            "median_variance=",
            median_variance,
        )

# ------------------------------------------------------------
# Empirical p-value
#
# Alternative hypothesis:
# observed variance is LOWER than null.
# ------------------------------------------------------------

null_medians = [
    x[2]
    for x in summary_rows
    if x[2] == x[2]
]

extreme = sum(
    x <= observed_median_variance
    for x in null_medians
)

empirical_p = (extreme + 1) / (len(null_medians) + 1)

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

with out_summary.open("w") as out:

    out.write(
        "permutation\tqualifying_units\tmedian_variance\n"
    )

    for row in summary_rows:
        out.write(
            f"{row[0]}\t{row[1]}\t{row[2]:.6f}\n"
        )

with out_unit.open("w") as out:

    out.write(
        "permutation\torthogroup\ttmd_index\t"
        "species_total\tspecies_with_feature\tvariance\n"
    )

    for row in unit_rows:
        out.write(
            f"{row[0]}\t{row[1]}\t{row[2]}\t"
            f"{row[3]}\t{row[4]}\t{row[5]:.6f}\n"
        )

print()
print("=== FINAL ===")
print("Permutations:", N_PERMUTATIONS)
print("Observed qualifying units:", len(observed_variances))
print("Observed median variance:", observed_median_variance)
print("Null median-of-medians:", statistics.median(null_medians))
print("Empirical one-sided p:", empirical_p)
print("Summary:", out_summary)
print("Unit table:", out_unit)