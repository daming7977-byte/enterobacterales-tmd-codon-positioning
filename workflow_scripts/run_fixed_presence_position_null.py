from pathlib import Path
from collections import defaultdict, Counter
import random
import statistics
import re

# ============================================================
# Files and parameters
# ============================================================

members_file = Path("work/codon/topology_pass225_members.tsv")
weights_file = Path("work/codon/species_codon_weights.tsv")
threshold_file = Path("work/codon/species_low_adaptation_thresholds.tsv")
tmd_table = Path("work/codon/tmd_segment_integration.tsv")
observed_segments_file = Path("work/codon/segment_nearest_tmd.tsv")
observed_conservation_file = Path(
    "work/codon/modal_tmd_cross_species_conservation.tsv"
)

out_summary = Path(
    "work/codon/fixed_presence_position_null_summary.tsv"
)

N_SHUFFLES = 1000
N_RESAMPLES = 1000
RANDOM_SEED = 20260825

random.seed(RANDOM_SEED)

# ============================================================
# Genetic code
# ============================================================

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

# ============================================================
# Load final target proteins
# ============================================================

targets = []

with members_file.open() as f:
    for line in f:
        if not line.strip():
            continue

        group, full_id = line.rstrip("\n").split("\t")
        genome, protein_id = full_id.split("|", 1)

        targets.append(
            (group, genome, protein_id, full_id)
        )

# ============================================================
# Load codon weights
# ============================================================

weights = {}

with weights_file.open() as f:
    next(f)

    for line in f:
        genome, aa, codon, count, weight = (
            line.rstrip("\n").split("\t")
        )

        weights[(genome, codon)] = float(weight)

# ============================================================
# Load low-adaptation thresholds
# ============================================================

thresholds = {}

with threshold_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")
        thresholds[fields[0]] = float(fields[2])

# ============================================================
# Load CDS sequences
# ============================================================

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
            cds_sequences[(genome, current_id)] = (
                "".join(seq_parts).upper()
            )

    with path.open() as f:

        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):

                save_record()

                header = line[1:]
                seq_parts = []

                pm = re.search(
                    r"\[protein_id=([^\]]+)\]",
                    header
                )

                if pm:
                    current_id = pm.group(1)
                else:
                    current_id = header.split()[0]

            else:
                seq_parts.append(line.strip())

        save_record()

# ============================================================
# Recover modal-topology proteins and TMD coordinates
# ============================================================

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

group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():
    group_counts[group].append(len(rows))

modal_count = {}

for group, counts in group_counts.items():
    modal_count[group] = Counter(counts).most_common(1)[0][0]

eligible_proteins = set()
tmds = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():

    if len(rows) != modal_count[group]:
        continue

    eligible_proteins.add(full_id)

    for row in rows:
        (
            genome,
            protein_id,
            tmd_index,
            tmd_start,
            tmd_end,
            tmd_center,
        ) = row

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

# ============================================================
# Build CDS codon lists
# ============================================================

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

# ============================================================
# Observed presence pattern
#
# key:
# (group, tmd_index, genome)
#
# We retain only primary observed units:
# species_total >= 8
# species_with_feature >= 3
# ============================================================

qualifying_units = set()

with observed_conservation_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        tmd_index = int(fields[1])
        species_total = int(fields[3])
        species_with_feature = int(fields[4])

        if species_total >= 8 and species_with_feature >= 3:
            qualifying_units.add((group, tmd_index))

observed_presence = defaultdict(set)

with observed_segments_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        tmd_index = int(fields[9])

        key = (group, tmd_index)

        if key in qualifying_units:
            observed_presence[key].add(genome)

# ============================================================
# Observed positional variance
# ============================================================

observed_variances = []

with observed_conservation_file.open() as f:
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

observed_median_variance = statistics.median(
    observed_variances
)

print(
    "Observed qualifying units:",
    len(observed_variances)
)

print(
    "Observed median variance:",
    observed_median_variance
)

# ============================================================
# Synonymous shuffle
# ============================================================

def synonymous_shuffle(codons):

    shuffled = list(codons)

    aa_positions = defaultdict(list)

    # Exclude initiation codon.
    for i in range(1, len(codons)):

        codon = codons[i]
        aa = CODON_TABLE.get(codon)

        if aa is None or aa == "*":
            continue

        aa_positions[aa].append(i)

    for aa, positions in aa_positions.items():

        pool = [codons[i] for i in positions]

        random.shuffle(pool)

        for pos, codon in zip(positions, pool):
            shuffled[pos] = codon

    return shuffled

# ============================================================
# Call low-adaptation segments
# ============================================================

def call_segments(genome, codons):

    low = set()

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
            low.add(idx)

    segments = []
    current = []

    for pos in range(2, len(codons) + 1):

        if pos in low:

            if current and pos != current[-1] + 1:

                if len(current) >= 3:
                    segments.append(
                        (current[0], current[-1])
                    )

                current = []

            current.append(pos)

        else:

            if current:

                if len(current) >= 3:
                    segments.append(
                        (current[0], current[-1])
                    )

                current = []

    if current and len(current) >= 3:
        segments.append(
            (current[0], current[-1])
        )

    return segments

# ============================================================
# Conditional positional pools
#
# For each observed-positive:
#
#   orthogroup x TMD x species
#
# collect null distances ONLY when a shuffled feature
# is assigned to that same TMD.
#
# This conditions explicitly on feature presence.
# ============================================================

conditional_pools = defaultdict(list)

for shuffle_i in range(1, N_SHUFFLES + 1):

    for full_id, info in protein_codons.items():

        genome = info["genome"]

        shuffled = synonymous_shuffle(
            info["codons"]
        )

        segments = call_segments(
            genome,
            shuffled
        )

        if not segments:
            continue

        # nearest segment per TMD in this shuffled protein
        assigned = {}

        for seg_start, seg_end in segments:

            seg_center = (
                seg_start + seg_end
            ) / 2

            candidates = []

            for row in tmds[full_id]:

                (
                    group,
                    genome2,
                    protein_id,
                    tmd_index,
                    tmd_start,
                    tmd_end,
                    tmd_center,
                ) = row

                d = seg_center - tmd_center

                candidates.append(
                    (
                        abs(d),
                        group,
                        tmd_index,
                        d,
                    )
                )

            if not candidates:
                continue

            candidates.sort()

            _, group, tmd_index, d = (
                candidates[0]
            )

            target_key = (
                group,
                tmd_index,
                genome
            )

            unit_key = (
                group,
                tmd_index
            )

            # Only collect distributions for species
            # that truly show the feature in observed data.
            if (
                genome
                not in observed_presence.get(
                    unit_key,
                    set()
                )
            ):
                continue

            previous = assigned.get(
                target_key
            )

            if (
                previous is None
                or abs(d) < abs(previous)
            ):
                assigned[target_key] = d

        for key, d in assigned.items():
            conditional_pools[key].append(d)

    if shuffle_i % 100 == 0:

        print(
            "Shuffle",
            shuffle_i,
            "/",
            N_SHUFFLES
        )

# ============================================================
# Diagnose conditional-pool coverage
# ============================================================

target_keys = []

for unit_key, genomes in observed_presence.items():

    group, tmd_index = unit_key

    for genome in genomes:

        target_keys.append(
            (group, tmd_index, genome)
        )

pool_sizes = [
    len(conditional_pools[key])
    for key in target_keys
]

print()
print("=== CONDITIONAL POOL QC ===")

print(
    "Observed positive species-unit targets:",
    len(target_keys)
)

print(
    "Targets with >=1 conditional draw:",
    sum(x >= 1 for x in pool_sizes)
)

print(
    "Targets with >=10 conditional draws:",
    sum(x >= 10 for x in pool_sizes)
)

print(
    "Targets with >=20 conditional draws:",
    sum(x >= 20 for x in pool_sizes)
)

print(
    "Targets with >=50 conditional draws:",
    sum(x >= 50 for x in pool_sizes)
)

print(
    "Median conditional pool size:",
    statistics.median(pool_sizes)
)

# ============================================================
# Retain units for which every observed-positive species
# has enough conditional null samples.
#
# Minimum 20 conditional samples per species-unit.
# ============================================================

usable_units = {}

MIN_POOL = 20

for unit_key, genomes in observed_presence.items():

    keys = [
        (
            unit_key[0],
            unit_key[1],
            genome
        )
        for genome in genomes
    ]

    if all(
        len(conditional_pools[key]) >= MIN_POOL
        for key in keys
    ):
        usable_units[unit_key] = list(genomes)

print()
print(
    "Usable observed units:",
    len(usable_units)
)

# ============================================================
# Observed median restricted to exactly the usable units
# ============================================================

observed_unit_variance = {}

with observed_conservation_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        key = (
            fields[0],
            int(fields[1])
        )

        if key not in usable_units:
            continue

        if fields[9] == "":
            continue

        observed_unit_variance[key] = float(
            fields[9]
        )

if not observed_unit_variance:

    raise RuntimeError(
        "No usable units after conditional-pool QC."
    )

observed_restricted_median = statistics.median(
    observed_unit_variance.values()
)

# ============================================================
# Conditional resampling
#
# For each observed unit:
# - preserve the exact observed-positive species
# - draw one random null position from each species-specific
#   conditional distribution
# - calculate variance
#
# Then take median variance across all usable units.
# ============================================================

null_medians = []

for resample in range(
    1,
    N_RESAMPLES + 1
):

    unit_variances = []

    for unit_key, genomes in usable_units.items():

        values = []

        for genome in genomes:

            pool_key = (
                unit_key[0],
                unit_key[1],
                genome
            )

            values.append(
                random.choice(
                    conditional_pools[
                        pool_key
                    ]
                )
            )

        if len(values) >= 3:

            unit_variances.append(
                statistics.variance(
                    values
                )
            )

    null_medians.append(
        statistics.median(
            unit_variances
        )
    )

# ============================================================
# Empirical test
# ============================================================

extreme = sum(
    x <= observed_restricted_median
    for x in null_medians
)

empirical_p = (
    extreme + 1
) / (
    len(null_medians) + 1
)

# ============================================================
# Write summary
# ============================================================

with out_summary.open("w") as out:

    out.write(
        "resample\tmedian_variance\n"
    )

    for i, value in enumerate(
        null_medians,
        start=1
    ):
        out.write(
            f"{i}\t{value:.6f}\n"
        )

# ============================================================
# Final
# ============================================================

print()
print("=== FIXED-PRESENCE POSITION NULL ===")

print(
    "Usable observed units:",
    len(usable_units)
)

print(
    "Observed restricted median variance:",
    observed_restricted_median
)

print(
    "Conditional-null median:",
    statistics.median(
        null_medians
    )
)

print(
    "Conditional-null minimum:",
    min(null_medians)
)

print(
    "Conditional-null maximum:",
    max(null_medians)
)

print(
    "Null <= observed:",
    extreme
)

print(
    "Empirical one-sided p:",
    empirical_p
)

print(
    "Output:",
    out_summary
)