from pathlib import Path
from collections import defaultdict, Counter
import random
import statistics

ALIGN_DIR = Path("work/topology/alignment/output")
TMD_TABLE = Path("work/codon/tmd_segment_integration.tsv")
SEGMENT_TABLE = Path("work/codon/low_adaptation_segments.tsv")
POOL_TABLE = Path("work/codon/non_tmd_candidate_pools.tsv")

OUT_SUMMARY = Path("work/codon/non_tmd_anchor_null_symmetric_1000perm.tsv")
OUT_UNITS = Path("work/codon/non_tmd_anchor_null_symmetric_units_1000perm.tsv")

N_PERM = 1000
MIN_POOL = 10
SEED = 20260826

random.seed(SEED)

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


def alignment_to_residue_map(aligned_seq):
    reverse = {}
    residue_pos = 0

    for aln_pos, aa in enumerate(aligned_seq, start=1):
        if aa != "-":
            residue_pos += 1
            reverse[aln_pos] = residue_pos

    return reverse


# ------------------------------------------------------------
# TMD architecture
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
            (genome, protein_id, tmd_index, start, end)
        )

group_counts = defaultdict(list)

for (group, full_id), rows in protein_tmd_rows.items():
    group_counts[group].append(len(rows))

modal_count = {}

for group, counts in group_counts.items():
    modal_count[group] = Counter(counts).most_common(1)[0][0]

eligible = {}

for (group, full_id), rows in protein_tmd_rows.items():

    if len(rows) == modal_count[group]:
        eligible[(group, full_id)] = sorted(
            rows,
            key=lambda x: x[2]
        )


# ------------------------------------------------------------
# Low-adaptation segments
# ------------------------------------------------------------

segments = defaultdict(list)

with SEGMENT_TABLE.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]

        start = int(fields[3])
        end = int(fields[4])

        center = (start + end) / 2

        segments[(group, genome, protein_id)].append(center)


# ------------------------------------------------------------
# Candidate pools
# ------------------------------------------------------------

candidate_pools = defaultdict(list)

with POOL_TABLE.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        aln_pos = int(fields[1])

        candidate_pools[group].append(aln_pos)


usable_groups = set()

for group, pool in candidate_pools.items():

    if len(pool) < MIN_POOL:
        continue

    if group not in modal_count:
        continue

    if len(pool) < modal_count[group]:
        continue

    usable_groups.add(group)


print("Usable orthogroups:", len(usable_groups))


# ------------------------------------------------------------
# Member alignment maps
# ------------------------------------------------------------

group_members = defaultdict(list)

for group in sorted(usable_groups):

    aln_path = ALIGN_DIR / f"{group}.faa"

    if not aln_path.exists():
        continue

    seqs = read_fasta(aln_path)

    for full_id, aligned_seq in seqs.items():

        key = (group, full_id)

        if key not in eligible:
            continue

        rows = eligible[key]

        genome = rows[0][0]
        protein_id = rows[0][1]

        reverse = alignment_to_residue_map(aligned_seq)

        group_members[group].append(
            (genome, protein_id, reverse)
        )


# ------------------------------------------------------------
# Observed TMD statistic restricted to same usable groups
# ------------------------------------------------------------

observed_table = Path(
    "work/codon/modal_tmd_cross_species_conservation.tsv"
)

observed_variances = []

with observed_table.open() as f:
    next(f)

    for line in f:

        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        species_total = int(fields[3])
        species_with_feature = int(fields[4])

        if group not in usable_groups:
            continue

        if species_total < 8:
            continue

        if species_with_feature < 3:
            continue

        if not fields[9]:
            continue

        observed_variances.append(float(fields[9]))


observed_median = statistics.median(observed_variances)

print("Observed qualifying units:", len(observed_variances))
print("Observed restricted median variance:", observed_median)


# ------------------------------------------------------------
# Symmetric pseudo-anchor null
# ------------------------------------------------------------

summary_rows = []
unit_rows = []

for perm in range(1, N_PERM + 1):

    variances_this_perm = []

    for group in sorted(usable_groups):

        n_tmd = modal_count[group]

        pseudo_columns = random.sample(
            candidate_pools[group],
            n_tmd
        )

        # stable ordinal order along alignment
        pseudo_columns.sort()

        members = group_members[group]

        # assignments[pseudo_index][genome] = nearest distance
        assignments = [
            {}
            for _ in range(n_tmd)
        ]

        species_total_per_anchor = [
            0
            for _ in range(n_tmd)
        ]

        for genome, protein_id, reverse in members:

            anchor_residues = []

            for pseudo_col in pseudo_columns:

                if pseudo_col in reverse:
                    anchor_residues.append(
                        reverse[pseudo_col]
                    )
                else:
                    anchor_residues.append(None)

            for i, anchor in enumerate(anchor_residues):
                if anchor is not None:
                    species_total_per_anchor[i] += 1

            feature_centers = segments.get(
                (group, genome, protein_id),
                []
            )

            # Each segment is assigned to exactly ONE nearest pseudo-anchor
            for feature_center in feature_centers:

                candidates = []

                for i, anchor in enumerate(anchor_residues):

                    if anchor is None:
                        continue

                    distance = feature_center - anchor

                    candidates.append(
                        (
                            abs(distance),
                            i,
                            distance
                        )
                    )

                if not candidates:
                    continue

                candidates.sort(
                    key=lambda x: (x[0], x[1])
                )

                _, best_i, best_distance = candidates[0]

                previous = assignments[best_i].get(genome)

                if (
                    previous is None
                    or abs(best_distance) < abs(previous)
                ):
                    assignments[best_i][genome] = best_distance


        for i in range(n_tmd):

            species_total = species_total_per_anchor[i]

            distances = list(
                assignments[i].values()
            )

            species_with_feature = len(distances)

            if (
                species_total >= 8
                and
                species_with_feature >= 3
            ):

                var = statistics.variance(distances)

                variances_this_perm.append(var)

                unit_rows.append(
                    (
                        perm,
                        group,
                        i + 1,
                        pseudo_columns[i],
                        species_total,
                        species_with_feature,
                        var,
                    )
                )


    median_var = (
        statistics.median(variances_this_perm)
        if variances_this_perm
        else float("nan")
    )

    summary_rows.append(
        (
            perm,
            len(variances_this_perm),
            median_var,
        )
    )

    if perm % 10 == 0:

        print(
            f"Permutation {perm}/{N_PERM}: "
            f"units={len(variances_this_perm)}, "
            f"median_variance={median_var}"
        )


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

with OUT_SUMMARY.open("w") as out:

    out.write(
        "permutation\t"
        "qualifying_units\t"
        "median_variance\n"
    )

    for row in summary_rows:
        out.write("\t".join(map(str, row)) + "\n")


with OUT_UNITS.open("w") as out:

    out.write(
        "permutation\t"
        "orthogroup\t"
        "pseudo_anchor_index\t"
        "pseudo_alignment_position\t"
        "species_total\t"
        "species_with_feature\t"
        "variance\n"
    )

    for row in unit_rows:
        out.write("\t".join(map(str, row)) + "\n")


null_vals = [
    row[2]
    for row in summary_rows
    if row[2] == row[2]
]

extreme = sum(
    x <= observed_median
    for x in null_vals
)

print()
print("=== SYMMETRIC NON-TMD ANCHOR NULL ===")
print("Permutations:", len(null_vals))
print("Observed qualifying units:", len(observed_variances))
print("Observed restricted median variance:", observed_median)
print("Null median:", statistics.median(null_vals))
print("Null minimum:", min(null_vals))
print("Null maximum:", max(null_vals))
print("Null <= observed:", extreme)
print(
    "Empirical one-sided p:",
    (extreme + 1) / (len(null_vals) + 1)
)
print("Summary:", OUT_SUMMARY)
print("Units:", OUT_UNITS)
