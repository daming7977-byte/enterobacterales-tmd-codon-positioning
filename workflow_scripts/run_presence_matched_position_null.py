from pathlib import Path
from collections import defaultdict
import random
import statistics

observed_file = Path(
    "work/codon/modal_tmd_cross_species_conservation.tsv"
)

null_units_file = Path(
    "work/codon/synonymous_shuffle_null_units_1000perm.tsv"
)

out_file = Path(
    "work/codon/presence_matched_position_null.tsv"
)

RANDOM_SEED = 20260825
random.seed(RANDOM_SEED)

# ------------------------------------------------------------
# Observed qualifying units
# Keep:
# - species_total >= 8
# - species_with_feature >= 3
# ------------------------------------------------------------

observed_units = {}

with observed_file.open() as f:
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

        variance = float(fields[9])

        observed_units[(group, tmd_index)] = {
            "species_total": species_total,
            "presence_n": species_with_feature,
            "variance": variance,
        }

observed_variances = [
    x["variance"]
    for x in observed_units.values()
]

observed_median = statistics.median(observed_variances)

print("Observed qualifying units:", len(observed_units))
print("Observed median variance:", observed_median)

# ------------------------------------------------------------
# Load null unit-level data
#
# permutation -> (group, tmd_index) -> record
# ------------------------------------------------------------

null_data = defaultdict(dict)

with null_units_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        perm = int(fields[0])
        group = fields[1]
        tmd_index = int(fields[2])
        species_total = int(fields[3])
        species_with_feature = int(fields[4])
        variance = float(fields[5])

        null_data[perm][(group, tmd_index)] = {
            "species_total": species_total,
            "presence_n": species_with_feature,
            "variance": variance,
        }

# ------------------------------------------------------------
# Presence-matched comparison
#
# Important:
# The existing null-unit table only contains units with
# >=3 null species carrying a feature.
#
# Therefore exact matching is possible only when:
# null presence_n == observed presence_n.
#
# We do NOT substitute a different presence count.
# We do NOT relax the rule silently.
# ------------------------------------------------------------

summary_rows = []

for perm in sorted(null_data):

    matched_variances = []
    exact_matches = 0

    for key, obs in observed_units.items():

        null_record = null_data[perm].get(key)

        if null_record is None:
            continue

        if (
            null_record["species_total"]
            != obs["species_total"]
        ):
            continue

        if (
            null_record["presence_n"]
            != obs["presence_n"]
        ):
            continue

        matched_variances.append(
            null_record["variance"]
        )

        exact_matches += 1

    if matched_variances:
        median_variance = statistics.median(
            matched_variances
        )
    else:
        median_variance = ""

    summary_rows.append(
        (
            perm,
            exact_matches,
            median_variance,
        )
    )

# ------------------------------------------------------------
# Write raw presence-matched results
# ------------------------------------------------------------

with out_file.open("w") as out:

    out.write(
        "permutation\t"
        "exact_matched_units\t"
        "median_variance\n"
    )

    for perm, n, median_variance in summary_rows:

        out.write(
            f"{perm}\t{n}\t{median_variance}\n"
        )

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

match_counts = [
    row[1]
    for row in summary_rows
]

valid_null_medians = [
    row[2]
    for row in summary_rows
    if row[2] != ""
]

print()
print("=== PRESENCE-MATCH DIAGNOSTICS ===")

print(
    "Permutations:",
    len(summary_rows)
)

print(
    "Median exact matched units:",
    statistics.median(match_counts)
)

print(
    "Min exact matched units:",
    min(match_counts)
)

print(
    "Max exact matched units:",
    max(match_counts)
)

print(
    "Permutations with >=10 exact matches:",
    sum(n >= 10 for n in match_counts)
)

print(
    "Permutations with >=20 exact matches:",
    sum(n >= 20 for n in match_counts)
)

print(
    "Permutations with >=30 exact matches:",
    sum(n >= 30 for n in match_counts)
)

# ------------------------------------------------------------
# Conservative test:
#
# Only evaluate permutations with at least 20 exact
# presence-matched orthogroup-TMD units.
#
# This threshold is a QC rule for stability, not a
# biological proximity cutoff.
# ------------------------------------------------------------

stable_null = [
    row[2]
    for row in summary_rows
    if row[1] >= 20 and row[2] != ""
]

print()
print("=== STABLE MATCHED NULL ===")

print(
    "Stable permutations:",
    len(stable_null)
)

if stable_null:

    extreme = sum(
        x <= observed_median
        for x in stable_null
    )

    empirical_p = (
        (extreme + 1)
        / (len(stable_null) + 1)
    )

    print(
        "Observed median variance:",
        observed_median
    )

    print(
        "Matched-null median:",
        statistics.median(stable_null)
    )

    print(
        "Matched-null minimum:",
        min(stable_null)
    )

    print(
        "Matched-null maximum:",
        max(stable_null)
    )

    print(
        "Null <= observed:",
        extreme
    )

    print(
        "Empirical one-sided p:",
        empirical_p
    )

else:

    print(
        "Insufficient exact presence matches "
        "for a stable positional null."
    )

print()
print("Output:", out_file)