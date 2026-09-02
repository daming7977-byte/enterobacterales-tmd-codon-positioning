from pathlib import Path
import statistics

observed_file = Path("work/codon/modal_tmd_cross_species_conservation.tsv")
null_file = Path("work/codon/synonymous_shuffle_null_units_1000perm.tsv")

# ------------------------------------------------------------
# Observed: species_total >=8, species_with_feature >=3, TMD index >1
# ------------------------------------------------------------

observed_variances = []

with observed_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        tmd_index = int(fields[1])
        species_total = int(fields[3])
        species_with_feature = int(fields[4])

        if tmd_index == 1:
            continue

        if species_total < 8:
            continue

        if species_with_feature < 3:
            continue

        if fields[9] == "":
            continue

        observed_variances.append(float(fields[9]))

observed_median = statistics.median(observed_variances)

# ------------------------------------------------------------
# Null: summarize each permutation after excluding TMD1
# ------------------------------------------------------------

perm_variances = {}

with null_file.open() as f:
    next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        perm = int(fields[0])
        tmd_index = int(fields[2])
        species_total = int(fields[3])
        species_with_feature = int(fields[4])
        variance = float(fields[5])

        if tmd_index == 1:
            continue

        if species_total < 8:
            continue

        if species_with_feature < 3:
            continue

        perm_variances.setdefault(perm, []).append(variance)

null_medians = []

for perm in sorted(perm_variances):
    values = perm_variances[perm]

    if values:
        null_medians.append(statistics.median(values))

extreme = sum(x <= observed_median for x in null_medians)

p = (extreme + 1) / (len(null_medians) + 1)

print("=== WITHOUT TMD1 ===")
print("Observed qualifying units:", len(observed_variances))
print("Observed median variance:", observed_median)
print("Null permutations:", len(null_medians))
print("Null median-of-medians:", statistics.median(null_medians))
print("Null minimum:", min(null_medians))
print("Null maximum:", max(null_medians))
print("Null <= observed:", extreme)
print("Empirical one-sided p:", p)