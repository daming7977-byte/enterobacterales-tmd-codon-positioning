from pathlib import Path
from collections import defaultdict
import csv
import statistics
import math
import random

MEM_START = Path(
    "work/codon/modal_tmd_start_cross_species_conservation.tsv"
)

MEM_END = Path(
    "work/codon/modal_tmd_end_cross_species_conservation.tsv"
)

SOL_START = Path(
    "work/codon/exact_species_final_soluble_pseudo_start_cross_species_conservation.tsv"
)

SOL_END = Path(
    "work/codon/exact_species_final_soluble_pseudo_end_cross_species_conservation.tsv"
)

OUT = Path(
    "work/codon/membrane_vs_soluble_family_paired.tsv"
)

N_PERM = 100000
SEED = 20260831


def load_membrane(path, label):
    values = defaultdict(list)

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if (
                int(row["species_total"]) < 8
                or int(row["species_with_feature"]) < 3
            ):
                continue

            x = row[f"variance_distance_{label}"]

            if x == "":
                continue

            values[row["orthogroup"]].append(float(x))

    return {
        g: statistics.median(v)
        for g, v in values.items()
    }


def load_soluble(path, label):
    values = defaultdict(list)

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if (
                int(row["species_total"]) < 8
                or int(row["species_with_feature"]) < 3
            ):
                continue

            x = row[f"variance_distance_{label}"]

            if x == "":
                continue

            # soluble control is indexed by its matched membrane family
            values[row["membrane_orthogroup"]].append(float(x))

    return {
        g: statistics.median(v)
        for g, v in values.items()
    }


def exact_sign_test_greater(diffs):
    nonzero = [d for d in diffs if d != 0]

    n = len(nonzero)
    k = sum(d > 0 for d in nonzero)

    if n == 0:
        return 1.0

    p = sum(
        math.comb(n, i)
        for i in range(k, n + 1)
    ) / (2 ** n)

    return p


def signflip_test(diffs):
    rng = random.Random(SEED)

    observed = statistics.median(diffs)

    count = 0

    for _ in range(N_PERM):
        permuted = [
            d if rng.random() < 0.5 else -d
            for d in diffs
        ]

        stat = statistics.median(permuted)

        if stat >= observed:
            count += 1

    return (
        count + 1
    ) / (
        N_PERM + 1
    )


rows_out = []

for label, mem_path, sol_path in [
    ("start", MEM_START, SOL_START),
    ("end", MEM_END, SOL_END),
]:

    mem = load_membrane(mem_path, label)
    sol = load_soluble(sol_path, label)

    shared = sorted(set(mem) & set(sol))

    pairs = [
        (
            g,
            mem[g],
            sol[g],
            sol[g] - mem[g],
        )
        for g in shared
    ]

    diffs = [x[3] for x in pairs]

    print()
    print(f"=== FAMILY-PAIRED {label.upper()} ===")
    print("Shared qualifying families:", len(pairs))

    if pairs:
        print(
            "Paired membrane median:",
            statistics.median(x[1] for x in pairs)
        )

        print(
            "Paired soluble median:",
            statistics.median(x[2] for x in pairs)
        )

        print(
            "Median paired difference (soluble - membrane):",
            statistics.median(diffs)
        )

        print(
            "Pairs soluble > membrane:",
            sum(d > 0 for d in diffs)
        )

        print(
            "Pairs soluble < membrane:",
            sum(d < 0 for d in diffs)
        )

        print(
            "Pairs tied:",
            sum(d == 0 for d in diffs)
        )

        print(
            "One-sided exact sign-test p:",
            exact_sign_test_greater(diffs)
        )

        print(
            "100000 family-level sign-flip p:",
            signflip_test(diffs)
        )

    for g, m, s, d in pairs:
        rows_out.append(
            (
                label,
                g,
                m,
                s,
                d,
            )
        )


with OUT.open("w") as out:
    out.write(
        "anchor_type\t"
        "membrane_orthogroup\t"
        "membrane_family_median_variance\t"
        "soluble_family_median_variance\t"
        "soluble_minus_membrane\n"
    )

    for row in rows_out:
        out.write(
            "\t".join(map(str, row))
            + "\n"
        )

print()
print("Output:", OUT)