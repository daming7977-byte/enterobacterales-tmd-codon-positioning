from pathlib import Path
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
    "work/codon/soluble_pseudo_start_cross_species_conservation.tsv"
)

SOL_END = Path(
    "work/codon/soluble_pseudo_end_cross_species_conservation.tsv"
)

OUT = Path(
    "work/codon/membrane_vs_soluble_paired_comparison.tsv"
)

N_PERM = 100000
SEED = 20260830


def load_membrane(path, label):

    data = {}

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            species_total = int(
                row["species_total"]
            )

            species_with = int(
                row["species_with_feature"]
            )

            var_text = row[
                f"variance_distance_{label}"
            ]

            if (
                species_total < 8
                or species_with < 3
                or var_text == ""
            ):
                continue

            key = (
                row["orthogroup"],
                int(row["tmd_index"]),
            )

            data[key] = float(var_text)

    return data


def load_soluble(path, label):

    data = {}

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            species_total = int(
                row["species_total"]
            )

            species_with = int(
                row["species_with_feature"]
            )

            var_text = row[
                f"variance_distance_{label}"
            ]

            if (
                species_total < 8
                or species_with < 3
                or var_text == ""
            ):
                continue

            key = (
                row["membrane_orthogroup"],
                int(row["pseudo_tmd_index"]),
            )

            data[key] = float(var_text)

    return data


def exact_sign_test_greater(differences):

    positive = sum(x > 0 for x in differences)
    negative = sum(x < 0 for x in differences)

    n = positive + negative

    if n == 0:
        return positive, negative, 1.0

    # One-sided:
    # P(X >= positive), X ~ Binomial(n, 0.5)
    p = sum(
        math.comb(n, k)
        for k in range(positive, n + 1)
    ) / (2 ** n)

    return positive, negative, p


def median_signflip_test(differences):

    observed = statistics.median(differences)

    rng = random.Random(SEED)

    extreme = 0

    for _ in range(N_PERM):

        perm = [
            x if rng.random() < 0.5 else -x
            for x in differences
        ]

        stat = statistics.median(perm)

        if stat >= observed:
            extreme += 1

    p = (
        extreme + 1
    ) / (
        N_PERM + 1
    )

    return observed, p


all_rows = []

for label, mem_path, sol_path in [
    (
        "start",
        MEM_START,
        SOL_START,
    ),
    (
        "end",
        MEM_END,
        SOL_END,
    ),
]:

    mem = load_membrane(
        mem_path,
        label
    )

    sol = load_soluble(
        sol_path,
        label
    )

    shared = sorted(
        set(mem) & set(sol)
    )

    differences = [
        sol[key] - mem[key]
        for key in shared
    ]

    ratios = [
        sol[key] / mem[key]
        for key in shared
        if mem[key] > 0
    ]

    print()
    print(
        f"=== PAIRED {label.upper()} ==="
    )

    print(
        "Membrane qualifying units:",
        len(mem)
    )

    print(
        "Soluble qualifying units:",
        len(sol)
    )

    print(
        "Shared qualifying units:",
        len(shared)
    )

    if not shared:
        continue

    print(
        "Paired membrane median variance:",
        statistics.median(
            mem[key]
            for key in shared
        )
    )

    print(
        "Paired soluble median variance:",
        statistics.median(
            sol[key]
            for key in shared
        )
    )

    print(
        "Median paired difference "
        "(soluble - membrane):",
        statistics.median(
            differences
        )
    )

    if ratios:

        print(
            "Median soluble/membrane ratio:",
            statistics.median(ratios)
        )

    positive, negative, sign_p = \
        exact_sign_test_greater(
            differences
        )

    print(
        "Pairs soluble > membrane:",
        positive
    )

    print(
        "Pairs soluble < membrane:",
        negative
    )

    print(
        "One-sided exact sign-test p:",
        sign_p
    )

    median_diff, perm_p = \
        median_signflip_test(
            differences
        )

    print(
        "100000-permutation "
        "median-difference p:",
        perm_p
    )

    for key in shared:

        group, tmd_index = key

        all_rows.append(
            (
                label,
                group,
                tmd_index,
                mem[key],
                sol[key],
                sol[key] - mem[key],
            )
        )


with OUT.open("w") as out:

    out.write(
        "anchor_type\t"
        "membrane_orthogroup\t"
        "tmd_index\t"
        "membrane_variance\t"
        "soluble_variance\t"
        "soluble_minus_membrane\n"
    )

    for row in all_rows:

        out.write(
            f"{row[0]}\t"
            f"{row[1]}\t"
            f"{row[2]}\t"
            f"{row[3]:.6f}\t"
            f"{row[4]:.6f}\t"
            f"{row[5]:.6f}\n"
        )

print()
print("Output:", OUT)