from pathlib import Path
from collections import defaultdict
import csv
import statistics
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
    "work/codon/membrane_vs_soluble_clustered_summary.tsv"
)

N_BOOT = 100000
SEED = 20260830


def load_membrane(path, label):

    by_group = defaultdict(list)

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            species_total = int(row["species_total"])
            species_with = int(row["species_with_feature"])

            var_text = row[
                f"variance_distance_{label}"
            ]

            if (
                species_total < 8
                or species_with < 3
                or var_text == ""
            ):
                continue

            group = row["orthogroup"]

            by_group[group].append(
                float(var_text)
            )

    return {
        group: statistics.median(vals)
        for group, vals in by_group.items()
    }


def load_soluble(path, label):

    by_group = defaultdict(list)

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            species_total = int(row["species_total"])
            species_with = int(row["species_with_feature"])

            var_text = row[
                f"variance_distance_{label}"
            ]

            if (
                species_total < 8
                or species_with < 3
                or var_text == ""
            ):
                continue

            # Important:
            # cluster by the matched membrane orthogroup
            group = row["membrane_orthogroup"]

            by_group[group].append(
                float(var_text)
            )

    return {
        group: statistics.median(vals)
        for group, vals in by_group.items()
    }


def percentile(values, q):

    values = sorted(values)

    pos = (len(values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)

    frac = pos - lo

    return (
        values[lo] * (1 - frac)
        + values[hi] * frac
    )


def bootstrap_difference(mem_values, sol_values):

    rng = random.Random(SEED)

    observed = (
        statistics.median(sol_values)
        - statistics.median(mem_values)
    )

    boot = []

    for _ in range(N_BOOT):

        mem_sample = [
            rng.choice(mem_values)
            for _ in range(len(mem_values))
        ]

        sol_sample = [
            rng.choice(sol_values)
            for _ in range(len(sol_values))
        ]

        diff = (
            statistics.median(sol_sample)
            - statistics.median(mem_sample)
        )

        boot.append(diff)

    ci_low = percentile(boot, 0.025)
    ci_high = percentile(boot, 0.975)

    # Bootstrap support for soluble > membrane
    p_nonpositive = (
        sum(x <= 0 for x in boot) + 1
    ) / (
        N_BOOT + 1
    )

    return (
        observed,
        ci_low,
        ci_high,
        p_nonpositive,
    )


rows_out = []

for label, mem_path, sol_path in [
    ("start", MEM_START, SOL_START),
    ("end", MEM_END, SOL_END),
]:

    mem = load_membrane(
        mem_path,
        label
    )

    sol = load_soluble(
        sol_path,
        label
    )

    mem_values = list(mem.values())
    sol_values = list(sol.values())

    shared_groups = set(mem) & set(sol)

    (
        observed_diff,
        ci_low,
        ci_high,
        bootstrap_p,
    ) = bootstrap_difference(
        mem_values,
        sol_values,
    )

    print()
    print(
        f"=== CLUSTERED {label.upper()} ==="
    )

    print(
        "Membrane qualifying orthogroups:",
        len(mem_values)
    )

    print(
        "Soluble qualifying orthogroups:",
        len(sol_values)
    )

    print(
        "Shared orthogroups:",
        len(shared_groups)
    )

    print(
        "Membrane orthogroup-level median variance:",
        statistics.median(mem_values)
    )

    print(
        "Soluble orthogroup-level median variance:",
        statistics.median(sol_values)
    )

    print(
        "Observed median difference "
        "(soluble - membrane):",
        observed_diff
    )

    print(
        "95% bootstrap CI:",
        ci_low,
        ci_high
    )

    print(
        "Bootstrap probability difference <= 0:",
        bootstrap_p
    )

    rows_out.append(
        (
            label,
            len(mem_values),
            len(sol_values),
            len(shared_groups),
            statistics.median(mem_values),
            statistics.median(sol_values),
            observed_diff,
            ci_low,
            ci_high,
            bootstrap_p,
        )
    )


with OUT.open("w") as out:

    out.write(
        "anchor_type\t"
        "membrane_qualifying_orthogroups\t"
        "soluble_qualifying_orthogroups\t"
        "shared_orthogroups\t"
        "membrane_median_variance\t"
        "soluble_median_variance\t"
        "soluble_minus_membrane\t"
        "bootstrap_CI_low\t"
        "bootstrap_CI_high\t"
        "bootstrap_probability_difference_le_zero\n"
    )

    for row in rows_out:

        out.write(
            "\t".join(map(str, row))
            + "\n"
        )

print()
print("Output:", OUT)