from pathlib import Path
import csv
import statistics

import matplotlib.pyplot as plt

infile = Path("results/membrane_vs_soluble_family_paired.tsv")
out_pdf = Path("figures/SUPPLEMENTARY_FIGURE_S5_family_level_differences.pdf")
out_png = Path("figures/SUPPLEMENTARY_FIGURE_S5_family_level_differences.png")

data = {
    "start": [],
    "end": [],
}

with infile.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        anchor = row["anchor_type"].strip().lower()

        if anchor not in data:
            continue

        diff = float(row["soluble_minus_membrane"])
        mem = float(row["membrane_family_median_variance"])
        sol = float(row["soluble_family_median_variance"])
        group = row["membrane_orthogroup"]

        data[anchor].append(
            {
                "orthogroup": group,
                "difference": diff,
                "membrane": mem,
                "soluble": sol,
            }
        )


def summarize(rows):
    diffs = [r["difference"] for r in rows]
    n = len(diffs)
    positive = sum(x > 0 for x in diffs)
    negative = sum(x < 0 for x in diffs)
    zero = sum(x == 0 for x in diffs)
    median_diff = statistics.median(diffs)
    return n, positive, negative, zero, median_diff


fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), constrained_layout=True)

panel_info = [
    ("start", "A  TMD-start"),
    ("end", "B  TMD-end"),
]

for ax, (anchor, title) in zip(axes, panel_info):
    rows = sorted(data[anchor], key=lambda x: x["difference"])
    diffs = [r["difference"] for r in rows]
    x = list(range(1, len(diffs) + 1))

    ax.scatter(x, diffs, s=45)
    ax.axhline(0, linestyle="--", linewidth=1.8)

    # symlog handles both large positive values and small/negative values
    ax.set_yscale("symlog", linthresh=10)

    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel("Matched family pair", fontsize=13)
    ax.set_ylabel("Soluble − membrane median variance (aa²)", fontsize=13)

    n, positive, negative, zero, median_diff = summarize(rows)

    ax.text(
        0.95,
        0.10,
        (
            f"Pairs = {n}\n"
            f"Positive = {positive}\n"
            f"Negative = {negative}\n"
            f"Zero = {zero}\n"
            f"Median difference = {median_diff:.2f} aa²"
        ),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
    )

    ax.set_xlim(0, len(diffs) + 1)


plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, dpi=300, bbox_inches="tight")

print("START pairs:", len(data["start"]))
print("END pairs:", len(data["end"]))
print("output PDF:", out_pdf)
print("output PNG:", out_png)