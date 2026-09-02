from pathlib import Path
import csv
import matplotlib.pyplot as plt

null_file = Path(
    "results/without_tmd1_null_medians_1000.tsv"
)

out_pdf = Path(
    "figures/SUPPLEMENTARY_FIGURE_S1_without_TMD1.pdf"
)

out_png = Path(
    "figures/SUPPLEMENTARY_FIGURE_S1_without_TMD1.png"
)

observed = 17.208333
null_median = 116.81249975
p_value = 0.000999000999000999

null_values = []

with null_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        null_values.append(
            float(row["median_variance_without_tmd1"])
        )

fig, ax = plt.subplots(figsize=(6.2, 4.5))

ax.hist(
    null_values,
    bins=35
)

ax.axvline(
    observed,
    linewidth=2,
    label=f"Observed = {observed:.2f} aa²"
)

ax.axvline(
    null_median,
    linestyle="--",
    linewidth=1.5,
    label=f"Null median = {null_median:.2f} aa²"
)

ax.set_xlabel(
    "Median cross-species positional variance (aa²)"
)

ax.set_ylabel("Permutation count")

ax.set_title(
    "Sensitivity analysis excluding first TMDs",
    fontsize=14
)

ax.legend(
    frameon=False,
    loc="upper right"
)

ax.text(
    0.97,
    0.60,
    "Qualifying units = 54\n"
    "Permutations = 1,000\n"
    "Empirical one-sided P = 9.99 × 10⁻⁴",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=10
)

ax.legend(frameon=False)

fig.tight_layout()

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")

print("null values:", len(null_values))
print("output PDF:", out_pdf)
print("output PNG:", out_png)