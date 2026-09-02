from pathlib import Path
import csv
import matplotlib.pyplot as plt

null_file = Path(
    "work/final_results/fixed_presence_position_null_summary_1000.tsv"
)

out_pdf = Path(
    "manuscript/SUPPLEMENTARY_FIGURE_S2_fixed_presence_null.pdf"
)

out_png = Path(
    "manuscript/SUPPLEMENTARY_FIGURE_S2_fixed_presence_null.png"
)

observed = 31.875
null_median = 93.1958

null_values = []

with null_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        null_values.append(float(row["median_variance"]))

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

ax.set_ylabel(
    "Resample count"
)

ax.set_title(
    "Fixed-presence conditional positional null",
    fontsize=14
)

ax.legend(
    frameon=False,
    loc="upper right"
)

ax.text(
    0.97,
    0.56,
    "Usable units = 80\n"
    "Conditional resamples = 1,000\n"
    "Empirical one-sided P = 9.99 × 10⁻⁴",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=10
)

fig.tight_layout()

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")

print("null values:", len(null_values))
print("output PDF:", out_pdf)
print("output PNG:", out_png)