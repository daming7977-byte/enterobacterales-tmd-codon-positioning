from pathlib import Path
import csv
import matplotlib.pyplot as plt

infile = Path(
    "results/leave_one_orthogroup_out_start.tsv"
)

out_pdf = Path(
    "figures/SUPPLEMENTARY_FIGURE_S3_leave_one_orthogroup_out.pdf"
)

out_png = Path(
    "figures/SUPPLEMENTARY_FIGURE_S3_leave_one_orthogroup_out.png"
)

full_median = 22.70

rows = []

with infile.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append({
            "orthogroup": row["removed_orthogroup"],
            "removed_units": int(row["removed_units"]),
            "retained_units": int(row["retained_units"]),
            "loo_median": float(row["loo_median_variance"]),
            "change": float(row["change_from_full_median"]),
        })

# Sort by leave-one-out median for visualization
rows.sort(key=lambda x: x["change"])

x = list(range(1, len(rows) + 1))
y = [r["change"] for r in rows]

fig, ax = plt.subplots(figsize=(6.4, 4.5))

ax.scatter(
    x,
    y,
    s=28
)

ax.axhline(
    0,
    linestyle="--",
    linewidth=1.5,
    label="Full-data reference"
)

ax.set_xlabel("Orthogroup removed")

ax.set_ylabel(
    "Change in median TMD-start-relative variance (aa²)"
)

ax.set_title(
    "Leave-one-orthogroup-out robustness",
    fontsize=14
)

ax.legend(
    frameon=False,
    loc="upper left"
)

ax.text(
    0.97,
    0.08,
    f"Orthogroups tested = {len(rows)}\n"
    f"LOO median range = {min(r['loo_median'] for r in rows):.2f}–"
    f"{max(r['loo_median'] for r in rows):.2f} aa²\n"
    f"Largest absolute change = {max(abs(r['change']) for r in rows):.2f} aa²",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=10
)

fig.tight_layout()

fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")

min_row = min(rows, key=lambda r: r["loo_median"])
max_row = max(rows, key=lambda r: r["loo_median"])
largest_change = max(rows, key=lambda r: abs(r["change"]))

print("orthogroups:", len(rows))
print(
    "minimum:",
    f"{min_row['loo_median']:.6f}",
    min_row["orthogroup"]
)
print(
    "maximum:",
    f"{max_row['loo_median']:.6f}",
    max_row["orthogroup"]
)
print(
    "largest absolute change:",
    f"{abs(largest_change['change']):.6f}",
    largest_change["orthogroup"]
)
print("output PDF:", out_pdf)
print("output PNG:", out_png)