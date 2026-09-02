from pathlib import Path
import csv
import statistics
import matplotlib.pyplot as plt

start_file = Path("data/processed/soluble_pseudo_anchor_start.tsv")
end_file = Path("data/processed/soluble_pseudo_anchor_end.tsv")

out_pdf = Path(
    "figures/SUPPLEMENTARY_FIGURE_S4_soluble_pseudo_anchor_QC.pdf"
)

out_png = Path(
    "figures/SUPPLEMENTARY_FIGURE_S4_soluble_pseudo_anchor_QC.png"
)

threshold = 0.05


def read_errors(path):
    values = []

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            values.append(
                float(row["absolute_relative_difference"])
            )

    return values


start = read_errors(start_file)
end = read_errors(end_file)

start_pass = sum(x <= threshold for x in start)
end_pass = sum(x <= threshold for x in end)

print("START anchors:", len(start))
print("START median error:", statistics.median(start))
print("START <= 0.05:", start_pass)

print("END anchors:", len(end))
print("END median error:", statistics.median(end))
print("END <= 0.05:", end_pass)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10.5, 4.3),
    sharey=True
)

# START
ax = axes[0]

ax.hist(
    start,
    bins=35
)

ax.axvline(
    threshold,
    linestyle="--",
    linewidth=1.5,
    label="QC threshold = 0.05"
)

ax.set_title("A  TMD-start pseudo-anchors")
ax.set_xlabel("Absolute relative-position difference")
ax.set_xlim(0, 0.06)
ax.set_ylabel("Pseudo-anchor count")

ax.text(
    0.97,
    0.76,
    f"Anchors = {len(start):,}\n"
    f"Median error = {statistics.median(start):.5f}\n"
    f"Passing QC = {start_pass:,}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=9
)
ax.set_xlim(0, 0.06)
ax.legend(
    frameon=False,
    loc="upper right"
)
ax.set_xlim(0, 0.06)
# END
ax = axes[1]

ax.hist(
    end,
    bins=35
)

ax.axvline(
    threshold,
    linestyle="--",
    linewidth=1.5,
    label="QC threshold = 0.05"
)

ax.set_title("B  TMD-end pseudo-anchors")
ax.set_xlabel("Absolute relative-position difference")
ax.set_xlim(0, 0.06)

ax.text(
    0.97,
    0.76,
    f"Anchors = {len(end):,}\n"
    f"Median error = {statistics.median(end):.5f}\n"
    f"Passing QC = {end_pass:,}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=9
)

ax.legend(
    frameon=False,
    loc="upper right"
)

fig.suptitle(
    "Matched soluble pseudo-anchor positional QC",
    fontsize=14
)

fig.tight_layout()

fig.savefig(
    out_pdf,
    bbox_inches="tight"
)

fig.savefig(
    out_png,
    dpi=300,
    bbox_inches="tight"
)

print("output PDF:", out_pdf)
print("output PNG:", out_png)