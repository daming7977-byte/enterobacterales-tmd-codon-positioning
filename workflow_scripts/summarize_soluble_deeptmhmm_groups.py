from pathlib import Path
from collections import defaultdict, Counter
import glob
import re

MEMBERS = Path(
    "work/codon/soluble_control_shortlist_members.tsv"
)

OUT_GROUPS = Path(
    "work/codon/soluble_control_group_topology_summary.tsv"
)

OUT_ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

# ------------------------------------------------------------
# full_id -> orthogroup
# ------------------------------------------------------------

id_to_group = {}

with MEMBERS.open() as f:
    next(f)

    for line in f:
        group, full_id = line.rstrip("\n").split("\t")
        id_to_group[full_id] = group

print("Members loaded:", len(id_to_group))

# ------------------------------------------------------------
# Read DeepTMHMM final classes
# ------------------------------------------------------------

predicted_class = {}

paths = sorted(
    glob.glob(
        "work/topology/deeptmhmm/soluble_batches/"
        "batch*/predicted_topologies.3line"
    )
)

for path in paths:

    with open(path) as f:

        for line in f:

            if not line.startswith(">"):
                continue

            text = line[1:].strip()

            if " | " not in text:
                continue

            full_id, cls = text.split(" | ", 1)

            predicted_class[full_id] = cls.strip()

print("Predictions loaded:", len(predicted_class))

# ------------------------------------------------------------
# Summarize by orthogroup
# ------------------------------------------------------------

group_classes = defaultdict(list)

for full_id, group in id_to_group.items():

    cls = predicted_class.get(full_id)

    if cls is None:
        group_classes[group].append("MISSING")
    else:
        group_classes[group].append(cls)

rows = []
all_glob_groups = []

for group in sorted(group_classes):

    classes = group_classes[group]

    counts = Counter(classes)

    n = len(classes)
    n_glob = counts["GLOB"]

    all_glob = (
        n > 0
        and n_glob == n
    )

    if all_glob:
        all_glob_groups.append(group)

    rows.append(
        (
            group,
            n,
            n_glob,
            counts["TM"],
            counts["BETA"],
            counts["SP"],
            counts["SP+TM"],
            counts["MISSING"],
            n_glob / n if n else 0,
            all_glob,
        )
    )

# ------------------------------------------------------------
# Write summary
# ------------------------------------------------------------

with OUT_GROUPS.open("w") as out:

    out.write(
        "orthogroup\t"
        "members\t"
        "GLOB\t"
        "TM\t"
        "BETA\t"
        "SP\t"
        "SP_plus_TM\t"
        "MISSING\t"
        "glob_fraction\t"
        "all_GLOB\n"
    )

    for row in rows:

        formatted = []

        for x in row:
            if isinstance(x, float):
                formatted.append(f"{x:.6f}")
            else:
                formatted.append(str(x))

        out.write("\t".join(formatted) + "\n")

with OUT_ALL_GLOB.open("w") as out:

    for group in all_glob_groups:
        out.write(group + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Orthogroups summarized:", len(rows))
print("All-GLOB orthogroups:", len(all_glob_groups))

print(
    "All-GLOB fraction:",
    len(all_glob_groups) / len(rows)
    if rows else 0
)

print("Output:", OUT_GROUPS)
print("All-GLOB list:", OUT_ALL_GLOB)