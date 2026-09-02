from pathlib import Path
from collections import defaultdict
from collections import deque
import csv
import glob

CANDIDATES = Path(
    "work/codon/exact_species_full_soluble_candidates.tsv"
)

NEW_MEMBERS = Path(
    "work/codon/exact_species_new_candidate_members.tsv"
)

OLD_ALL_GLOB = Path(
    "work/codon/soluble_control_all_glob_groups.txt"
)

OUT_NEW_SUMMARY = Path(
    "work/codon/exact_species_new_group_topology_summary.tsv"
)

OUT_ALL_GLOB = Path(
    "work/codon/exact_species_full_all_glob_groups.txt"
)

OUT_MATCHES = Path(
    "work/codon/exact_species_final_one_to_one_min_cost_matches.tsv"
)

# ============================================================
# 1. Map new protein IDs -> orthogroups
# ============================================================

id_to_group = {}

with NEW_MEMBERS.open() as f:
    next(f)

    for line in f:
        group, full_id = line.rstrip("\n").split("\t")
        id_to_group[full_id] = group

print("New candidate protein members:", len(id_to_group))

# ============================================================
# 2. Read new DeepTMHMM classifications
# ============================================================

pred_class = {}

paths = sorted(
    glob.glob(
        "work/topology/deeptmhmm/"
        "exact_species_new_batches/"
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

            pred_class[full_id] = cls.strip()

print("New DeepTMHMM predictions:", len(pred_class))

# ============================================================
# 3. Determine new all-GLOB orthogroups
# ============================================================

group_classes = defaultdict(list)

for full_id, group in id_to_group.items():
    cls = pred_class.get(full_id, "MISSING")
    group_classes[group].append(cls)

new_all_glob = set()

with OUT_NEW_SUMMARY.open("w") as out:

    out.write(
        "orthogroup\tmembers\tGLOB\t"
        "non_GLOB\tall_GLOB\n"
    )

    for group in sorted(group_classes):

        classes = group_classes[group]

        n = len(classes)
        n_glob = sum(x == "GLOB" for x in classes)
        n_non = n - n_glob

        all_glob = (
            n > 0
            and n_glob == n
        )

        if all_glob:
            new_all_glob.add(group)

        out.write(
            f"{group}\t"
            f"{n}\t"
            f"{n_glob}\t"
            f"{n_non}\t"
            f"{all_glob}\n"
        )

print("New orthogroups evaluated:", len(group_classes))
print("New all-GLOB orthogroups:", len(new_all_glob))

# ============================================================
# 4. Combine with old all-GLOB set
# ============================================================

with OLD_ALL_GLOB.open() as f:
    old_all_glob = {
        x.strip()
        for x in f
        if x.strip()
    }

all_glob = old_all_glob | new_all_glob

with OUT_ALL_GLOB.open("w") as out:
    for g in sorted(all_glob):
        out.write(g + "\n")

print("Old all-GLOB:", len(old_all_glob))
print("Combined all-GLOB universe:", len(all_glob))

# ============================================================
# 5. Keep only exact-species ±5% all-GLOB candidate edges
# ============================================================

rows = []
membranes_with_any_edge = set()

with CANDIDATES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        soluble = row["candidate_orthogroup"]

        if soluble not in all_glob:
            continue

        r = {
            "membrane_orthogroup":
                row["membrane_orthogroup"],

            "species_coverage":
                int(row["species_coverage"]),

            "membrane_median_length":
                float(row["membrane_median_length"]),

            "candidate_rank":
                int(row["candidate_rank"]),

            "candidate_orthogroup":
                soluble,

            "candidate_median_length":
                float(row["candidate_median_length"]),

            "absolute_length_difference":
                float(row["absolute_length_difference"]),

            "relative_length_difference":
                float(row["relative_length_difference"]),
        }

        rows.append(r)
        membranes_with_any_edge.add(
            r["membrane_orthogroup"]
        )

print(
    "Exact-species membrane groups with >=1 all-GLOB candidate:",
    len(membranes_with_any_edge)
)

print(
    "Eligible all-GLOB candidate edges:",
    len(rows)
)

# ============================================================
# 6. Maximum-cardinality + minimum-cost bipartite matching
# ============================================================

membranes = sorted(
    {r["membrane_orthogroup"] for r in rows}
)

solubles = sorted(
    {r["candidate_orthogroup"] for r in rows}
)

M = len(membranes)
S = len(solubles)

source = 0
m_offset = 1
s_offset = 1 + M
sink = 1 + M + S
N = sink + 1

m_node = {
    g: m_offset + i
    for i, g in enumerate(membranes)
}

s_node = {
    g: s_offset + i
    for i, g in enumerate(solubles)
}

graph = [[] for _ in range(N)]

def add_edge(u, v, cap, cost, payload=None):
    graph[u].append(
        [v, cap, cost, len(graph[v]), payload]
    )

    graph[v].append(
        [u, 0, -cost, len(graph[u]) - 1, None]
    )

for m in membranes:
    add_edge(
        source,
        m_node[m],
        1,
        0
    )

for s in solubles:
    add_edge(
        s_node[s],
        sink,
        1,
        0
    )

for r in rows:

    # Primary objective:
    # minimize relative length difference.
    #
    # Rank only acts as deterministic tie-breaker.
    base = round(
        r["relative_length_difference"]
        * 1_000_000_000
    )

    cost = (
        base * 100
        + r["candidate_rank"]
    )

    add_edge(
        m_node[r["membrane_orthogroup"]],
        s_node[r["candidate_orthogroup"]],
        1,
        cost,
        payload=r
    )

flow = 0

while True:

    INF = 10**30

    dist = [INF] * N
    parent = [None] * N
    in_queue = [False] * N

    dist[source] = 0

    q = deque([source])
    in_queue[source] = True

    while q:

        u = q.popleft()
        in_queue[u] = False

        for ei, e in enumerate(graph[u]):

            v, cap, cost, rev, payload = e

            if cap <= 0:
                continue

            nd = dist[u] + cost

            if nd < dist[v]:

                dist[v] = nd
                parent[v] = (u, ei)

                if not in_queue[v]:
                    q.append(v)
                    in_queue[v] = True

    if parent[sink] is None:
        break

    v = sink

    while v != source:

        u, ei = parent[v]
        e = graph[u][ei]

        e[1] -= 1

        rev = e[3]

        graph[v][rev][1] += 1

        v = u

    flow += 1

# ============================================================
# 7. Recover final pairs
# ============================================================

final_rows = []

for m in membranes:

    u = m_node[m]

    for e in graph[u]:

        v, cap, cost, rev, payload = e

        if payload is not None and cap == 0:
            final_rows.append(payload)

final_rows.sort(
    key=lambda r:
        r["membrane_orthogroup"]
)

print()
print(
    "FINAL maximum-cardinality minimum-cost pairs:",
    len(final_rows)
)

print(
    "Unique membrane groups:",
    len({
        r["membrane_orthogroup"]
        for r in final_rows
    })
)

print(
    "Unique soluble groups:",
    len({
        r["candidate_orthogroup"]
        for r in final_rows
    })
)

# ============================================================
# 8. QC matching quality
# ============================================================

diffs = sorted(
    r["relative_length_difference"]
    for r in final_rows
)

ranks = defaultdict(int)

for r in final_rows:
    ranks[r["candidate_rank"]] += 1

print("Candidate-rank distribution:")

for rank in sorted(ranks):
    print(
        f"  rank {rank}: {ranks[rank]}"
    )

if diffs:

    n = len(diffs)

    if n % 2:
        median = diffs[n // 2]
    else:
        median = (
            diffs[n // 2 - 1]
            + diffs[n // 2]
        ) / 2

    print(
        "Median relative length difference:",
        median
    )

    print(
        "Mean relative length difference:",
        sum(diffs) / n
    )

    print(
        "Maximum relative length difference:",
        max(diffs)
    )

# ============================================================
# 9. Write final table
# ============================================================

fields = [
    "membrane_orthogroup",
    "species_coverage",
    "membrane_median_length",
    "candidate_rank",
    "candidate_orthogroup",
    "candidate_median_length",
    "absolute_length_difference",
    "relative_length_difference",
]

with OUT_MATCHES.open("w") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for r in final_rows:

        writer.writerow(
            {
                field: r[field]
                for field in fields
            }
        )

print()
print("Final output:", OUT_MATCHES)