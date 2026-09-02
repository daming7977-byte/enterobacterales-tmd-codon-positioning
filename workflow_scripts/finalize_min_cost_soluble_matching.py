from pathlib import Path
from collections import defaultdict, deque
import csv

MATCHES = Path("work/codon/soluble_control_top5_matches.tsv")
ALL_GLOB = Path("work/codon/soluble_control_all_glob_groups.txt")
OUT = Path("work/codon/strict_soluble_one_to_one_min_cost_matches.tsv")

with ALL_GLOB.open() as f:
    all_glob = {x.strip() for x in f if x.strip()}

rows = []
membranes = set()
solubles = set()

with MATCHES.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["candidate_orthogroup"] not in all_glob:
            continue

        r = {
            "membrane_orthogroup": row["membrane_orthogroup"],
            "species_coverage": int(row["species_coverage"]),
            "membrane_median_length": float(row["membrane_median_length"]),
            "candidate_rank": int(row["candidate_rank"]),
            "candidate_orthogroup": row["candidate_orthogroup"],
            "candidate_median_length": float(row["candidate_median_length"]),
            "absolute_length_difference": float(row["absolute_length_difference"]),
            "relative_length_difference": float(row["relative_length_difference"]),
        }

        rows.append(r)
        membranes.add(r["membrane_orthogroup"])
        solubles.add(r["candidate_orthogroup"])

membranes = sorted(membranes)
solubles = sorted(solubles)

M = len(membranes)
S = len(solubles)

source = 0
m_offset = 1
s_offset = 1 + M
sink = 1 + M + S
N = sink + 1

membrane_node = {
    g: m_offset + i
    for i, g in enumerate(membranes)
}

soluble_node = {
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
    add_edge(source, membrane_node[m], 1, 0)

for s in solubles:
    add_edge(soluble_node[s], sink, 1, 0)

# Integer cost.
# Primary objective: relative length difference.
# Rank and IDs only break near-ties deterministically.
for r in rows:
    base = round(r["relative_length_difference"] * 1_000_000_000)
    rank_penalty = r["candidate_rank"]
    cost = base * 10 + rank_penalty

    add_edge(
        membrane_node[r["membrane_orthogroup"]],
        soluble_node[r["candidate_orthogroup"]],
        1,
        cost,
        payload=r,
    )

flow = 0
total_cost = 0

while True:
    INF = 10**30
    dist = [INF] * N
    inq = [False] * N
    parent = [None] * N

    dist[source] = 0
    q = deque([source])
    inq[source] = True

    while q:
        u = q.popleft()
        inq[u] = False

        for ei, e in enumerate(graph[u]):
            v, cap, cost, rev, payload = e

            if cap <= 0:
                continue

            nd = dist[u] + cost

            if nd < dist[v]:
                dist[v] = nd
                parent[v] = (u, ei)

                if not inq[v]:
                    q.append(v)
                    inq[v] = True

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
    total_cost += dist[sink]

final_rows = []

for m in membranes:
    u = membrane_node[m]

    for e in graph[u]:
        v, cap, cost, rev, payload = e

        if payload is not None and cap == 0:
            final_rows.append(payload)

final_rows.sort(
    key=lambda r: r["membrane_orthogroup"]
)

print("Maximum-cardinality minimum-cost matches:", len(final_rows))
print(
    "Unique membrane groups:",
    len({r["membrane_orthogroup"] for r in final_rows})
)
print(
    "Unique soluble groups:",
    len({r["candidate_orthogroup"] for r in final_rows})
)

unmatched = sorted(
    set(membranes)
    - {r["membrane_orthogroup"] for r in final_rows}
)

print("Unmatched among candidate-bearing groups:", len(unmatched))

for x in unmatched:
    print(" ", x)

rank_counts = defaultdict(int)

for r in final_rows:
    rank_counts[r["candidate_rank"]] += 1

print("Candidate rank distribution:")

for k in sorted(rank_counts):
    print(f"  rank {k}: {rank_counts[k]}")

diffs = sorted(
    r["relative_length_difference"]
    for r in final_rows
)

n = len(diffs)

if n % 2:
    median = diffs[n // 2]
else:
    median = (diffs[n // 2 - 1] + diffs[n // 2]) / 2

print("Median relative length difference:", median)
print("Mean relative length difference:", sum(diffs) / n)
print("Maximum relative length difference:", max(diffs))
print("Total relative length difference:", sum(diffs))

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

with OUT.open("w") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()

    for r in final_rows:
        writer.writerow({k: r[k] for k in fields})

print("Output:", OUT)