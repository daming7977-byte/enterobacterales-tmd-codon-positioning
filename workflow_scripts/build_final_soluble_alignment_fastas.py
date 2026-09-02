from pathlib import Path
from collections import defaultdict

MEMBERS = Path(
    "work/codon/final_soluble_matched_members.tsv"
)

ALL_FASTA = Path(
    "work/codon/final_soluble_matched_all.faa"
)

OUTDIR = Path(
    "work/codon/soluble_alignments/input"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    seqs = {}
    current = None
    parts = []

    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):

                if current is not None:
                    seqs[current] = "".join(parts)

                current = line[1:].split()[0]
                parts = []

            else:
                parts.append(line.strip())

        if current is not None:
            seqs[current] = "".join(parts)

    return seqs


seqs = read_fasta(ALL_FASTA)

groups = defaultdict(list)

with MEMBERS.open() as f:
    next(f)

    for line in f:
        group, full_id = line.rstrip("\n").split("\t")
        groups[group].append(full_id)


written_sequences = 0

for group, ids in sorted(groups.items()):

    outpath = OUTDIR / f"{group}.faa"

    with outpath.open("w") as out:

        for full_id in sorted(ids):

            seq = seqs[full_id]

            out.write(f">{full_id}\n")

            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

            written_sequences += 1


print("Orthogroups:", len(groups))
print("Sequences written:", written_sequences)
print("Output directory:", OUTDIR)