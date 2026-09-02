from pathlib import Path
from collections import defaultdict

genome_fa = Path(
    "data/unpacked/GCF_000006945.1/ncbi_dataset/data/"
    "GCF_000006945.1/GCF_000006945.1_ASM694v1_genomic.fna"
)

gff_file = Path(
    "data/unpacked/GCF_000006945.1/ncbi_dataset/data/"
    "GCF_000006945.1/genomic.gff"
)

out_file = Path(
    "data/unpacked/GCF_000006945.1/ncbi_dataset/data/"
    "GCF_000006945.1/cds_from_genomic.fna"
)


def read_fasta(path):
    seqs = {}
    current = None
    parts = []

    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(parts).upper()

                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())

        if current is not None:
            seqs[current] = "".join(parts).upper()

    return seqs


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


genome = read_fasta(genome_fa)

features = defaultdict(list)

with gff_file.open() as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attrs = fields

        if feature != "CDS":
            continue

        attr_dict = {}

        for item in attrs.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                attr_dict[k] = v

        protein_id = attr_dict.get("protein_id")

        if not protein_id:
            continue

        features[protein_id].append(
            (
                seqid,
                int(start),
                int(end),
                strand,
                int(phase) if phase != "." else 0,
            )
        )


written = 0

with out_file.open("w") as out:
    for protein_id, parts in sorted(features.items()):

        strands = {x[3] for x in parts}

        if len(strands) != 1:
            raise RuntimeError(f"Mixed strands for {protein_id}")

        strand = parts[0][3]

        if strand == "+":
            parts = sorted(parts, key=lambda x: x[1])
        else:
            parts = sorted(parts, key=lambda x: x[1], reverse=True)

        cds_parts = []

        for seqid, start, end, strand, phase in parts:
            seq = genome[seqid][start - 1:end]

            if strand == "-":
                seq = revcomp(seq)

            if phase:
                seq = seq[phase:]

            cds_parts.append(seq)

        cds = "".join(cds_parts)

        out.write(f">{protein_id}\n")

        for i in range(0, len(cds), 60):
            out.write(cds[i:i + 60] + "\n")

        written += 1


print("CDS written:", written)
print("Output:", out_file)