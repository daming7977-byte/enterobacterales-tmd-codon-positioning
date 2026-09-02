from pathlib import Path
from collections import defaultdict, Counter
import re

OUT = Path("work/codon/species_codon_weights.tsv")

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",

    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",

    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",

    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def read_fasta(path):
    records = []
    seq = []

    with path.open() as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if seq:
                    records.append("".join(seq).upper())
                seq = []
            else:
                seq.append(line)

        if seq:
            records.append("".join(seq).upper())

    return records


cds_files = list(Path("data").rglob("cds_from_genomic.fna"))

print("CDS files:", len(cds_files))

rows = []

for path in sorted(cds_files):

    m = re.search(r"(GCF_\d+\.\d+)", str(path))
    if not m:
        continue

    genome = m.group(1)

    codon_counts = Counter()

    for seq in read_fasta(path):

        if len(seq) < 6:
            continue

        # Exclude initiation codon.
        start = 3

        for i in range(start, len(seq) - 2, 3):

            codon = seq[i:i+3]

            aa = CODON_TABLE.get(codon)

            if aa is None or aa == "*":
                continue

            codon_counts[codon] += 1

    aa_max = defaultdict(int)

    for codon, count in codon_counts.items():
        aa = CODON_TABLE[codon]
        aa_max[aa] = max(aa_max[aa], count)

    for codon in sorted(codon_counts):

        aa = CODON_TABLE[codon]
        count = codon_counts[codon]

        weight = count / aa_max[aa]

        rows.append(
            (genome, aa, codon, count, weight)
        )

    print(genome, "sense codons counted:", sum(codon_counts.values()))


with OUT.open("w") as out:

    out.write("genome\tamino_acid\tcodon\tcount\trelative_weight\n")

    for genome, aa, codon, count, weight in rows:

        out.write(
            f"{genome}\t{aa}\t{codon}\t{count}\t{weight:.6f}\n"
        )


print("Rows written:", len(rows))
print("Output:", OUT)