from pathlib import Path
from collections import defaultdict

infile = Path("work/codon/codon_position_adaptation.tsv")
outfile = Path("work/codon/low_adaptation_segments.tsv")

proteins = defaultdict(list)

with infile.open() as f:
    header = next(f)

    for line in f:
        fields = line.rstrip("\n").split("\t")

        group = fields[0]
        genome = fields[1]
        protein_id = fields[2]
        codon_position = int(fields[3])
        codon = fields[4]
        aa = fields[5]
        weight = float(fields[6])
        low = fields[7] == "True"

        proteins[(group, genome, protein_id)].append(
            (codon_position, codon, aa, weight, low)
        )


segments = []

for (group, genome, protein_id), records in proteins.items():

    records.sort()

    current = []

    def save_segment(seg):
        if len(seg) < 3:
            return

        start = seg[0][0]
        end = seg[-1][0]
        length = len(seg)

        mean_weight = sum(x[3] for x in seg) / length

        segments.append(
            (
                group,
                genome,
                protein_id,
                start,
                end,
                length,
                mean_weight,
            )
        )

    previous_pos = None

    for record in records:

        pos = record[0]
        low = record[4]

        if low:

            if (
                current
                and previous_pos is not None
                and pos != previous_pos + 1
            ):
                save_segment(current)
                current = []

            current.append(record)

        else:

            if current:
                save_segment(current)
                current = []

        previous_pos = pos

    if current:
        save_segment(current)


segments.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[2],
        x[3],
    )
)


with outfile.open("w") as out:

    out.write(
        "orthogroup\tgenome\tprotein_id\t"
        "start_codon\tend_codon\tlength_codons\tmean_relative_weight\n"
    )

    for row in segments:

        out.write(
            f"{row[0]}\t"
            f"{row[1]}\t"
            f"{row[2]}\t"
            f"{row[3]}\t"
            f"{row[4]}\t"
            f"{row[5]}\t"
            f"{row[6]:.6f}\n"
        )


print("Proteins evaluated:", len(proteins))
print("Segments >=3 codons:", len(segments))

proteins_with_segments = len(
    set((x[0], x[1], x[2]) for x in segments)
)

print("Proteins with >=1 segment:", proteins_with_segments)
print("Output:", outfile)